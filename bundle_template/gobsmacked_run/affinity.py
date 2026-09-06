"""Stage 5: what would this thing bind at?

Docking scores rank; they do not predict. PandaDock's number orders ten poses of
one ligand against one another and means nothing across targets. This stage
answers a different question, on the one structure in the pipeline that has been
through explicit solvent at 300 K.

The pattern is Boltzina's (Furui and Ohue 2025). Boltz-2 predicts affinity in
three parts: a trunk that builds pairwise representations from the sequence and
an MSA, a structure module that diffuses coordinates from them, and an affinity
head that reads both. The structure module is the expensive part and it is
redundant when a pose already exists, so an existing pose goes straight to the
affinity head. They used rigid Vina poses. Ours have been docked, merged,
solvated, minimised, equilibrated and run for a nanosecond, which is a better
pose fed to the same head.

Two things are scored, always: the docked pose before MD and the relaxed complex
after it. The pair is the point, and it is nearly free, because the trunk and
the MSA are shared and a second pose costs one more pass of the affinity head.
It answers whether relaxation improved the complex or merely moved it, which
neither number says alone.

Nothing here enters the GOBSMACK score. Every graded metric on that card has a
crystal to be right or wrong about and this has none.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Optional

from .console import bar_for

# Shared with the GNN checkpoint's cache, and for the same reason: an MSA and a
# 1 GB model should be fetched once per machine, not once per bundle.
CACHE = Path.home() / ".gobsmacked"
MSA_CACHE = CACHE / "msa"
LIGAND_ID = "LIG"
DEFAULT_N_FRAMES = 5
DEFAULT_WINDOW = 0.2
DEFAULT_RECYCLING = 3
# How far the prediction may deviate from the template, in Angstroms. The point
# of the template is to hold the protein at the conformation MD produced, so it
# is tight: 2 A is the scale the scorecard already calls a good pocket backbone
# agreement. Loose enough and the template stops meaning anything; tighter and
# the potential has nothing to relax into.
DEFAULT_THRESHOLD_A = 2.0


def run(campaign: dict, work: Path, results: Path, log) -> dict[str, Any]:
    cfg = campaign.get("affinity") or {}
    out_dir = results / "affinity"
    warnings: list[str] = []

    if not cfg.get("include", True):
        write(out_dir, {"schema": "1.0", "requested": False, "ran": False,
                        "reason": "the campaign did not ask for it"})
        log("affinity: not requested by the campaign")
        return {"warnings": warnings, "headline": "not requested"}

    started = time.time()
    sequence = ((campaign.get("protein") or {}).get("sequence") or "").strip()
    smiles = ((campaign.get("ligand") or {}).get("smiles") or "").strip()
    if not sequence or not smiles:
        return decline(out_dir, log, warnings,
                       "the campaign carries no sequence or no SMILES")

    if shutil.which("pixi") is None and os.environ.get("GOBSMACKED_BOLTZ") is None:
        return decline(out_dir, log, warnings,
                       "pixi is not on PATH, so the affinity environment cannot be entered")

    poses = collect_poses(results, cfg, log)
    if not poses:
        return decline(out_dir, log, warnings, "no pose to score")

    frames_dir = out_dir / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    try:
        written = [(name, to_cif(path, frames_dir / f"{name}.cif")) for name, path in poses]
    except Exception as exc:                       # noqa: BLE001 - decorative stage, never fatal
        return decline(out_dir, log, warnings, f"the poses could not be converted to mmCIF: {exc}")

    # From the first pose's own coordinates, not from the campaign: see
    # sequence_of. Every pose in one run is the same chain, so one lookup does.
    scored_sequence = sequence_of(poses[0][1]) or sequence
    if scored_sequence != sequence:
        log(f"affinity: scoring the {len(scored_sequence)} residues in the structure, "
            f"not the {len(sequence)} in the campaign")

    msa = ensure_msa(scored_sequence, log, warnings)
    scored, error = predict(written, scored_sequence, smiles, msa, cfg, out_dir, log)
    if error:
        return decline(out_dir, log, warnings, error)

    block = summarise(scored, cfg, msa, poses, time.time() - started)
    write(out_dir, block)
    headline = describe(block)
    log(f"affinity: {headline}")
    return {"warnings": warnings, "headline": headline,
            "pic50_post": (block.get("post_md") or {}).get("pic50_mean")}


# ---------------------------------------------------------------------------
# Poses
# ---------------------------------------------------------------------------

def collect_poses(results: Path, cfg: dict, log) -> list[tuple[str, Path]]:
    """The pre-MD pose, then the post-MD frame or frames, all as PDB on disk.

    The pre-MD pose is scored exactly once: it is a single structure and there
    is nothing to sample. The mode governs the post-MD side alone. Trajectory
    frames are written out here rather than passed along as an index, because
    everything downstream reads a file and a DCD frame is not one.
    """
    poses: list[tuple[str, Path]] = []
    pose1 = results / "complex_pose1.pdb"
    if pose1.exists():
        poses.append(("pre_md", pose1))

    final = results / "complex_md_final.pdb"
    mode = (cfg.get("frames") or "cluster").lower()
    dcd = results / "traj" / "traj.dcd"
    topology = results / "traj" / "topology.pdb"

    if mode == "single" or not dcd.exists() or not topology.exists():
        if final.exists():
            poses.append(("post_md_final", final))
        return poses

    extracted = extract_frames(topology, dcd, cfg, results / "affinity" / "frames")
    if not extracted:
        if final.exists():
            poses.append(("post_md_final", final))
        return poses
    poses.extend(extracted)
    log(f"affinity: the docked pose and {len(extracted)} frames from the last "
        f"{float(cfg.get('window_fraction', DEFAULT_WINDOW)) * 100:.0f}% of the run")
    return poses


def extract_frames(topology: Path, dcd: Path, cfg: dict,
                   dest_dir: Path) -> list[tuple[str, Path]]:
    """Write the chosen trajectory frames out as PDB.

    The trajectory holds the solute only and its ligand is already named LIG, so
    a frame is a complete complex and needs no reassembly. MDTraj announces
    itself from C on every DCD open, which is hushed for the same reason it is
    hushed in summarise.
    """
    import mdtraj as md

    from .noise import hush_c_stdout

    dest_dir.mkdir(parents=True, exist_ok=True)
    out: list[tuple[str, Path]] = []
    with hush_c_stdout():
        traj = md.load(str(dcd), top=str(topology))
        # The trajectory's own frame count, read here rather than from a summary
        # file that this stage runs before.
        for index in frame_indices(traj.n_frames, cfg):
            if index >= traj.n_frames:
                continue
            path = dest_dir / f"post_md_{index:04d}.pdb"
            traj[index].save_pdb(str(path))
            out.append((f"post_md_{index:04d}", path))
    return out


def frame_indices(n_frames: int, cfg: dict) -> list[int]:
    """Evenly spaced frames across the tail of the production run.

    Takes a count rather than a directory, because the count used to come from
    `traj/summary.json` and that file does not exist yet: summarise runs AFTER
    this stage. The lookup failed, the cluster silently fell back to the single
    final frame, and a run that asked for five frames reported one with a spread
    of zero. Nothing errored; the number was simply less than it claimed.

    Evenly spaced rather than clustered by RMSD: the question is whether the
    prediction is steady over the end of the run, and a clustering that picks
    the most distinct frames answers a different one.
    """
    if n_frames < 2:
        return []
    wanted = max(1, int(cfg.get("n_frames", DEFAULT_N_FRAMES) or DEFAULT_N_FRAMES))
    window = float(cfg.get("window_fraction", DEFAULT_WINDOW) or DEFAULT_WINDOW)
    first = max(0, int(n_frames * (1.0 - window)))
    span = n_frames - first
    if span <= wanted:
        return list(range(first, n_frames))
    step = span / wanted
    return sorted({min(n_frames - 1, int(first + i * step)) for i in range(wanted)})


def to_cif(source: Path, dest: Path) -> Path:
    """The pose's protein, as mmCIF, with the ligand left out.

    Two things had to be got right here and neither is guessable.

    The ligand is removed. Boltz's template reader parses templates as polymers,
    and a ligand chain in one walks its residue index off the end of the
    sequence: `res_name = sequence[j]`, IndexError, inside parse_polymer. The
    ligand is declared separately in the YAML as SMILES, which is where Boltz
    wants it.

    mmCIF rather than PDB because that is what Boltz's template reader takes.
    """
    import gemmi

    structure = gemmi.read_structure(str(source))
    # setup_entities() first: remove_ligands_and_waters() reads the entity_type
    # that only setup assigns, and raises "missing entity_type in chain A"
    # without it. Called again afterwards to renumber what is left.
    structure.setup_entities()
    structure.remove_ligands_and_waters()
    structure.remove_empty_chains()
    structure.setup_entities()
    fill_entity_sequences(structure)
    structure.assign_label_seq_id()
    dest.parent.mkdir(parents=True, exist_ok=True)
    structure.make_mmcif_document().write_file(str(dest))
    return dest


def fill_entity_sequences(structure) -> None:
    """Give each polymer entity the canonical sequence gemmi leaves empty.

    Without this the mmCIF carries no `_entity_poly_seq` loop, because gemmi
    only writes one for an entity whose full_sequence is set and `setup_entities`
    does not set it from the coordinates. Boltz builds its residue list from
    exactly that loop, so a template without it gives an empty list and the
    first residue lookup raises `IndexError: list index out of range` inside
    parse_polymer. The template is silently useless and the error names neither
    the file nor the reason.
    """
    import gemmi

    for entity in structure.entities:
        if entity.entity_type != gemmi.EntityType.Polymer or entity.full_sequence:
            continue
        for model in structure:
            for chain in model:
                polymer = chain.get_polymer()
                if len(polymer):
                    entity.full_sequence = [residue.name for residue in polymer]
                    break
            if entity.full_sequence:
                break


def sequence_of(path: Path) -> str:
    """The one-letter sequence of the structure being scored.

    Taken from the structure and never from the campaign. The campaign carries
    UniProt's full precursor -- 1,210 residues for EGFR -- while the pose is the
    trimmed domain, 253 of them. Boltz maps a template onto the input sequence
    by index, so handing it the precursor beside a domain template walks off the
    end of the list and raises IndexError deep inside its parser, which is
    exactly what happened on the first real run.

    Deriving it from the file also makes the two agree by construction, rather
    than by two pieces of code happening to trim the same way.
    """
    import gemmi

    structure = gemmi.read_structure(str(path))
    structure.setup_entities()
    structure.remove_ligands_and_waters()
    for model in structure:
        for chain in model:
            polymer = chain.get_polymer()
            if len(polymer):
                return gemmi.one_letter_code(polymer.extract_sequence())
    return ""


# ---------------------------------------------------------------------------
# MSA
# ---------------------------------------------------------------------------

def msa_key(sequence: str) -> str:
    """The cache key is the PREPARED sequence, not the accession.

    A campaign may trim to a domain, and the alignment for residues 714-966 of
    EGFR is not the alignment for the 1,210-residue precursor. Keying on the
    accession would hand the second to a campaign that asked for the first.
    """
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()[:16]


def ensure_msa(sequence: str, log, warnings: list[str]) -> dict[str, Any]:
    """Look in the cache before asking the network. Say which one answered."""
    MSA_CACHE.mkdir(parents=True, exist_ok=True)
    key = msa_key(sequence)
    path = MSA_CACHE / f"{key}.a3m"
    for path in (MSA_CACHE / f"{key}.csv", MSA_CACHE / f"{key}.a3m"):
        if path.exists() and path.stat().st_size > 0:
            age_days = (time.time() - path.stat().st_mtime) / 86400
            log(f"affinity: MSA cache hit for {key}, written {age_days:.0f} days ago")
            return {"cached": True, "key": key, "path": str(path),
                    "depth": depth_of(path), "seconds": 0.0,
                    "age_days": round(age_days, 1)}
    return {"cached": False, "key": key, "path": None, "depth": None, "seconds": None}


def capture_msa(work: Path, msa: dict, log) -> dict:
    """Keep the alignment the first pose paid for.

    Boltz writes the MSA it generated to `msa/<name>_0.csv`, in exactly the
    `key,sequence` form its own `msa:` field accepts. Without this every pose
    calls the public ColabFold server again for the same protein: five frames
    plus the docked pose is six identical queries, six times the wait, against a
    free service. Copied into the shared cache, the second pose onwards reads it
    from disk and every later run on the same trimmed sequence does too.
    """
    if msa.get("cached"):
        return msa
    candidates = sorted(work.rglob("msa/*_0.csv")) + sorted(work.rglob("msa/*_0.a3m"))
    for candidate in candidates:
        if not candidate.exists() or candidate.stat().st_size == 0:
            continue
        MSA_CACHE.mkdir(parents=True, exist_ok=True)
        dest = MSA_CACHE / f"{msa['key']}{candidate.suffix}"
        # Copied then moved, so a half-written cache entry is never found by a
        # run that starts while this one is still copying 3 MB.
        partial = dest.with_suffix(candidate.suffix + ".part")
        shutil.copyfile(candidate, partial)
        partial.replace(dest)
        depth = depth_of(dest)
        log(f"affinity: cached the MSA for {msa['key']}, {depth or '?'} sequences, "
            f"so the remaining poses do not ask the server again")
        return {**msa, "cached": True, "path": str(dest), "depth": depth,
                "generated_here": True}
    return msa


def depth_of(path: Path) -> Optional[int]:
    """Sequences in the alignment: header lines in a3m, data rows in CSV."""
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            if path.suffix == ".csv":
                return max(0, sum(1 for _ in fh) - 1)      # minus the header row
            return sum(1 for line in fh if line.startswith(">"))
    except OSError:
        return None


# ---------------------------------------------------------------------------
# Boltz
# ---------------------------------------------------------------------------

def predict(written: list[tuple[str, Path]], sequence: str, smiles: str,
            msa: dict, cfg: dict, out_dir: Path, log) -> tuple[list[dict], Optional[str]]:
    """Score every pose, one subprocess per pose, and never raise."""
    scored: list[dict] = []
    log_path = out_dir / "boltz.log"
    log_path.write_text("", encoding="utf-8")
    with bar_for(log, "scoring poses with the Boltz-2 affinity head",
                 total=len(written)) as bar:
        for index, (name, cif) in enumerate(written):
            bar.update(index, note=name.replace("_", " "))
            result, error = predict_one(name, cif, sequence, smiles, msa, cfg, out_dir, log_path)
            if error:
                return scored, error
            scored.append(result)
            # The first pose pays for the alignment; the rest read it off disk.
            msa = capture_msa(out_dir / "boltz" / name, msa, log)
        bar.update(len(written), note="done")
    return scored, None


def predict_one(name: str, cif: Path, sequence: str, smiles: str, msa: dict,
                cfg: dict, out_dir: Path, log_path: Path) -> tuple[dict, Optional[str]]:
    yaml_path = out_dir / "frames" / f"{name}.yaml"
    threshold = float(cfg.get("template_threshold_a", DEFAULT_THRESHOLD_A)
                      or DEFAULT_THRESHOLD_A)
    yaml_path.write_text(boltz_input(sequence, smiles, cif, msa, threshold), encoding="utf-8")
    work = out_dir / "boltz" / name
    work.mkdir(parents=True, exist_ok=True)

    cmd = boltz_command(yaml_path, work, cfg, msa)
    # encoding named, not inherited. text=True decodes the child's output with
    # the LOCALE's encoding, which inside a pixi task is ASCII, and boltz prints
    # UTF-8: a progress bar's box characters raised UnicodeDecodeError and took
    # the stage down after MD had already run. Same root cause as the log crash,
    # a different call.
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as fh:
        fh.write(f"$ {' '.join(cmd)}\n{proc.stdout}\n{proc.stderr}\n")
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-3:]
        return {}, "boltz failed: " + " / ".join(tail)

    values = read_prediction(work)
    if values is None:
        return {}, f"boltz wrote no affinity for {name}"
    values["name"] = name
    return values, None


def boltz_command(yaml_path: Path, work: Path, cfg: dict, msa: dict) -> list[str]:
    """The subprocess, entering the affinity environment.

    `--use_potentials` IS passed, and it has to be: Boltz enforces a forced
    template through a steering potential, so without the flag `force: true`
    is silently inert and the pose this pipeline spent an hour producing does
    not constrain anything.

    It is not free. On Apple silicon `--use_potentials` can kill the whole
    predict process with a shape mismatch, from a masked assignment on MPS
    returning different counts for the same mask. It is trajectory-dependent,
    so it looks intermittent and is not, and it is a crash rather than an
    exception: the subprocess dies and this stage records why. That is survivable
    here precisely because this stage is one subprocess per pose and never fatal
    to the run.
    """
    override = os.environ.get("GOBSMACKED_BOLTZ")
    prefix = override.split() if override else ["pixi", "run", "-e", "affinity", "boltz"]
    cmd = [*prefix, "predict", str(yaml_path),
           "--out_dir", str(work),
           "--recycling_steps", str(int(cfg.get("recycling_steps", DEFAULT_RECYCLING) or 3)),
           "--diffusion_samples", "1",
           "--output_format", "mmcif",
           "--use_potentials"]
    if not msa.get("cached"):
        cmd += ["--use_msa_server"]
    return cmd


def boltz_input(sequence: str, smiles: str, cif: Path, msa: dict,
                threshold: float = DEFAULT_THRESHOLD_A) -> str:
    """Boltz's YAML: the complex as a template, and the affinity head asked for."""
    lines = [
        "version: 1",
        "sequences:",
        "  - protein:",
        "      id: A",
        f"      sequence: {sequence}",
    ]
    if msa.get("path"):
        lines.append(f"      msa: {msa['path']}")
    lines += [
        "  - ligand:",
        f"      id: {LIGAND_ID}",
        f"      smiles: '{smiles}'",
        "templates:",
        f"  - cif: {cif}",
        # Forced, with a threshold in Angstroms: the prediction may not deviate
        # from this structure by more than that. The whole point of the stage is
        # that the affinity head reads the conformation MD produced rather than
        # one the model invented, and an unforced template only conditions the
        # trunk.
        #
        # Boltz implements this as a steering potential, so it does nothing
        # unless --use_potentials is passed. See boltz_command.
        "    force: true",
        f"    threshold: {threshold}",
        "properties:",
        "  - affinity:",
        f"      binder: {LIGAND_ID}",
    ]
    return "\n".join(lines) + "\n"


def read_prediction(work: Path) -> Optional[dict]:
    """Boltz writes affinity_*.json somewhere under its output directory."""
    for path in sorted(work.rglob("affinity*.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            continue
        if not isinstance(data, dict):
            continue
        value = data.get("affinity_pred_value")
        if value is None:
            continue
        return {
            "affinity_pred_value": round(float(value), 4),
            "pic50": pic50(float(value)),
            "affinity_probability_binary": (
                round(float(data["affinity_probability_binary"]), 4)
                if data.get("affinity_probability_binary") is not None else None),
        }
    return None


def pic50(value: float) -> float:
    """pIC50 from Boltz's log10(IC50 in micromolar).

    The model reports log10 of an IC50 in micromolar, where LOWER is stronger.
    pIC50 is the negative log10 of the same quantity in molar, where HIGHER is
    stronger. They differ by a sign as well as an offset, so a mix-up gives a
    plausible number rather than an obvious one:

        IC50 = 10^value micromolar = 10^(value - 6) molar
        pIC50 = -log10(10^(value - 6)) = 6 - value
    """
    return round(6.0 - value, 3)


# ---------------------------------------------------------------------------
# The block
# ---------------------------------------------------------------------------

def summarise(scored: list[dict], cfg: dict, msa: dict,
              poses: list[tuple[str, Path]], seconds: float) -> dict[str, Any]:
    import statistics

    pre = next((s for s in scored if s["name"] == "pre_md"), None)
    post = [s for s in scored if s["name"].startswith("post_md")]

    block: dict[str, Any] = {
        "schema": "1.0",
        "requested": True,
        "ran": True,
        # Not "boltzina": the diffusion module still runs. What it does not do
        # is wander, because the template is forced through a steering potential
        # and the prediction may not leave it by more than the threshold. Named
        # for what happened rather than for the pattern it approximates.
        "route": "boltz-forced-template",
        "engine": {"recycling_steps": int(cfg.get("recycling_steps", DEFAULT_RECYCLING) or 3),
                   "template_threshold_a": float(cfg.get("template_threshold_a",
                                                         DEFAULT_THRESHOLD_A)
                                                 or DEFAULT_THRESHOLD_A),
                   "potentials": True,
                   "diffusion": "one sample, constrained to the template"},
        "frames": {"mode": (cfg.get("frames") or "cluster").lower(),
                   "scored": [name for name, _ in poses]},
        "msa": msa,
        "unit": ("affinity_pred_value is log10(IC50) with IC50 in micromolar, lower is "
                 "stronger; pIC50 = 6 - affinity_pred_value, higher is stronger"),
        "seconds": round(seconds, 1),
        "warnings": [],
    }
    if pre:
        block["pre_md"] = {k: v for k, v in pre.items() if k != "name"}
    if post:
        values = [s["affinity_pred_value"] for s in post]
        pics = [s["pic50"] for s in post]
        probs = [s["affinity_probability_binary"] for s in post
                 if s["affinity_probability_binary"] is not None]
        block["post_md"] = {
            "per_frame": post,
            "affinity_pred_value_mean": round(statistics.fmean(values), 4),
            "affinity_pred_value_sd": round(statistics.pstdev(values), 4) if len(values) > 1 else 0.0,
            "affinity_pred_value_range": [min(values), max(values)],
            "pic50_mean": round(statistics.fmean(pics), 3),
            "pic50_sd": round(statistics.pstdev(pics), 3) if len(pics) > 1 else 0.0,
            "probability_binary_mean": round(statistics.fmean(probs), 4) if probs else None,
            "probability_binary_sd": (round(statistics.pstdev(probs), 4)
                                      if len(probs) > 1 else 0.0 if probs else None),
        }
    if pre and post:
        block["delta"] = {
            "affinity_pred_value": round(
                block["post_md"]["affinity_pred_value_mean"] - pre["affinity_pred_value"], 4),
            "pic50": round(block["post_md"]["pic50_mean"] - pre["pic50"], 3),
            "probability_binary": (
                round(block["post_md"]["probability_binary_mean"]
                      - pre["affinity_probability_binary"], 4)
                if block["post_md"]["probability_binary_mean"] is not None
                and pre["affinity_probability_binary"] is not None else None),
            "note": ("post-MD minus pre-MD; a positive pIC50 means MD improved the "
                     "predicted binding"),
        }
    return block


def describe(block: dict) -> str:
    post = block.get("post_md") or {}
    pre = block.get("pre_md") or {}
    if not post:
        return "scored the docked pose only"
    text = f"pIC50 {post['pic50_mean']}"
    if post.get("pic50_sd"):
        text += f" ± {post['pic50_sd']}"
    if pre:
        text += f" after MD, {pre['pic50']} before"
    return text


def decline(out_dir: Path, log, warnings: list[str], reason: str) -> dict[str, Any]:
    """Never fatal. A missing affinity is a missing panel, not a failed run."""
    write(out_dir, {"schema": "1.0", "requested": True, "ran": False, "reason": reason})
    warnings.append(f"No affinity prediction: {reason}.")
    log(f"affinity: skipped, {reason}")
    return {"warnings": warnings, "headline": f"skipped: {reason}"}


def write(out_dir: Path, block: dict) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "affinity.json").write_text(json.dumps(block, indent=2), encoding="utf-8")
