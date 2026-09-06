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

    msa = ensure_msa(sequence, log, warnings)
    scored, error = predict(written, sequence, smiles, msa, cfg, out_dir, log)
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

    indices = frame_indices(results, cfg)
    if not indices:
        if final.exists():
            poses.append(("post_md_final", final))
        return poses

    extracted = extract_frames(topology, dcd, indices, results / "affinity" / "frames")
    if not extracted:
        if final.exists():
            poses.append(("post_md_final", final))
        return poses
    poses.extend(extracted)
    log(f"affinity: the docked pose and {len(extracted)} frames from the last "
        f"{float(cfg.get('window_fraction', DEFAULT_WINDOW)) * 100:.0f}% of the run")
    return poses


def extract_frames(topology: Path, dcd: Path, indices: list[int],
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
        for index in indices:
            if index >= traj.n_frames:
                continue
            path = dest_dir / f"post_md_{index:04d}.pdb"
            traj[index].save_pdb(str(path))
            out.append((f"post_md_{index:04d}", path))
    return out


def frame_indices(results: Path, cfg: dict) -> list[int]:
    """Evenly spaced frames across the tail of the production run.

    Evenly spaced rather than clustered by RMSD: the question is whether the
    prediction is steady over the end of the run, and a clustering that picks
    the most distinct frames answers a different one.
    """
    summary_path = results / "traj" / "summary.json"
    if not summary_path.exists():
        return []
    try:
        frames = int(json.loads(summary_path.read_text(encoding="utf-8")).get("frames") or 0)
    except (json.JSONDecodeError, TypeError, ValueError, OSError):
        return []
    if frames < 2:
        return []
    wanted = max(1, int(cfg.get("n_frames", DEFAULT_N_FRAMES) or DEFAULT_N_FRAMES))
    window = float(cfg.get("window_fraction", DEFAULT_WINDOW) or DEFAULT_WINDOW)
    first = max(0, int(frames * (1.0 - window)))
    span = frames - first
    if span <= wanted:
        return list(range(first, frames))
    step = span / wanted
    return sorted({min(frames - 1, int(first + i * step)) for i in range(wanted)})


def to_cif(source: Path, dest: Path) -> Path:
    """PDB in, mmCIF out.

    mmCIF rather than PDB because the ligand is a full residue with bond orders
    that a HETATM block cannot carry, and because Boltz's own reader is happiest
    there.
    """
    import gemmi

    structure = gemmi.read_structure(str(source))
    structure.setup_entities()
    dest.parent.mkdir(parents=True, exist_ok=True)
    structure.make_mmcif_document().write_file(str(dest))
    return dest


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
    if path.exists() and path.stat().st_size > 0:
        age_days = (time.time() - path.stat().st_mtime) / 86400
        log(f"affinity: MSA cache hit for {key}, written {age_days:.0f} days ago")
        return {"cached": True, "key": key, "path": str(path),
                "depth": depth_of(path), "seconds": 0.0,
                "age_days": round(age_days, 1)}
    return {"cached": False, "key": key, "path": None, "depth": None, "seconds": None}


def depth_of(path: Path) -> Optional[int]:
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
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
        bar.update(len(written), note="done")
    return scored, None


def predict_one(name: str, cif: Path, sequence: str, smiles: str, msa: dict,
                cfg: dict, out_dir: Path, log_path: Path) -> tuple[dict, Optional[str]]:
    yaml_path = out_dir / "frames" / f"{name}.yaml"
    yaml_path.write_text(boltz_input(sequence, smiles, cif, msa), encoding="utf-8")
    work = out_dir / "boltz" / name
    work.mkdir(parents=True, exist_ok=True)

    cmd = boltz_command(yaml_path, work, cfg, msa)
    proc = subprocess.run(cmd, capture_output=True, text=True)
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

    `--use_potentials` is never passed. It steers diffusion, which is bypassed
    here, and on Apple silicon it kills the whole predict process with a shape
    mismatch that is trajectory-dependent, so it looks intermittent and is not.
    """
    override = os.environ.get("GOBSMACKED_BOLTZ")
    prefix = override.split() if override else ["pixi", "run", "-e", "affinity", "boltz"]
    cmd = [*prefix, "predict", str(yaml_path),
           "--out_dir", str(work),
           "--recycling_steps", str(int(cfg.get("recycling_steps", DEFAULT_RECYCLING) or 3)),
           "--diffusion_samples", "1",
           "--output_format", "mmcif"]
    if not msa.get("cached"):
        cmd += ["--use_msa_server"]
    return cmd


def boltz_input(sequence: str, smiles: str, cif: Path, msa: dict) -> str:
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
        "    force: true",
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
        "route": "boltz-template",
        "engine": {"recycling_steps": int(cfg.get("recycling_steps", DEFAULT_RECYCLING) or 3),
                   "diffusion": "one sample, template forced"},
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
