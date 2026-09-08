"""Stage 1: fold, co-fold, or skip.

Skipped whenever `model_apo.pdb` is already in the bundle, which is the usual
case: the server fetches an AlphaFold DB, ESM Atlas, PDB or user-supplied
structure and ships it. Folding here is the fallback for a sequence nothing has
a model for, and it is the reason this bundle wants a GPU at all.

Two ways to fold, chosen by `fold.method` in the campaign:

`esmfold`   sequence in, apo structure out. Fast, and knows nothing about the
            ligand that is about to be docked into it.

`boltz2`    the protein and the ligand folded TOGETHER, then the ligand thrown
            away and the protein kept. Measured across five starting structures
            on EGFR plus erlotinib, this was the only one whose docked pose came
            back at rank 1 (1.40 A against 4.29 to 8.77 A for a crystal, an
            AlphaFold model and an ESMFold model), and the only run to grade B.

            The reason it is the protein and not the pose that is kept: in the
            same experiment Boltz-2's own ligand placement was 4.5 to 4.8 A from
            the crystal while the pocket it built was the closest of any starting
            structure. Its protein is better than its ligand. So the pocket is
            taken and the ligand is re-docked into it. The co-folded pose is
            still written out, as `cofold_ligand.sdf`, because a reader should be
            able to see what was discarded.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any, Optional

from . import console
from .console import bar_for

# ESMFold's memory use grows with sequence length; chunking trades speed for
# peak memory and is the difference between running and an out-of-memory abort
# on a 12 GB card. These are the thresholds the ESM authors recommend.
CHUNK_SIZES = ((700, None), (1000, 128), (1600, 64), (10_000, 32))

POCKET_PLDDT_WARNING = 70.0


def chunk_size_for(length: int) -> Optional[int]:
    for limit, chunk in CHUNK_SIZES:
        if length <= limit:
            return chunk
    return 32


def method_of(campaign: dict) -> str:
    """`esmfold` unless the campaign asked for co-folding. Unknown values fall
    back rather than raising: a bundle that reaches a machine running an older
    runner should still produce a structure."""
    wanted = str(((campaign.get("fold") or {}).get("method") or "esmfold")).lower()
    return wanted if wanted in ("esmfold", "boltz2") else "esmfold"


def run(campaign: dict, work: Path, results: Path, log) -> dict[str, Any]:
    """Write `model_apo.pdb` and, when folded, `plddt.json`."""
    protein = campaign.get("protein") or {}
    sequence = (protein.get("sequence") or "").strip().upper()
    supplied = work.parent / "model_apo.pdb"
    target = results / "model_apo.pdb"
    warnings: list[str] = []

    if supplied.exists():
        target.write_text(supplied.read_text(encoding="utf-8"), encoding="utf-8")
        log(f"fold: skipped, using the supplied {protein.get('source_structure')} model "
            f"{protein.get('source_id') or ''}".rstrip())
        return {"folded": False, "warnings": warnings,
                "headline": "the supplied model was used"}

    if not sequence:
        raise RuntimeError("No structure in the bundle and no sequence to fold.")

    if method_of(campaign) == "boltz2":
        return cofold(campaign, sequence, work, results, target, log, warnings)

    return _esmfold(campaign, sequence, results, target, log, warnings)


def per_residue_plddt(pdb_text: str) -> dict[str, float]:
    """ESMFold writes pLDDT into the B-factor column; CA carries the residue's."""
    out: dict[str, float] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        # Fixed-width PDB: columns are the format, not whitespace. Splitting on
        # whitespace breaks the moment a coordinate runs into its neighbour.
        if line[12:16].strip() != "CA":
            continue
        try:
            out[str(int(line[22:26]))] = float(line[60:66])
        except ValueError:
            continue
    return out


def pocket_residue_numbers(campaign: dict) -> list[int]:
    residues = ((campaign.get("pocket") or {}).get("residues") or [])
    numbers = []
    for item in residues:
        text = str(item)
        tail = text.split(":")[-1]
        try:
            numbers.append(int(tail))
        except ValueError:
            continue
    return numbers


# ---------------------------------------------------------------------------
# Co-folding with Boltz-2
# ---------------------------------------------------------------------------

def cofold(campaign: dict, sequence: str, work: Path, results: Path,
           target: Path, log, warnings: list[str]) -> dict[str, Any]:
    """Fold the protein WITH the ligand, then keep only the protein.

    Everything Boltz needs is already in this package for the affinity stage:
    the subprocess into the `affinity` environment, the MSA cache keyed on the
    sequence, and the YAML shape. Reused rather than reimplemented, so the MSA
    computed here is the one the affinity stage will find later and neither asks
    the public server twice for the same protein.

    Failure is not fatal. A co-fold that does not produce a structure falls back
    to ESMFold with the reason recorded, because the alternative is a bundle that
    has spent an hour and produced nothing to dock.
    """
    from . import affinity as aff

    out_dir = results / "cofold"
    out_dir.mkdir(parents=True, exist_ok=True)
    smiles = ((campaign.get("ligand") or {}).get("smiles") or "").strip()
    if not smiles:
        warnings.append("Co-folding was requested but the campaign carries no SMILES; "
                        "ESMFold was used instead.")
        log("fold: co-fold asked for but there is no ligand, falling back to ESMFold")
        return _esmfold(campaign, sequence, results, target, log, warnings)

    # Co-fold the domain that is going to be docked, not the whole precursor.
    # The first co-folded run handed Boltz all 1,210 residues of EGFR when the
    # campaign's own residue_range asked for 253, and prep then threw away
    # four fifths of what had just been predicted. It cost 54 minutes against an
    # estimate of 15, and it shapes the pocket in the context of a precursor
    # nobody is docking into.
    folded_sequence, offset = _domain_of(campaign, sequence, log)
    msa = aff.ensure_msa(folded_sequence, log, warnings)
    yaml_path = out_dir / "cofold.yaml"
    yaml_path.write_text(cofold_input(folded_sequence, smiles, msa), encoding="utf-8")

    boltz_work = out_dir / "boltz"
    boltz_work.mkdir(parents=True, exist_ok=True)
    cmd = aff.boltz_command(yaml_path, boltz_work, campaign.get("fold") or {}, msa)
    # No template here, so the structure module is doing the work rather than
    # being steered, and --use_potentials is left as boltz_command sets it.
    log(f"fold: co-folding {len(sequence)} residues with the ligand, Boltz-2")
    # run_with_progress, not subprocess.run. This blocks for about a quarter of
    # an hour, and the bar only redraws when something calls update: under
    # subprocess.run it painted "0s" once and then sat there, which reads as a
    # hang and gets the run killed. On a first co-fold the silence is real work
    # -- boltz fetches its weights before it prints anything -- and now that
    # shows as its own line rather than as nothing at all.
    output, returncode = console.run_with_progress(
        cmd, log, "co-folding the protein and the ligand", estimate_s=900.0)
    try:
        (out_dir / "boltz.log").write_text(
            f"$ {' '.join(str(c) for c in cmd)}\n{output}\n", encoding="utf-8")
    except OSError:
        pass

    predicted = _newest_cif(boltz_work)
    if returncode != 0 or predicted is None:
        tail = (output or "").strip().splitlines()[-2:]
        reason = " / ".join(t.strip() for t in tail) or "no structure was written"
        warnings.append(f"Co-folding failed ({reason}); ESMFold was used instead.")
        log(f"fold: co-fold failed, falling back to ESMFold: {reason}")
        return _esmfold(campaign, sequence, results, target, log, warnings)

    aff.capture_msa(boltz_work, msa, log)
    kept, ligand_atoms = split_cofold(predicted, target, out_dir / "cofold_ligand.sdf",
                                      first_residue=offset)
    if not kept:
        warnings.append("The co-folded structure had no protein chain; ESMFold was used instead.")
        return _esmfold(campaign, sequence, results, target, log, warnings)

    # The docking box has to be moved into this structure's frame, and this is
    # the whole reason the first co-folded run put its ligand 87 A from the
    # site. The campaign's box centre was computed at Prepare time from the
    # FETCHED structure's coordinates; Boltz-2 returns its own origin, 59 A
    # away, and only 150 of 9,391 receptor atoms fell inside the box that
    # resulted. Docking then searched empty space next to the protein and
    # reported a confident pose in it.
    #
    # The centre is taken from the co-folded LIGAND, not from the pocket
    # residue list. Two reasons. The residue numbers do not transfer: a campaign
    # can carry a pocket in a crystal's numbering (1M17 counts 24 lower than
    # UniProt) while the co-folded model is numbered from the sequence, so those
    # numbers name different residues here. And the ligand needs no mapping at
    # all -- Boltz just placed this exact molecule in the pocket it built, in
    # this frame, which is the one thing that is true by construction.
    #
    # Its POSE is still discarded. Placing a box is a question about where the
    # site is, which co-folding answers well; where the ligand sits inside it is
    # the question docking is for, and the one co-folding answers at 4.5 to
    # 4.8 A.
    centre = _ligand_centre(out_dir / "cofold_ligand.sdf")
    if centre is not None:
        was = list((campaign.get("pocket") or {}).get("center") or [])
        campaign.setdefault("pocket", {})["center"] = centre
        moved = _distance(was, centre)
        log(f"fold: box centre moved into the co-folded frame, "
            f"{[round(v, 2) for v in centre]}"
            + (f" ({moved:.1f} A from the campaign's)" if moved is not None else ""))
    else:
        warnings.append(
            "The co-folded ligand could not be read, so the docking box kept the centre "
            "computed from the fetched structure. Those are different coordinate frames "
            "and the pose is very unlikely to be meaningful.")
        log("fold: WARNING could not recentre the box on the co-folded site")

    log(f"fold: co-folded, kept {kept} protein residues and set aside the "
        f"{ligand_atoms}-atom predicted pose")
    warnings.append(
        "The receptor was co-folded with the ligand and the predicted ligand pose was "
        "then discarded, so docking starts from a pocket shaped around this ligand "
        "rather than an apo one. The discarded pose is in cofold/cofold_ligand.sdf.")
    return {"folded": True, "method": "boltz2", "warnings": warnings,
            "headline": f"co-folded with the ligand, {kept} residues"}


def cofold_input(sequence: str, smiles: str, msa: dict) -> str:
    """Protein and ligand together, and deliberately no template.

    The affinity stage forces a template because it must score the exact pose MD
    produced. Here the point is the opposite: let the structure module build the
    pocket around the ligand, because that pocket is the thing being kept.
    """
    lines = ["version: 1", "sequences:", "  - protein:", "      id: A",
             f"      sequence: {sequence}"]
    if msa.get("path"):
        lines.append(f"      msa: {msa['path']}")
    lines += ["  - ligand:", "      id: L", f"      smiles: '{smiles}'"]
    return "\n".join(lines) + "\n"


def _newest_cif(work: Path) -> Optional[Path]:
    """Boltz nests its output under a results directory named for the input."""
    hits = sorted(work.rglob("*_model_0.cif")) or sorted(work.rglob("*.cif"))
    return hits[0] if hits else None


def split_cofold(cif: Path, protein_dest: Path, ligand_dest: Path,
                 first_residue: int = 1) -> tuple[int, int]:
    """Protein to model_apo.pdb, ligand to its own SDF, from one predicted complex.

    `first_residue` renumbers the output. Boltz numbers whatever sequence it was
    given from 1, so a domain folded on its own comes back as 1..253 while the
    campaign refers to it as 714..966. Renumbering here means prep's trim is a
    no-op rather than a second, conflicting opinion about which residues these
    are.
    """
    import gemmi

    structure = gemmi.read_structure(str(cif))
    structure.setup_entities()
    ligand_atoms = 0
    for model in structure:
        for chain in model:
            for residue in list(chain):
                if residue.name in ("HOH", "WAT"):
                    chain.remove_residue(residue) if hasattr(chain, "remove_residue") else None
    # Write the ligand first, while the full complex is still in hand.
    with open(ligand_dest, "w", encoding="utf-8") as fh:
        fh.write("cofold_ligand\n     GOBSMACKED\n\n")
        atoms = [(a, r) for m in structure for c in m for r in c
                 for a in r if not _is_polymer(r)]
        ligand_atoms = len(atoms)
        fh.write(f"{ligand_atoms:>3}  0  0  0  0  0  0  0  0  0999 V2000\n")
        for atom, _ in atoms:
            fh.write(f"{atom.pos.x:10.4f}{atom.pos.y:10.4f}{atom.pos.z:10.4f} "
                     f"{atom.element.name:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n")
        fh.write("M  END\n$$$$\n")

    kept = 0
    with open(protein_dest, "w", encoding="utf-8") as fh:
        serial = 0
        for model in structure:
            for chain in model:
                for residue in chain:
                    if not _is_polymer(residue):
                        continue
                    kept += 1
                    for atom in residue:
                        serial += 1
                        fh.write(
                            f"ATOM  {serial:>5d} {atom.name:<4s}{residue.name:>4s} A"
                            f"{residue.seqid.num + first_residue - 1:>4d}    "
                            f"{atom.pos.x:8.3f}{atom.pos.y:8.3f}{atom.pos.z:8.3f}"
                            f"  1.00{atom.b_iso:6.2f}          {atom.element.name:>2s}\n")
            break
        fh.write("TER\nEND\n")
    return kept, ligand_atoms


def _is_polymer(residue) -> bool:
    import gemmi
    return gemmi.find_tabulated_residue(residue.name) is not None and \
        gemmi.find_tabulated_residue(residue.name).is_amino_acid()


def _esmfold(campaign: dict, sequence: str, results: Path, target: Path,
             log, warnings: list[str]) -> dict[str, Any]:
    """ESMFold, reached directly or as the co-fold fallback."""
    log(f"fold: ESMFold on {len(sequence)} residues")
    with bar_for(log, "loading the ESMFold weights"):
        import torch
        import esm

        model = esm.pretrained.esmfold_v1()
        model = model.eval()
    device = "cuda" if torch.cuda.is_available() else ("mps" if torch.backends.mps.is_available() else "cpu")
    if device == "cpu":
        warnings.append("ESMFold ran on the CPU: this is minutes to hours rather than seconds.")
    model = model.to(device)
    chunk = chunk_size_for(len(sequence))
    if chunk:
        model.set_chunk_size(chunk)
        log(f"fold: chunk size {chunk}")

    # No progress to report from inside a single forward pass: the spinner says
    # the process is alive, which on a CPU fold is the only question being asked.
    with bar_for(log, f"folding {len(sequence)} residues on {device}"):
        with torch.no_grad():
            pdb_text = model.infer_pdb(sequence)
    target.write_text(pdb_text, encoding="utf-8")

    plddt = per_residue_plddt(pdb_text)
    (results / "plddt.json").write_text(json.dumps({
        "mean": round(sum(plddt.values()) / len(plddt), 2) if plddt else None,
        "per_residue": plddt,
    }, indent=2), encoding="utf-8")

    pocket = pocket_residue_numbers(campaign)
    low = [n for n in pocket if plddt.get(str(n), 100.0) < POCKET_PLDDT_WARNING]
    if low:
        warnings.append(
            f"{len(low)} pocket residues have pLDDT below {POCKET_PLDDT_WARNING:.0f} "
            f"({', '.join(str(n) for n in low[:8])}): the pocket geometry this run docks "
            f"into is a guess, and the scorecard will show it."
        )
    log(f"fold: done, mean pLDDT {round(sum(plddt.values()) / len(plddt), 1) if plddt else '?'}")
    mean_plddt = round(sum(plddt.values()) / len(plddt), 1) if plddt else None
    return {"folded": True, "warnings": warnings,
            "headline": f"mean pLDDT {mean_plddt}" if mean_plddt else ""}


def _ligand_centre(sdf: Path) -> Optional[list]:
    """Centroid of the co-folded ligand, in the co-folded structure's frame.

    The atom block is read by position rather than by "any line with three
    numbers on it". A molfile's counts line is `  2  0  0 ...`, which parses
    perfectly well as the coordinate (2, 0, 0) and drags the centroid toward the
    origin: with two real atoms it moved the box centre by 5 A, and the whole
    point of this function is that the box lands in the right place.
    """
    if not sdf.exists():
        return None
    lines = sdf.read_text(encoding="utf-8").splitlines()
    if len(lines) < 5:
        return None
    try:                                   # line 4 is the counts line, V2000
        count = int(lines[3][:3])
    except (ValueError, IndexError):
        return None
    points = []
    for line in lines[4:4 + count]:
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            points.append((float(parts[0]), float(parts[1]), float(parts[2])))
        except ValueError:
            continue
    if len(points) != count or not points:
        return None
    n = len(points)
    return [round(sum(p[i] for p in points) / n, 3) for i in range(3)]


def _distance(a, b) -> Optional[float]:
    if not (a and b and len(a) == 3 and len(b) == 3):
        return None
    return sum((float(x) - float(y)) ** 2 for x, y in zip(a, b)) ** 0.5


def _domain_of(campaign: dict, sequence: str, log) -> tuple[str, int]:
    """The slice of the sequence that will actually be docked, and where it starts.

    Returns the full sequence and an offset of 1 when the campaign asks for no
    range, which is the ordinary case.
    """
    span = ((campaign.get("protein") or {}).get("residue_range") or None)
    try:
        first, last = int(span[0]), int(span[1])
    except (TypeError, ValueError, IndexError):
        return sequence, 1
    if not (1 <= first <= last <= len(sequence)):
        log(f"fold: residue_range {span} does not fit a {len(sequence)}-residue "
            f"sequence, co-folding all of it")
        return sequence, 1
    domain = sequence[first - 1:last]
    log(f"fold: co-folding residues {first}-{last} ({len(domain)} of "
        f"{len(sequence)}), which is what gets docked")
    return domain, first
