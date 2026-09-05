"""Stage 1: fold, or skip.

Skipped whenever `model_apo.pdb` is already in the bundle, which is the usual
case: the server fetches an AlphaFold DB, ESM Atlas, PDB or user-supplied
structure and ships it. Folding here is the fallback for a sequence nothing has
a model for, and it is the reason this bundle wants a GPU at all.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional

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
