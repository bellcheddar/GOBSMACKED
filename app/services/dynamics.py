"""Trajectory panels and the two numbers the scorecard takes from them.

Most of this is reading `traj/summary.json`, which the bundle wrote on the
machine that had the trajectory. Two things are computed here instead:

* **drift**, the ligand's mean RMSD over the last 200 ps minus its mean over the
  first 200 ps: whether the pose held.
* **ligand RMSD to the crystal ligand, per frame**, which needs the reference
  structure and therefore could not be computed in the bundle. This is the
  induced-fit trace: a prediction that starts 3 A out and relaxes to 1.5 A is a
  different result from one that starts at 1.5 A and stays there.

The second one loads the trajectory on the droplet, so it is bounded: an atom x
frame budget, and a straight refusal past it rather than an OOM kill.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import numpy as np

# 3 million atom-frames is 300 frames of 10,000 atoms: more than a 1 ns run at
# the default 10 ps frame interval ever packs, and small enough that MDTraj's
# copies fit alongside a gunicorn worker on a 3.8 GB box.
BUDGET_ATOM_FRAMES = 3_000_000
DRIFT_WINDOW_PS = 200.0


def summarise(summary: dict) -> dict[str, Any]:
    """The dynamics panels, with the derived numbers the scorecard reads."""
    times = summary.get("times_ps") or []
    lig = summary.get("ligand_rmsd_pose1") or []
    out: dict[str, Any] = {
        "frames": len(times),
        "times_ps": times,
        "ligand_rmsd_pose1": lig,
        "protein_ca_rmsd": summary.get("protein_ca_rmsd") or [],
        "pocket_ca_rmsd": summary.get("pocket_ca_rmsd") or [],
        "pocket_volume": summary.get("pocket_volume") or [],
        "rmsf": summary.get("rmsf") or {},
        "contacts": summary.get("contacts") or {},
        "pocket_residues": summary.get("pocket_residues") or [],
    }
    out["drift"] = drift(times, lig)
    out["contact_persistence"] = contact_persistence(out["contacts"])
    return out


def drift(times_ps: list[float], values: list[float],
          window: float = DRIFT_WINDOW_PS) -> Optional[float]:
    """Mean over the last `window` ps minus the mean over the first `window` ps.

    Positive means the ligand moved away from where docking put it. Returns None
    when the run is too short for two non-overlapping windows: a number computed
    from windows that share frames would flatter every short run.
    """
    if not times_ps or len(times_ps) != len(values) or len(times_ps) < 4:
        return None
    t = np.asarray(times_ps, dtype=float)
    v = np.asarray(values, dtype=float)
    if (t[-1] - t[0]) < 2 * window:
        return None
    first = v[t <= t[0] + window]
    last = v[t >= t[-1] - window]
    if len(first) == 0 or len(last) == 0:
        return None
    return round(float(last.mean() - first.mean()), 3)


def contact_persistence(contacts: dict) -> list[dict]:
    """Fraction of frames each pocket residue is in contact with the ligand."""
    residues = contacts.get("residues") or []
    matrix = contacts.get("matrix") or []
    if not residues or not matrix:
        return []
    m = np.asarray(matrix, dtype=float)
    if m.ndim != 2 or m.shape[0] != len(residues):
        return []
    fractions = m.mean(axis=1)
    return [{"residue": int(r), "persistence": round(float(f), 3)}
            for r, f in zip(residues, fractions)]


def rescue(pocket_rmsd_model: Optional[float],
           pocket_rmsd_md_final: Optional[float]) -> Optional[float]:
    """Pocket Ca RMSD-to-reference before MD minus after MD.

    Positive means MD moved the pocket toward the crystal, which is the whole
    induced-fit claim: an ESMFold pocket is apo-like, and the question is
    whether docking plus relaxation recovers the holo shape.
    """
    if pocket_rmsd_model is None or pocket_rmsd_md_final is None:
        return None
    return round(pocket_rmsd_model - pocket_rmsd_md_final, 3)


def ligand_rmsd_to_reference(topology: str | Path, trajectory: str | Path,
                             reference_path: str | Path, rot, trans,
                             reference_ccd: Optional[str] = None,
                             ligand_resname: str = "LIG",
                             smiles: str = "") -> dict[str, Any]:
    """Per-frame heavy-atom RMSD of the docked ligand to the crystal ligand.

    Each frame is moved by the pocket superposition computed once on the final
    frame, then compared to the reference ligand. Symmetry is handled the same
    way as the single-structure measurement, by matching the two graphs once and
    reusing the atom correspondence for every frame: re-running the graph match
    per frame would cost more than the trajectory read and give the same answer.
    """
    import mdtraj as md

    from .superpose import apply_transform, _ligand_mol

    topology, trajectory = Path(topology), Path(trajectory)
    if not trajectory.exists() or not topology.exists():
        return {"error": "No trajectory in the archive, so there is no per-frame trace."}

    top = md.load_topology(str(topology))
    n_atoms = top.n_atoms
    try:
        n_frames = int(md.load(str(trajectory), top=top).n_frames)
    except Exception as exc:
        return {"error": f"Could not read the trajectory: {exc}"}
    if n_atoms * n_frames > BUDGET_ATOM_FRAMES:
        return {"error": f"That trajectory is {n_atoms * n_frames:,} atom-frames, past this "
                         f"server's {BUDGET_ATOM_FRAMES:,} budget. Re-pack it with a longer "
                         f"frame interval."}

    selection = top.select(f"resname {ligand_resname} and not element H")
    if len(selection) == 0:
        return {"error": f"No residue named {ligand_resname} in the trajectory topology."}

    traj = md.load(str(trajectory), top=top, atom_indices=selection)
    coords = traj.xyz * 10.0                      # MDTraj works in nm, this app in A

    ref_mol = _ligand_mol(reference_path, reference_ccd, smiles)
    if ref_mol is None:
        return {"error": "No reference ligand to measure against."}
    ref_xyz = np.array([list(ref_mol.GetConformer().GetAtomPosition(i))
                        for i in range(ref_mol.GetNumAtoms())])
    if ref_xyz.shape[0] != coords.shape[1]:
        # Different heavy-atom counts: compare centroids, and say so rather than
        # silently reporting a number computed on a different definition.
        series = [round(float(np.linalg.norm(
            apply_transform(frame, rot, trans).mean(axis=0) - ref_xyz.mean(axis=0))), 3)
            for frame in coords]
        return {"series": series, "approximate": True,
                "note": "Atom counts differ between the docked and crystal ligands, so this "
                        "trace is the distance between their centroids."}

    series = []
    for frame in coords:
        moved = apply_transform(frame, rot, trans)
        series.append(round(float(np.sqrt(((moved - ref_xyz) ** 2).sum(axis=1).mean())), 3))
    return {"series": series, "approximate": False}
