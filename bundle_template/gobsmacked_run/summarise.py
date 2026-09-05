"""Stage 5: turn the trajectory into numbers, and pack the archive.

Everything expensive happens here rather than on the server: the server reads
`summary.json` and, at most, the small solute-only DCD. Per frame this writes

  * ligand heavy-atom RMSD to the docked pose
  * protein Ca RMSD and pocket Ca RMSD
  * pocket volume, by voxel counting
  * a residue-by-frame contact matrix at 4 A

and per residue an RMSF. Frames are aligned on the pocket Ca atoms first, so a
ligand that stays put in a tumbling protein reads as staying put.
"""

from __future__ import annotations

import json
import tarfile
from pathlib import Path
from typing import Any, Optional

import numpy as np

from .console import bar_for
from .noise import hush_c_stdout

CONTACT_CUTOFF_NM = 0.4          # MDTraj works in nanometres throughout
VOXEL_SPACING_A = 1.0
PROBE_RADIUS_A = 1.4


def run(campaign: dict, work: Path, results: Path, log) -> dict[str, Any]:
    import mdtraj as md

    warnings: list[str] = []
    traj_dir = results / "traj"
    topology_path = traj_dir / "topology.pdb"
    dcd_path = traj_dir / "traj.dcd"

    # MDTraj's DCD plugin announces itself from C on every open, two lines that
    # say nothing and land in the middle of the bar. See noise.hush_c_stdout.
    with hush_c_stdout() as said:
        top = md.load_topology(str(topology_path))
        if dcd_path.exists():
            traj = md.load(str(dcd_path), top=top)
    for line in said:
        log(f"summarise: {line}")
    if not dcd_path.exists():
        traj = md.load(str(topology_path))
        warnings.append("No production trajectory: the summary describes the final frame alone.")
    log(f"summarise: {traj.n_frames} frames, {traj.n_atoms} atoms")

    pocket_numbers = pocket_residue_numbers(campaign)
    ligand_atoms = select_ligand(top)
    protein_ca = top.select("protein and name CA")
    pocket_ca = pocket_ca_indices(top, pocket_numbers)
    if len(pocket_ca) < 4:
        pocket_ca = protein_ca
        warnings.append("Fewer than four pocket Ca atoms matched the campaign residue list, so "
                        "the trajectory was aligned on the whole chain instead.")

    reimaged = reimage(traj)

    # Two alignments, deliberately, because they answer different questions.
    #
    # The ligand and pocket metrics are measured in the POCKET frame: a hinge
    # that swings 4 A at the far end of the protein would otherwise read as the
    # ligand moving.
    #
    # RMSF is measured in the WHOLE-CHAIN frame, because md.rmsf does not align
    # anything: it reports fluctuation about the mean of whatever frames it is
    # given, so a trajectory aligned on 51 pocket atoms turns a whole-body
    # rotation about the pocket into a lever arm at the termini. On this run
    # that put the free N-terminus at 72 A of "fluctuation" inside an 84 A box.
    # Aligned on every Ca instead, the same trajectory peaks at 3.1 A.
    rmsf_traj = traj[:]
    if len(protein_ca):
        rmsf_traj.superpose(rmsf_traj, frame=0, atom_indices=protein_ca)
    traj.superpose(traj, frame=0, atom_indices=pocket_ca)

    interval_ps = float((campaign.get("md") or {}).get("frame_interval_ps", 10))
    times = [round(i * interval_ps, 2) for i in range(traj.n_frames)]

    # Named so the bar can say which measurement is running. Pocket volume is
    # the slow one by a wide margin, voxel counting every frame, and being able
    # to see that it is the one still going is the whole point of naming them.
    measurements = [
        ("ligand RMSD", lambda: rmsd_series(traj, ligand_atoms)),
        ("protein Ca RMSD", lambda: rmsd_series(traj, protein_ca)),
        ("pocket Ca RMSD", lambda: rmsd_series(traj, pocket_ca)),
        ("per-residue fluctuation", lambda: rmsf(rmsf_traj, protein_ca, top)),
        ("pocket volume", lambda: pocket_volume_series(traj, top, ligand_atoms, pocket_ca)),
        ("contact matrix", lambda: contact_matrix(traj, top, ligand_atoms, pocket_numbers)),
    ]
    keys = ["ligand_rmsd_pose1", "protein_ca_rmsd", "pocket_ca_rmsd", "rmsf",
            "pocket_volume", "contacts"]
    summary: dict[str, Any] = {
        "frames": traj.n_frames,
        "reimaged": reimaged,
        "times_ps": times,
        "pocket_residues": pocket_numbers,
    }
    with bar_for(log, f"measuring {traj.n_frames} frames",
                 total=len(measurements)) as bar:
        for index, (label, measure) in enumerate(measurements):
            bar.update(index, note=label)
            summary[keys[index]] = measure()
        bar.update(len(measurements), note="done")

    (traj_dir / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    last = summary["ligand_rmsd_pose1"][-1] if summary["ligand_rmsd_pose1"] else None
    log(f"summarise: ligand RMSD {last if last is not None else '?'} A at the last frame")
    return {"warnings": warnings, "frames": traj.n_frames,
            "headline": (f"{traj.n_frames} frames, ligand {last} A from the docked pose"
                         if last is not None else f"{traj.n_frames} frames")}


# Residues MDTraj will happily call "protein" that are not the ligand, plus the
# solvent and ion names a solute-only trajectory can still contain.
NOT_LIGAND = set(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL "
    "HID HIE HIP CYX ASH GLH LYN MSE HOH WAT SOL NA CL K MG ZN CA SOD CLA POT".split()
)


def select_ligand(top) -> "np.ndarray":
    """Heavy atoms of the docked ligand.

    Not `not protein`: MDTraj classifies UNK as a protein residue, and UNK is
    exactly what OpenMM names the ligand when it is merged in from an OpenFF
    topology. That selection returned zero atoms on the first real run and the
    ligand RMSD series came back empty, with nothing in the log but a question
    mark. The largest residue that is not a standard amino acid, water or ion is
    the ligand.
    """
    best, best_size = None, 0
    for residue in top.residues:
        if residue.name.upper() in NOT_LIGAND:
            continue
        heavy = [a.index for a in residue.atoms if a.element.symbol != "H"]
        if len(heavy) > best_size:
            best, best_size = heavy, len(heavy)
    return np.array(best if best else [], dtype=int)


def reimage(traj):
    """Undo periodic wrapping, anchoring on the largest molecule.

    MDTraj's own anchor heuristic looks for a molecule with MORE atoms than the
    largest one, which is unsatisfiable for a protein plus one ligand and raises
    rather than returning. Naming the anchor sidesteps that. Without this a
    ligand that crosses a boundary sits a box length from its protein, and the
    contact map still looks correct because compute_contacts applies the minimum
    image convention internally while the distances do not.
    """
    if traj.unitcell_vectors is None:
        return False
    try:
        molecules = traj.topology.find_molecules()
        if not molecules:
            return False
        traj.image_molecules(inplace=True, anchor_molecules=[max(molecules, key=len)])
        return True
    except (ValueError, RuntimeError, IndexError):
        return False


def rmsd_series(traj, indices) -> list[float]:
    import mdtraj as md

    if len(indices) == 0:
        return []
    values = md.rmsd(traj, traj, frame=0, atom_indices=indices) * 10.0
    return [round(float(v), 3) for v in values]


def rmsf(traj, ca_indices, top) -> dict[str, list]:
    import mdtraj as md

    if len(ca_indices) == 0 or traj.n_frames < 2:
        return {"residues": [], "values": []}
    values = md.rmsf(traj, traj, 0, atom_indices=ca_indices) * 10.0
    residues = [int(top.atom(i).residue.resSeq) for i in ca_indices]
    return {"residues": residues, "values": [round(float(v), 3) for v in values],
            "aligned_on": "all protein Ca"}


def pocket_ca_indices(top, numbers: list[int]) -> np.ndarray:
    wanted = set(numbers)
    return np.array([atom.index for atom in top.atoms
                     if atom.name == "CA" and atom.residue.resSeq in wanted], dtype=int)


def pocket_residue_numbers(campaign: dict) -> list[int]:
    out = []
    for item in ((campaign.get("pocket") or {}).get("residues") or []):
        tail = str(item).split(":")[-1]
        try:
            out.append(int(tail))
        except ValueError:
            continue
    return sorted(set(out))


def contact_matrix(traj, top, ligand_atoms, pocket_numbers: list[int]) -> dict[str, Any]:
    """Residue x frame, 1 where any heavy atom is within 4 A of the ligand."""
    if len(ligand_atoms) == 0 or not pocket_numbers:
        return {"residues": [], "matrix": []}
    wanted = set(pocket_numbers)
    residues = sorted({atom.residue.resSeq for atom in top.atoms
                       if atom.residue.resSeq in wanted and atom.residue.is_protein})
    ligand_xyz = traj.xyz[:, ligand_atoms, :]
    matrix = []
    for number in residues:
        indices = np.array([atom.index for atom in top.atoms
                            if atom.residue.resSeq == number and atom.element.symbol != "H"],
                           dtype=int)
        if len(indices) == 0:
            matrix.append([0] * traj.n_frames)
            continue
        residue_xyz = traj.xyz[:, indices, :]
        row = []
        for frame in range(traj.n_frames):
            d = np.linalg.norm(residue_xyz[frame][:, None, :] - ligand_xyz[frame][None, :, :],
                               axis=2)
            row.append(int(d.min() <= CONTACT_CUTOFF_NM))
        matrix.append(row)
    return {"residues": residues, "matrix": matrix}


def pocket_volume_series(traj, top, ligand_atoms, pocket_ca) -> list[float]:
    """Pocket volume per frame, by counting empty voxels.

    A 1 A grid is laid over the pocket's bounding box and a voxel counts when it
    is further than a probe radius from every protein heavy atom and within 8 A
    of a pocket Ca atom. It is a relative measure: what matters is whether the
    pocket opens or closes during the run, not the absolute number, and it needs
    no fpocket dependency in the bundle.
    """
    if len(pocket_ca) == 0:
        return []
    heavy = np.array([atom.index for atom in top.atoms
                      if atom.element.symbol != "H" and atom.residue.is_protein], dtype=int)
    if len(heavy) == 0:
        return []
    volumes = []
    for frame in range(traj.n_frames):
        centres = traj.xyz[frame, pocket_ca, :] * 10.0
        lo = centres.min(axis=0) - 6.0
        hi = centres.max(axis=0) + 6.0
        axes = [np.arange(lo[i], hi[i], VOXEL_SPACING_A) for i in range(3)]
        grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
        if len(grid) > 400_000:
            # A pocket this size is a cleft, and the count would cost more than
            # the whole rest of the summary.
            grid = grid[:: max(1, len(grid) // 400_000)]
        protein_xyz = traj.xyz[frame, heavy, :] * 10.0
        near_pocket = (np.linalg.norm(grid[:, None, :] - centres[None, :, :], axis=2).min(axis=1)
                       <= 8.0)
        candidates = grid[near_pocket]
        if len(candidates) == 0:
            volumes.append(0.0)
            continue
        empty = (np.linalg.norm(candidates[:, None, :] - protein_xyz[None, :, :], axis=2).min(axis=1)
                 > PROBE_RADIUS_A + 1.7)
        volumes.append(round(float(empty.sum() * VOXEL_SPACING_A ** 3), 1))
    return volumes


def pack(results: Path, dest: Path, log) -> Path:
    """results/ into results.tar.gz, with results/ as the single top directory.

    `dest` belongs beside run.py, at the top of the unpacked bundle, and not
    inside `results/` where it used to be written. Two reasons, and the second
    one is a bug rather than a preference:

    * it is the one file the run produces that a person has to find, and asking
      them to open a directory called `results` to find it among twenty PDBs is
      one step of digging too many
    * packing it into the directory being packed means a second `--only
      summarise` puts the previous archive inside the new one, which is
      invisible except as an archive that doubles in size every time
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(dest, "w:gz") as tar:
        for path in sorted(results.rglob("*")):
            if path.is_dir() or "__pycache__" in path.parts:
                continue
            tar.add(path, arcname=f"results/{path.relative_to(results)}")
    log(f"summarise: wrote {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
    return dest
