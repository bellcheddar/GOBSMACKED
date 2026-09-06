"""Every pose docking produced, not just the one that went forward.

The archive has always carried all of them, in `poses/poses.sdf`, and the page
has always shown one: pose 1, merged into the complex that goes on to MD. The
other nine were scored, ranked and then never looked at again.

That is a waste of the most interesting thing docking produces. The question a
verification tool should be able to answer is not only "is the top pose right"
but "was the right pose ever found and then mis-ranked", and those are different
failures with different fixes. A scoring function that ranks the crystal pose
seventh is a scoring problem; a search that never visited it at all is a
sampling problem, and no amount of rescoring will help.

So every pose is measured against the crystal ligand, in the same superposed
frame as everything else on the page, and the best-by-RMSD is named alongside
the best-by-score.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import numpy as np


def overlay(poses_sdf: Path, scores: list[dict], reference_path: Optional[Path],
            rot, trans, dest_dir: Path, reference_ccd: Optional[str] = None,
            smiles: str = "", receptor_path: Optional[Path] = None) -> dict[str, Any]:
    """Superpose every pose, measure each against the crystal, write the files.

    Two files rather than ten: the top pose alone, and the rest together. The
    distinction the eye needs is "the one that went forward" against "the ones
    that did not", and ten colours in one panel is a rainbow nobody can read.
    """
    from rdkit import Chem

    if not poses_sdf.exists():
        return {"error": "This archive holds no poses.sdf, so there is nothing to overlay."}

    supplier = Chem.SDMolSupplier(str(poses_sdf), sanitize=False, removeHs=False)
    molecules = [mol for mol in supplier if mol is not None]
    if not molecules:
        return {"error": "poses.sdf held no readable molecules."}

    reference_mol = None
    if reference_path is not None:
        from .superpose import _ligand_mol

        reference_mol = _ligand_mol(reference_path, reference_ccd, smiles)

    # The receptor, once, for the closest-contact check. Read from the
    # superposed copy so it is in the same frame the poses are moved into.
    protein = _protein_atoms(receptor_path)

    dest_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    moved: list = []
    for index, mol in enumerate(molecules):
        placed = _transformed(mol, rot, trans)
        moved.append(placed)
        # PandaDock's own numbers, read from the SD tags rather than from
        # scores.csv: the tags carry more than the CSV does, they are in every
        # archive already written, and they need no rerun to recover.
        tags = _tags(mol)
        fallback = scores[index] if index < len(scores) else {}
        rows.append({
            "rank": int(tags.get("rank") or index + 1),
            "pose_id": fallback.get("pose_id") or f"pose{index + 1}",
            "score": _number(tags.get("score_kcal_per_mol"), fallback.get("score")),
            "gnn_energy": _number(tags.get("energy_gnn_energy")),
            "vina_energy": _number(tags.get("energy_vina_energy")),
            "gnn_pec50": _number(tags.get("energy_gnn_pec50"), fallback.get("gnn_affinity")),
            "confidence": _number(tags.get("confidence")),
            "rmsd": _rmsd_to(placed, reference_mol),
            # Two more, because RMSD alone cannot say WHICH failure this is.
            # Centroid distance separates "wrong pocket" from "right pocket";
            # best-fit RMSD, which superposes the two ligands freely, separates
            # "wrong orientation" from "wrong conformation". A pose can sit on
            # top of the crystal ligand's centre, fold the same way, and still
            # score 8 A because it is turned around.
            "centroid_distance": _centroid_distance(placed, reference_mol),
            "shape_rmsd": _rmsd_to(placed, reference_mol, best_fit=True),
            "closest_contact": _closest_contact(placed, protein),
        })

    # The crystal ligand alone. Loading reference.pdb here draws the whole
    # 1M17 chain in amber, and the panel becomes two proteins with the poses
    # lost between them: the comparison this panel exists to make is
    # ligand against ligand.
    reference_file = None
    if reference_mol is not None:
        reference_file = dest_dir / "reference_ligand.pdb"
        _write([reference_mol], reference_file)
        reference_file = reference_file.name

    top = dest_dir / "poses_top.pdb"
    rest = dest_dir / "poses_rest.pdb"
    _write(moved[:1], top)
    if len(moved) > 1:
        _write(moved[1:], rest, first_rank=2)

    measured = [r for r in rows if r["rmsd"] is not None]
    best_by_rmsd = min(measured, key=lambda r: r["rmsd"]) if measured else None
    return {
        "rows": rows,
        "top_file": top.name,
        "reference_file": reference_file,
        "rest_file": rest.name if len(moved) > 1 else None,
        "count": len(rows),
        "best_by_rmsd": best_by_rmsd,
        "verdict": _verdict(rows, best_by_rmsd),
    }


def _tags(mol) -> dict:
    try:
        return mol.GetPropsAsDict()
    except Exception:                              # noqa: BLE001 - a pose with no tags
        return {}


def _number(*values) -> Optional[float]:
    """The first of these that is a number, rounded for display."""
    for value in values:
        if value is None or value == "":
            continue
        try:
            return round(float(value), 3)
        except (TypeError, ValueError):
            continue
    return None


def _protein_atoms(receptor_path: Optional[Path]):
    """Heavy atoms of the receptor, for the closest-contact check."""
    if receptor_path is None or not Path(receptor_path).exists():
        return None
    import gemmi

    from .superpose import is_ligand_residue

    structure = gemmi.read_structure(str(receptor_path))
    coords = [[atom.pos.x, atom.pos.y, atom.pos.z]
              for model in structure for chain in model for residue in chain
              if not is_ligand_residue(residue) and residue.name != "HOH"
              for atom in residue if atom.element.name != "H"]
    return np.array(coords) if coords else None


def _closest_contact(mol, protein) -> Optional[float]:
    """Nearest heavy-atom approach to the receptor.

    The cheap half of a validity check, and the half that catches the failure
    that matters: anything under about 2.2 A is a clash, and a pose that clashes
    is not a pose however well it scored. Computed per pose because the
    scorecard's own check runs on one complex and these are ten loose ligands.
    """
    if protein is None:
        return None
    conformer = mol.GetConformer()
    coords = np.array([list(conformer.GetAtomPosition(i))
                       for i in range(mol.GetNumAtoms())
                       if mol.GetAtomWithIdx(i).GetAtomicNum() > 1])
    if not len(coords):
        return None
    distances = np.linalg.norm(coords[:, None, :] - protein[None, :, :], axis=-1)
    return round(float(distances.min()), 2)


def _transformed(mol, rot, trans):
    """The same rotation and translation the rest of the page was measured in.

    A pose drawn in the model's own frame would sit metres from the crystal
    ligand it is supposed to be compared with, which reads as a catastrophic
    prediction rather than as an unapplied transform.
    """
    from rdkit.Chem import Mol

    placed = Mol(mol)
    conformer = placed.GetConformer()
    coords = np.array([list(conformer.GetAtomPosition(i))
                       for i in range(placed.GetNumAtoms())])
    if rot is not None and trans is not None:
        coords = coords @ np.asarray(rot).T + np.asarray(trans)
    for i, xyz in enumerate(coords):
        conformer.SetAtomPosition(i, [float(v) for v in xyz])
    return placed


def _centroid_distance(mol, reference_mol) -> Optional[float]:
    """How far the pose's centre is from the crystal ligand's.

    The question this answers is "is it even in the right pocket", which RMSD
    conflates with "is it the right way round".
    """
    if reference_mol is None:
        return None
    try:
        a = _heavy_coords(mol)
        b = _heavy_coords(reference_mol)
        return round(float(np.linalg.norm(a.mean(axis=0) - b.mean(axis=0))), 2)
    except Exception:                              # noqa: BLE001
        return None


def _heavy_coords(mol) -> np.ndarray:
    conformer = mol.GetConformer()
    return np.array([list(conformer.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())
                     if mol.GetAtomWithIdx(i).GetAtomicNum() > 1])


def _rmsd_to(mol, reference_mol, best_fit: bool = False) -> Optional[float]:
    """Symmetry-corrected, and never re-aligned.

    CalcRMS rather than GetBestRMS for the reason it is used everywhere else
    here: GetBestRMS superposes the two molecules first, which would throw away
    the protein superposition the whole measurement rests on and report how
    similar the two shapes are rather than whether the pose is in the right
    place.
    """
    if reference_mol is None:
        return None
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    from rdkit.Chem import AllChem

    try:
        probe = Chem.Mol(mol)
        # Bond orders first. A pose read from SDF with sanitize=False carries no
        # aromaticity or bond orders, so CalcRMS finds no substructure match
        # between it and the crystal ligand and every pose comes back
        # unmeasured: "No sub-structure match found between the reference and
        # probe mol". Templating from the reference gives the two graphs the
        # same chemistry, which is what the automorphism search needs.
        probe = Chem.RemoveHs(probe, sanitize=False)
        ref = Chem.RemoveHs(Chem.Mol(reference_mol), sanitize=False)
        try:
            probe = AllChem.AssignBondOrdersFromTemplate(ref, probe)
        except Exception:                          # noqa: BLE001 - fall back to as-read
            pass
        if best_fit:
            # GetBestRMS superposes the two molecules before measuring, which
            # throws away the placement entirely and reports how alike the two
            # SHAPES are. That is exactly the question here, and exactly why it
            # must never be used for the main RMSD.
            return round(float(rdMolAlign.GetBestRMS(Chem.Mol(probe), Chem.Mol(ref))), 3)
        return round(float(rdMolAlign.CalcRMS(probe, ref)), 3)
    except Exception:                              # noqa: BLE001 - a pose that will not match
        return None


def _write(molecules: list, dest: Path, first_rank: int = 1) -> None:
    """Poses into one PDB, each as its own chain, so a single load shows them all.

    Chains rather than models: Mol* shows model 1 of a multi-model file and
    hides the rest, which is exactly the opposite of what an overlay is for.
    """
    from rdkit import Chem

    blocks = []
    for offset, mol in enumerate(molecules):
        rank = first_rank + offset
        chain = chr(ord("A") + (rank - 1) % 26)
        block = Chem.MolToPDBBlock(mol, flavor=2)
        for line in block.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                # Column 22 is the chain, and 18-20 the residue name: one
                # chain per pose, all named LIG so the viewer treats them
                # alike.
                line = line[:17] + "LIG" + line[20:21] + chain + line[22:]
                blocks.append(line)
        blocks.append("TER")
    dest.write_text("\n".join(blocks) + "\nEND\n", encoding="utf-8")


def _verdict(rows: list[dict], best: Optional[dict]) -> str:
    """Which failure this is, which is not the same question as how big it is.

    An earlier version read RMSD alone and told a run whose poses sat 2.4 A from
    the crystal ligand's centre, fully inside the docking box, that "the search
    never visited the right place, check the box centre". The box was correct
    and the advice was wrong. Three numbers are needed to tell these apart:

    * centroid distance -- is it in the right pocket at all
    * in-place RMSD     -- is it the right pose
    * best-fit RMSD     -- is the molecule folded the same way, ignoring where
                           it sits

    Right pocket with a high in-place RMSD and a low best-fit one is an
    orientation failure, and that is a sampling and scoring problem inside the
    site rather than a box problem outside it.
    """
    measured = [r for r in rows if r["rmsd"] is not None]
    if not measured:
        return "No crystal ligand to measure these against."
    top = rows[0]
    if top["rmsd"] is not None and top["rmsd"] <= 2.0:
        return (f"The top pose is {top['rmsd']:.2f} A from the crystal ligand: the search "
                f"found it and the scoring function agreed.")

    if best and best["rmsd"] <= 2.0:
        return (f"The search found the crystal pose and ranked it {best['rank']} of "
                f"{len(rows)}: pose {best['rank']} is {best['rmsd']:.2f} A out while the "
                f"top-ranked one is {top['rmsd']:.2f} A. That is a scoring failure rather "
                f"than a sampling one, so rescoring is the thing to try.")

    centres = [r["centroid_distance"] for r in rows if r["centroid_distance"] is not None]
    shapes = [r["shape_rmsd"] for r in rows if r["shape_rmsd"] is not None]
    nearest_centre = min(centres) if centres else None
    best_shape = min(shapes) if shapes else None

    if nearest_centre is not None and nearest_centre <= 4.0:
        text = (f"Every pose is in the right pocket and none is the right pose: the closest "
                f"centre is {nearest_centre:.1f} A from the crystal ligand's, while the best "
                f"in-place RMSD is {min(r['rmsd'] for r in measured):.1f} A.")
        if best_shape is not None and best_shape <= 3.0:
            text += (f" Superposed freely, the ligand matches the crystal to "
                     f"{best_shape:.1f} A, so it is folded about right and turned the wrong "
                     f"way round. Raise the exhaustiveness and tighten the box before "
                     f"suspecting the scoring function; a flexible-side-chain run is the "
                     f"next thing after that.")
        return text

    return (f"No pose came within 2 A of the crystal ligand, and the closest centre is "
            f"{nearest_centre:.1f} A away, so these are not in the crystal's pocket at all. "
            f"Check the box centre before anything else." if nearest_centre is not None else
            f"No pose came within 2 A of the crystal ligand: the closest of {len(rows)} is "
            f"pose {best['rank']} at {best['rmsd']:.2f} A.")
