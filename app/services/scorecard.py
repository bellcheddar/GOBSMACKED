"""Grades, the composite GOBSMACK score, and one sentence saying what to do.

Every metric is graded against a fixed threshold table rather than a curve, so a
score means the same thing in January and in June and across two different
targets. The composite is a weighted mean of the graded metrics; a metric that
could not be measured (no reference structure, no trajectory) is dropped and the
remaining weights are renormalised, which is stated on the card rather than
hidden in the arithmetic.

Pose validity is a gate, not a weight: a pose with a 1.8 A clash between ligand
and protein is not a B whatever else it scores, so a validity failure caps the
composite at 40.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Optional

# Grade -> points for the weighted mean. A is worth 100 so that a structure
# scored against itself comes out at exactly 100, which is the sanity check the
# test suite runs and the only self-evidently correct anchor available.
GRADE_POINTS = {"A": 100.0, "B": 85.0, "C": 70.0, "D": 55.0, "F": 25.0}
# Composite -> letter. The boundaries sit midway between adjacent grade points,
# so a run that graded B on everything reads as B.
GRADE_BOUNDARIES = (("A", 92.5), ("B", 77.5), ("C", 62.5), ("D", 47.5))
VALIDITY_FAIL_CAP = 40.0


@dataclass
class Metric:
    key: str
    label: str
    unit: str
    weight: float
    # Thresholds are read in order A, B, C, D; anything past the last is F.
    thresholds: tuple[float, float, float, float]
    lower_is_better: bool = True
    explain: Callable[[Optional[float], str], str] = lambda v, g: ""

    def grade(self, value: Optional[float]) -> Optional[str]:
        if value is None:
            return None
        for letter, limit in zip("ABCD", self.thresholds):
            if (value <= limit) if self.lower_is_better else (value >= limit):
                return letter
        return "F"


def _ligand_note(v, g):
    if v is None:
        return "No reference structure, so there is nothing to measure the pose against."
    if g == "A":
        return "The pose reproduces the crystal ligand. Nothing to fix."
    if g == "B":
        return "Close enough to be the same binding mode. Check the hinge H-bonds before trusting the details."
    if g == "C":
        return "Right pocket, wrong detail: often a flipped ring or a shifted solvent-exposed arm. Compare the interaction fingerprints."
    if g == "D":
        return "The ligand is in the pocket but the pose is not the crystal one. Re-dock with a larger box or flexible side chains."
    return "The pose is not the crystal pose. Check the box centre first: this usually means the wrong site, not a scoring failure."


def _pocket_note(v, g):
    if v is None:
        return "Not measured."
    if g in ("A", "B"):
        return "The pocket backbone matches the crystal, so the ligand had somewhere correct to bind."
    if g == "C":
        return "The pocket backbone has moved. Check whether the model was apo-like and whether MD closed the gap (see rescue)."
    return "The pocket backbone is wrong. A ligand RMSD measured against it means little: fix the model before reading anything else."


def _chi1_note(v, g):
    if v is None:
        return "Not measured."
    if g in ("A", "B"):
        return "Pocket side chains are in the crystal rotamers."
    if g == "C":
        return "Roughly a third of the pocket side chains are in a different rotamer. Flexible docking on those residues is the usual fix."
    return "Most pocket side chains are in the wrong rotamer, which is what an apo-like predicted pocket looks like. Try flex docking plus a longer equilibration."


def _jaccard_note(v, g):
    if v is None:
        return "Not measured: PLIP needs both a predicted complex and a reference holo structure."
    if g in ("A", "B"):
        return "The predicted complex makes the crystal's interactions."
    if g == "C":
        return "About half the crystal contacts are reproduced. Look at which type is missing: a lost hydrogen bond and a lost hydrophobic contact mean different things."
    return "The interaction pattern does not match the crystal even where the ligand overlaps. Check protonation and the tautomer you docked."


def _stability_note(v, g):
    if v is None:
        return "Not measured: no trajectory in the archive."
    if g in ("A", "B"):
        return "The ligand stayed where docking put it for the whole run."
    if g == "C":
        return "The ligand drifted during the run. Extend production before drawing conclusions from the final frame."
    return "The ligand left its docked pose. Either the pose was wrong or the parameters are: check the ligand force field warnings in the run log."


def _rescue_note(v, g):
    if v is None:
        return "Not measured."
    if g == "A":
        return "MD moved the pocket toward the crystal: induced fit was recovered from an apo-like start."
    if g == "B":
        return "MD improved the pocket slightly."
    if g == "C":
        return "MD left the pocket where it started."
    return "MD moved the pocket away from the crystal. Check restraint release and whether the box is large enough."


METRICS = [
    Metric("ligand_rmsd", "Ligand RMSD", "A", 30, (1.0, 2.0, 3.0, 4.0), True, _ligand_note),
    Metric("plip_jaccard", "PLIP overlap", "", 20, (0.75, 0.55, 0.40, 0.25), False, _jaccard_note),
    Metric("pocket_ca_rmsd", "Pocket Ca RMSD", "A", 15, (0.8, 1.2, 1.8, 2.5), True, _pocket_note),
    Metric("chi1_agreement", "chi1 agreement", "", 10, (0.85, 0.70, 0.55, 0.40), False, _chi1_note),
    Metric("md_drift", "Drift, last 200 ps", "A", 10, (0.5, 1.0, 1.5, 2.5), True, _stability_note),
    Metric("rescue", "MD rescue", "A", 10, (0.5, 0.2, 0.0, -0.3), False, _rescue_note),
]

VALIDITY_WEIGHT = 5.0

# PoseBusters-style checks. Each is a boolean; all must pass for the gate.
VALIDITY_CHECKS = {
    "no_clash": "No ligand-protein heavy-atom contact under 2.2 A",
    "bond_lengths": "Ligand bond lengths within 25 % of the ideal",
    "chirality": "Stereocentres unchanged from the input SMILES",
    "inside_box": "Ligand centroid inside the docking box",
}


def composite(values: dict[str, Optional[float]], validity: dict[str, bool] | None = None) -> dict[str, Any]:
    """Grade every metric, then combine what could be measured.

    `values` keys are the metric keys above. A missing or None value means the
    metric was not measurable, and it takes its weight out of the mean with it.
    """
    rows = []
    total_weight = 0.0
    total_points = 0.0
    for m in METRICS:
        value = values.get(m.key)
        grade = m.grade(value)
        rows.append({
            "key": m.key,
            "label": m.label,
            "unit": m.unit,
            "value": value,
            "grade": grade,
            "weight": m.weight,
            "note": m.explain(value, grade or ""),
            # Bar fill for the gauge: grade points, so the bar and the letter
            # never disagree.
            "fill": GRADE_POINTS.get(grade or "", 0.0),
        })
        if grade is not None:
            total_weight += m.weight
            total_points += m.weight * GRADE_POINTS[grade]

    validity = validity or {}
    checked = {k: bool(validity.get(k)) for k in VALIDITY_CHECKS if k in validity}
    validity_pass = all(checked.values()) if checked else None
    if validity_pass is not None:
        total_weight += VALIDITY_WEIGHT
        total_points += VALIDITY_WEIGHT * (100.0 if validity_pass else 0.0)

    # Without a reference structure there is nothing to verify, and a
    # composite built from the two metrics that survive (drift and validity)
    # would read as a grade for something never measured. Say unverified
    # instead of scoring 95 for a run nobody checked.
    verified = values.get("ligand_rmsd") is not None
    score = round(total_points / total_weight, 1) if (total_weight and verified) else None
    capped = False
    if score is not None and validity_pass is False and score > VALIDITY_FAIL_CAP:
        score = VALIDITY_FAIL_CAP
        capped = True

    return {
        "metrics": rows,
        "score": score,
        "grade": score_to_grade(score),
        "validity": {
            "checks": [{"key": k, "label": VALIDITY_CHECKS[k], "pass": v} for k, v in checked.items()],
            "pass": validity_pass,
            "capped": capped,
        },
        "verified": verified,
        "measured": int(total_weight),
        "unmeasured": [r["label"] for r in rows if r["grade"] is None],
        "label": _score_label(score, validity_pass, capped, verified),
    }


def score_to_grade(score: Optional[float]) -> Optional[str]:
    if score is None:
        return None
    for letter, floor in GRADE_BOUNDARIES:
        if score >= floor:
            return letter
    return "F"


def _score_label(score, validity_pass, capped, verified=True) -> str:
    if not verified:
        return "Unverified: no reference structure, so there is no ligand RMSD to grade"
    if score is None:
        return "Nothing measurable in this archive"
    if capped:
        return "GOBSMACK score, capped by a validity failure"
    if validity_pass is None:
        return "GOBSMACK score, validity not checked"
    return "GOBSMACK score, validity gate passed" if validity_pass else "GOBSMACK score, validity gate failed"


# ---------------------------------------------------------------------------
# Pose validity
# ---------------------------------------------------------------------------

def check_validity(complex_path, ligand_smiles: str = "", box_centre=None, box=None,
                   ligand_ccd: str | None = None) -> dict[str, bool]:
    """PoseBusters-style checks on a docked complex, using RDKit and geometry.

    This is not PoseBusters (which is Apache-2.0 but pulls a large dependency
    tree); it is the four of its checks that catch the failures docking actually
    produces here, computed from the same files everything else reads.
    """
    import gemmi
    import numpy as np
    from rdkit import Chem

    from .superpose import ligand_atoms

    out: dict[str, bool] = {}
    lig = ligand_atoms(complex_path, ligand_ccd)
    if not lig:
        return out
    lig_xyz = np.array([xyz for _, xyz in lig])

    st = gemmi.read_structure(str(complex_path))
    st.setup_entities()
    st.remove_waters()
    from .superpose import is_ligand_residue

    prot = []
    for ch in st[0]:
        for res in ch:
            if is_ligand_residue(res) or res.is_water():
                continue
            for a in res:
                if a.element != gemmi.Element("H"):
                    prot.append([a.pos.x, a.pos.y, a.pos.z])
    if prot:
        prot_xyz = np.array(prot)
        d = np.linalg.norm(lig_xyz[:, None, :] - prot_xyz[None, :, :], axis=2)
        out["no_clash"] = bool(d.min() >= 2.2)

    # Bond lengths: any heavy-atom bond outside 0.9 to 1.9 A is broken geometry
    # whatever the element pair.
    dists = np.linalg.norm(lig_xyz[:, None, :] - lig_xyz[None, :, :], axis=2)
    np.fill_diagonal(dists, np.inf)
    nearest = dists.min(axis=1)
    out["bond_lengths"] = bool((nearest > 0.9).all() and (nearest < 1.9).all())

    if box_centre is not None and box is not None:
        centre = np.array(box_centre, dtype=float)
        half = np.array(box, dtype=float) / 2.0
        out["inside_box"] = bool((np.abs(lig_xyz.mean(axis=0) - centre) <= half).all())

    if ligand_smiles:
        template = Chem.MolFromSmiles(ligand_smiles)
        if template is not None:
            wanted = Chem.FindMolChiralCenters(template, useLegacyImplementation=False)
            # No stereocentres in the input means nothing to preserve, which
            # passes rather than being unmeasurable.
            out["chirality"] = True if not wanted else _chirality_preserved(
                complex_path, ligand_smiles, ligand_ccd, wanted)
    return out


def _chirality_preserved(complex_path, smiles: str, ccd, wanted) -> bool:
    from rdkit import Chem

    from .superpose import _ligand_mol

    mol = _ligand_mol(complex_path, ccd, smiles)
    if mol is None:
        return False
    try:
        Chem.AssignStereochemistryFrom3D(mol)
        got = dict(Chem.FindMolChiralCenters(mol, useLegacyImplementation=False))
    except Exception:
        return False
    return all(got.get(idx) == tag for idx, tag in wanted)
