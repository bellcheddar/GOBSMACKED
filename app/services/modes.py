"""Binding-mode classification: GATEKEEPER.

Two families get a real answer and everything else gets an honest one.

**Kinases** are labelled from the KLIFS 85-residue pocket: DFG-in / out / inter
by the Modi and Dunbrack distance criteria, alphaC-in / out from the beta3-Lys to
alphaC-Glu salt bridge, and a Type I / I-and-a-half / II / allosteric label from
which subpockets the ligand occupies.

**GPCRs** are labelled from GPCRdb generic numbering: which site the ligand sits
in, and whether the receptor's own microswitches read active-like or
inactive-like.

Both classifiers run on the predicted complex and on the reference crystal with
the same code, so the comparison is like for like: a difference in the label is
a difference in the structure, not a difference in method.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import gemmi
import numpy as np

from .superpose import Chain, is_amino_acid, load_chain, ligand_atoms

# --- Kinase thresholds -----------------------------------------------------
# Modi and Dunbrack (2019, PNAS, doi:10.1073/pnas.1814279116) classify the DFG
# state from two distances to the DFG-Phe ring centre atom (CZ):
#   D1: alphaC-Glu(+4) CA to DFG-Phe CZ
#   D2: beta3-Lys CA to DFG-Phe CZ
# In DFG-in the Phe sits in the back pocket, far from the beta3 lysine and close
# to the alphaC; in DFG-out it has swapped into the ATP site.
DFG_D1_CUT = 11.0
DFG_D2_CUT = 11.0
# The beta3-Lys NZ to alphaC-Glu carboxylate salt bridge. 4.0 A is a formed
# ionic bond; beyond it the helix has swung out.
SALT_BRIDGE_CUT = 4.0
CONTACT_CUT = 4.0
# A subpocket counts as occupied only when the ligand touches at least two of
# its residues. One residue at 4 A is a graze: erlotinib in 1M17 reaches exactly
# one alphaC residue (Glu at KLIFS 24) and a one-contact rule labels a textbook
# type I inhibitor as type I-and-a-half. Single-residue subpockets (the
# gatekeeper) necessarily use a threshold of one.
MIN_SUBPOCKET_CONTACTS = 2

# --- GPCR thresholds -------------------------------------------------------
# TM3-TM6, measured CA(3.50) to CA(6.30): the ionic lock opening as TM6 swings
# out is the clearest single number for class A activation. The cut was placed
# by measuring six structures with this code rather than taken from a paper
# using a different atom pair:
#
#   inactive   rhodopsin 1GZM 8.7   A2A 3EML 9.7   b2AR 2RH1 11.2
#   active     metarhodopsin 3PQR 14.7   A2A 5G53 18.5   b2AR 3SN6 19.0
#
# 13 A sits in the empty band between the two groups. Note that quoted "ionic
# lock" distances of 3 to 4 A are guanidinium-to-carboxylate, not CA-to-CA, and
# are not comparable to these.
TM3_TM6_CUT = 13.0

ORTHOSTERIC_GENERIC = {
    "3.32", "3.33", "3.36", "3.37", "5.42", "5.43", "5.46", "5.47",
    "6.48", "6.51", "6.52", "6.55", "7.38", "7.39", "7.42", "7.43",
}
INTRACELLULAR_GENERIC = {
    "3.50", "3.53", "3.54", "6.29", "6.30", "6.33", "6.36", "6.37",
    "7.53", "7.55", "7.56", "8.47", "8.49",
}
VESTIBULE_GENERIC = {
    "2.60", "2.63", "3.28", "3.29", "5.39", "6.58", "6.59", "7.32", "7.35",
}

# Subpockets, as sets of KLIFS pocket positions. Regions alone are too coarse:
# every ATP-site ligand touches the xDFG region and the beta3 strand, so a
# subpocket defined as a whole region calls erlotinib a type I-and-a-half. What
# separates the types is whether the ligand passes the gatekeeper (back pocket
# I, lined by the alphaC positions) and whether it reaches the deep allosteric
# pocket that only opens when DFG swings out (back pocket II).
SUBPOCKETS = {
    "adenine": [15, 17, 45, 46, 47, 48, 51, 75, 77, 81],
    "hinge": [46, 47, 48],
    "gatekeeper": [45],
    "front pocket": [6, 7, 8, 9, 49, 50, 51, 52],
    "back pocket I": [24, 25, 26, 27, 28, 29, 30, 31],
    "back pocket II": [32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
    "alphaC": list(range(20, 31)),
    "DFG": [81, 82, 83],
    "catalytic loop": list(range(68, 76)),
}


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def map_sequence_to_structure(chain: Chain, sequence: str) -> dict[int, int]:
    """Position in `sequence` (1-based) -> residue number in `chain`.

    KLIFS and GPCRdb both number against the canonical UniProt sequence, while a
    crystal is numbered however its depositor chose. Aligning once here means
    every position downstream is looked up rather than assumed.
    """
    from biotite.sequence import ProteinSequence
    from biotite.sequence.align import SubstitutionMatrix, align_optimal

    if not sequence or not chain.sequence:
        return {}
    try:
        s1 = ProteinSequence(sequence.replace("X", "A"))
        s2 = ProteinSequence(chain.sequence.replace("X", "A"))
    except Exception:
        return {}
    matrix = SubstitutionMatrix.std_protein_matrix()
    alignment = align_optimal(s1, s2, matrix, gap_penalty=(-10, -1), terminal_penalty=False)[0]
    out: dict[int, int] = {}
    for i, j in alignment.trace:
        if i < 0 or j < 0:
            continue
        out[i + 1] = chain.numbers[j]
    return out


def _atom(chain: Chain, number: Optional[int], name: str) -> Optional[np.ndarray]:
    if number is None:
        return None
    res = chain.residues.get(number)
    if res is None:
        return None
    atom = res.find_atom(name, "*")
    return np.array([atom.pos.x, atom.pos.y, atom.pos.z]) if atom else None


def _residue_heavy(chain: Chain, number: Optional[int]) -> np.ndarray:
    res = chain.residues.get(number) if number else None
    if res is None:
        return np.zeros((0, 3))
    return np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res
                     if a.element != gemmi.Element("H")])


def _contacts(chain: Chain, numbers: list[int], lig_xyz: np.ndarray,
              cutoff: float = CONTACT_CUT) -> list[int]:
    """Which of `numbers` have a heavy atom within `cutoff` of the ligand."""
    if len(lig_xyz) == 0:
        return []
    hits = []
    for n in numbers:
        coords = _residue_heavy(chain, n)
        if len(coords) == 0:
            continue
        d = np.linalg.norm(coords[:, None, :] - lig_xyz[None, :, :], axis=2)
        if d.min() <= cutoff:
            hits.append(n)
    return hits


def _dist(a: Optional[np.ndarray], b: Optional[np.ndarray]) -> Optional[float]:
    if a is None or b is None:
        return None
    return round(float(np.linalg.norm(a - b)), 2)


def solvent_exposed_fraction(chain: Chain, lig_xyz: np.ndarray, cutoff: float = 4.5,
                             min_neighbours: int = 2) -> Optional[float]:
    """Fraction of ligand heavy atoms with almost no protein around them.

    A proxy for the solvent-exposed part of the ligand, computed from the same
    coordinates as everything else rather than from an SASA library: an atom
    with fewer than two protein heavy atoms within 4.5 A is sticking out.
    """
    if len(lig_xyz) == 0:
        return None
    prot = np.array([[a.pos.x, a.pos.y, a.pos.z]
                     for res in chain.residues.values() for a in res
                     if a.element != gemmi.Element("H")])
    if len(prot) == 0:
        return None
    d = np.linalg.norm(lig_xyz[:, None, :] - prot[None, :, :], axis=2)
    neighbours = (d <= cutoff).sum(axis=1)
    return round(float((neighbours < min_neighbours).mean()), 3)


# ---------------------------------------------------------------------------
# Kinase
# ---------------------------------------------------------------------------

def classify_kinase(structure_path: str | Path, sequence: str, pocket_map: dict[str, int],
                    chain_name: Optional[str] = None, ligand_ccd: Optional[str] = None,
                    hinge_hbonds: Optional[int] = None) -> dict[str, Any]:
    """DFG state, alphaC state, subpocket occupancy and a type label."""
    chain = load_chain(structure_path, chain_name)
    if chain is None:
        return {"error": "No protein chain in that structure."}
    seq_to_struct = map_sequence_to_structure(chain, sequence)

    def at(klifs_pos: int) -> Optional[int]:
        seq_num = pocket_map.get(str(klifs_pos))
        return seq_to_struct.get(seq_num) if seq_num else None

    lys, glu = at(17), at(24)
    gatekeeper, dfg_asp, dfg_phe = at(45), at(81), at(82)
    lig = ligand_atoms(structure_path, ligand_ccd, chain_name)
    lig_xyz = np.array([xyz for _, xyz in lig]) if lig else np.zeros((0, 3))

    # alphaC-Glu(+4) is the residue four along the helix, which is where the
    # Modi and Dunbrack D1 distance is measured from.
    glu_plus4 = glu + 4 if glu is not None else None

    d1 = _dist(_atom(chain, glu_plus4, "CA"), _atom(chain, dfg_phe, "CZ"))
    d2 = _dist(_atom(chain, lys, "CA"), _atom(chain, dfg_phe, "CZ"))
    dfg = _dfg_label(d1, d2)

    salt = None
    lys_nz = _atom(chain, lys, "NZ")
    for name in ("OE1", "OE2", "CD"):
        candidate = _dist(lys_nz, _atom(chain, glu, name))
        if candidate is not None and (salt is None or candidate < salt):
            salt = candidate
    alphac = None if salt is None else ("in" if salt <= SALT_BRIDGE_CUT else "out")

    occupancy: dict[str, bool] = {}
    contacts_by_region: dict[str, list[int]] = {}
    for subpocket, positions in SUBPOCKETS.items():
        numbers = [n for n in (at(p) for p in positions) if n is not None]
        hits = _contacts(chain, sorted(set(numbers)), lig_xyz)
        threshold = 1 if len(positions) == 1 else MIN_SUBPOCKET_CONTACTS
        occupancy[subpocket] = len(hits) >= threshold
        contacts_by_region[subpocket] = hits

    if hinge_hbonds is None:
        hinge_hbonds = _hinge_hbonds(chain, [at(p) for p in range(46, 49)], lig)

    exposed = solvent_exposed_fraction(chain, lig_xyz)
    label, reason = _kinase_type(occupancy, dfg, bool(lig))

    return {
        "family": "kinase",
        "label": label,
        "reason": reason,
        "dfg": dfg,
        "dfg_d1": d1,
        "dfg_d2": d2,
        "alphac": alphac,
        "salt_bridge": salt,
        "occupancy": occupancy,
        "contacts": contacts_by_region,
        "hinge_hbonds": hinge_hbonds,
        "gatekeeper_residue": _residue_label(chain, gatekeeper),
        "gatekeeper_contact": bool(contacts_by_region.get("gatekeeper")),
        "solvent_exposed_fraction": exposed,
        "positions": {"beta3_lys": lys, "alphaC_glu": glu, "gatekeeper": gatekeeper,
                      "dfg_asp": dfg_asp, "dfg_phe": dfg_phe},
    }


def _dfg_label(d1: Optional[float], d2: Optional[float]) -> Optional[str]:
    if d1 is None or d2 is None:
        return None
    if d1 <= DFG_D1_CUT and d2 >= DFG_D2_CUT:
        return "in"
    if d1 > DFG_D1_CUT and d2 < DFG_D2_CUT:
        return "out"
    return "inter"


def _kinase_type(occupancy: dict[str, bool], dfg: Optional[str], has_ligand: bool) -> tuple[str, str]:
    if not has_ligand:
        return "apo", "No ligand in this structure, so there is no binding mode to label."
    if not occupancy.get("hinge"):
        return "allosteric", "The ligand makes no hinge contact, so it is not competing with ATP."
    back = occupancy.get("back pocket I") or occupancy.get("back pocket II")
    if dfg == "out":
        where = "into the allosteric pocket DFG-out opens" if back else "while holding the hinge"
        return "II", f"DFG is out and the ligand reaches {where}: type II."
    if back and dfg == "in":
        return "I 1/2", "Hinge plus the back pocket with DFG in: type I and a half."
    return "I", "Hinge contact, DFG in, and the ligand stays in front of the gatekeeper: type I."


def _hinge_hbonds(chain: Chain, hinge_numbers: list[Optional[int]],
                  ligand: list[tuple[str, np.ndarray]], cutoff: float = 3.5) -> int:
    """Hydrogen bonds between the ligand and the hinge, by heavy-atom geometry.

    Donor and acceptor are not distinguished: a nitrogen or oxygen pair inside
    3.5 A across the hinge is the interaction every kinase paper counts, and
    PLIP's angle-aware count replaces this one when it is available.
    """
    polar_lig = [xyz for name, xyz in ligand if name[:1] in ("N", "O")]
    if not polar_lig:
        return 0
    count = 0
    for num in hinge_numbers:
        res = chain.residues.get(num) if num else None
        if res is None:
            continue
        for atom in res:
            if atom.element.name not in ("N", "O"):
                continue
            pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
            if any(np.linalg.norm(pos - p) <= cutoff for p in polar_lig):
                count += 1
    return count


def _residue_label(chain: Chain, number: Optional[int]) -> Optional[str]:
    res = chain.residues.get(number) if number else None
    if res is None:
        return None
    return f"{res.name.title()}{number}"


# ---------------------------------------------------------------------------
# GPCR
# ---------------------------------------------------------------------------

def classify_gpcr(structure_path: str | Path, sequence: str, generic: dict[str, int],
                  chain_name: Optional[str] = None,
                  ligand_ccd: Optional[str] = None) -> dict[str, Any]:
    """Which site the ligand occupies, and where the microswitches sit."""
    chain = load_chain(structure_path, chain_name)
    if chain is None:
        return {"error": "No protein chain in that structure."}
    seq_to_struct = map_sequence_to_structure(chain, sequence)

    struct_generic: dict[str, int] = {}
    for gn, seq_num in generic.items():
        struct_num = seq_to_struct.get(seq_num)
        if struct_num is not None:
            struct_generic[gn] = struct_num
    inverse = {v: k for k, v in struct_generic.items()}

    lig = ligand_atoms(structure_path, ligand_ccd, chain_name)
    lig_xyz = np.array([xyz for _, xyz in lig]) if lig else np.zeros((0, 3))
    contact_numbers = _contacts(chain, sorted(chain.residues), lig_xyz)
    contact_generic = sorted({inverse[n] for n in contact_numbers if n in inverse})

    orthosteric = len(set(contact_generic) & ORTHOSTERIC_GENERIC)
    vestibule = len(set(contact_generic) & VESTIBULE_GENERIC)
    intracellular = len(set(contact_generic) & INTRACELLULAR_GENERIC)

    if not lig:
        site, reason = "apo", "No ligand in this structure."
    elif orthosteric >= 3:
        site, reason = "orthosteric", f"{orthosteric} contacts to consensus orthosteric positions."
    elif intracellular >= 2:
        site, reason = "intracellular", f"{intracellular} contacts at the cytoplasmic ends of TM6 and TM7."
    elif vestibule >= 2:
        site, reason = "vestibule", f"{vestibule} contacts in the extracellular vestibule."
    elif contact_generic:
        site, reason = "lipid-facing", "Contacts are outside every named site: the ligand is on the helix bundle's outer face."
    else:
        site, reason = "unassigned", "No contacts to any generically numbered residue."

    tm3_tm6 = _dist(_atom(chain, struct_generic.get("3.50"), "CA"),
                    _atom(chain, struct_generic.get("6.30"), "CA"))
    state = None if tm3_tm6 is None else ("active-like" if tm3_tm6 >= TM3_TM6_CUT else "inactive-like")

    toggle = struct_generic.get("6.48")
    toggle_chi1 = chain.chi1(toggle) if toggle else None
    pif = _dist(_atom(chain, struct_generic.get("3.40"), "CA"),
                _atom(chain, struct_generic.get("6.44"), "CA"))
    npxxy = [struct_generic.get(gn) for gn in ("7.49", "7.50", "7.51", "7.52", "7.53")]
    sodium = struct_generic.get("2.50")
    sodium_contact = bool(_contacts(chain, [sodium], lig_xyz)) if sodium else False

    return {
        "family": "gpcr",
        "label": site,
        "reason": reason,
        "state": state,
        "tm3_tm6": tm3_tm6,
        "toggle_chi1": round(toggle_chi1, 1) if toggle_chi1 is not None else None,
        "toggle_residue": _residue_label(chain, toggle),
        "pif_distance": pif,
        "sodium_contact": sodium_contact,
        "sodium_residue": _residue_label(chain, sodium),
        "contact_generic": contact_generic,
        "npxxy_residues": [n for n in npxxy if n],
        "positions": struct_generic,
    }


def npxxy_rmsd(model_path: str | Path, reference_path: str | Path,
               model_generic: dict[str, int], reference_generic: dict[str, int],
               model_chain: Optional[str] = None,
               reference_chain: Optional[str] = None) -> Optional[float]:
    """Backbone RMSD over NPxxY (7.49 to 7.53) after superposing on TM3 and TM6.

    Measured between the two structures rather than against a stored template:
    a template would need a class and a numbering scheme of its own, and the
    reference crystal is already the thing the model is being judged against.
    """
    from .superpose import kabsch, apply_transform, rmsd

    m = load_chain(model_path, model_chain)
    r = load_chain(reference_path, reference_chain)
    if m is None or r is None:
        return None
    anchor = [gn for gn in ("3.46", "3.50", "3.54", "6.30", "6.34", "6.40", "6.44", "6.48")
              if gn in model_generic and gn in reference_generic]
    mob = [_atom(m, model_generic[gn], "CA") for gn in anchor]
    tgt = [_atom(r, reference_generic[gn], "CA") for gn in anchor]
    pairs = [(a, b) for a, b in zip(mob, tgt) if a is not None and b is not None]
    if len(pairs) < 4:
        return None
    rot, trans = kabsch(np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs]))

    loop = [gn for gn in ("7.49", "7.50", "7.51", "7.52", "7.53")
            if gn in model_generic and gn in reference_generic]
    m_xyz, r_xyz = [], []
    for gn in loop:
        for name in ("N", "CA", "C", "O"):
            a = _atom(m, model_generic[gn], name)
            b = _atom(r, reference_generic[gn], name)
            if a is not None and b is not None:
                m_xyz.append(a)
                r_xyz.append(b)
    if len(m_xyz) < 8:
        return None
    moved = apply_transform(np.array(m_xyz), rot, trans)
    value = rmsd(moved, np.array(r_xyz))
    return round(value, 2) if value is not None else None


def compare_modes(predicted: dict, reference: Optional[dict]) -> dict[str, Any]:
    """The verdict card: does the predicted binding mode match the crystal's?"""
    if not reference or reference.get("error"):
        return {"match": None, "verdict": "Binding mode unverified",
                "detail": "No reference structure, so there is nothing to compare the label to."}
    if predicted.get("error"):
        return {"match": None, "verdict": "Binding mode unverified",
                "detail": predicted["error"]}

    p_label = predicted.get("label")
    r_label = reference.get("label")
    if p_label == r_label:
        detail = f"Both the prediction and the crystal are {p_label}."
        extras = _mode_differences(predicted, reference)
        if extras:
            detail += " " + extras
        return {"match": True, "verdict": "Binding mode matches", "detail": detail}

    detail = f"The prediction is {p_label}; the crystal is {r_label}."
    extras = _mode_differences(predicted, reference)
    if extras:
        detail += " " + extras
    return {"match": False, "verdict": "Binding mode differs", "detail": detail}


def _mode_differences(predicted: dict, reference: dict) -> str:
    bits = []
    if predicted.get("family") == "kinase":
        if predicted.get("dfg") != reference.get("dfg"):
            bits.append(f"DFG is {predicted.get('dfg')} in the prediction and {reference.get('dfg')} in the crystal.")
        if predicted.get("alphac") != reference.get("alphac"):
            bits.append(f"alphaC is {predicted.get('alphac')} against {reference.get('alphac')}.")
        lost = [k for k, v in (reference.get("occupancy") or {}).items()
                if v and not (predicted.get("occupancy") or {}).get(k)]
        if lost:
            bits.append("The prediction does not reach " + ", ".join(lost) + ".")
    else:
        if predicted.get("state") != reference.get("state"):
            bits.append(f"The receptor reads {predicted.get('state')} against {reference.get('state')}.")
    return " ".join(bits)
