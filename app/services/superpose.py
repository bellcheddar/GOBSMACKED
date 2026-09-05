"""Superposition and the geometric half of the scorecard.

Everything here is measured after superposing on **pocket Ca atoms**, not on the
whole chain. A model can be excellent at the binding site and 6 A out at a
disordered terminus; whole-chain superposition spreads that error into the
pocket and inflates the ligand RMSD, which is the one number this app exists to
report honestly. Whole-chain TM-score is reported alongside, as context.

Residue numbering is never assumed to match. The reference chain and the model
are aligned by sequence first, and every measurement walks that mapping, so an
AFDB model numbered by UniProt and a crystal numbered by its own construct
compare correctly.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

import gemmi
import numpy as np

from .. import config

# Chi1 is defined by these four atoms; the third depends on the residue type.
CHI1_THIRD = {
    "ARG": "CG", "ASN": "CG", "ASP": "CG", "CYS": "SG", "GLN": "CG", "GLU": "CG",
    "HIS": "CG", "ILE": "CG1", "LEU": "CG", "LYS": "CG", "MET": "CG", "PHE": "CG",
    "PRO": "CG", "SER": "OG", "THR": "OG1", "TRP": "CG", "TYR": "CG", "VAL": "CG1",
}
CHI1_TOLERANCE_DEG = 40.0


def is_amino_acid(res: gemmi.Residue) -> bool:
    """Is this residue a standard amino acid?

    gemmi puts this on the tabulated component rather than on the residue
    itself, so it goes through the chemical component dictionary lookup.
    """
    return gemmi.find_tabulated_residue(res.name) is not None and \
        gemmi.find_tabulated_residue(res.name).is_amino_acid()


@dataclass
class Chain:
    """One polymer chain, indexed the way the measurements need it."""
    name: str
    sequence: str
    numbers: list[int]                          # residue number per sequence position
    residues: dict[int, gemmi.Residue] = field(default_factory=dict)

    def ca(self, number: int) -> Optional[np.ndarray]:
        res = self.residues.get(number)
        if res is None:
            return None
        atom = res.find_atom("CA", "*")
        return np.array([atom.pos.x, atom.pos.y, atom.pos.z]) if atom else None

    def heavy_sidechain(self, number: int) -> list[tuple[str, np.ndarray]]:
        res = self.residues.get(number)
        if res is None:
            return []
        out = []
        for atom in res:
            if atom.element == gemmi.Element("H"):
                continue
            if atom.name in ("N", "CA", "C", "O", "OXT"):
                continue
            out.append((atom.name, np.array([atom.pos.x, atom.pos.y, atom.pos.z])))
        return out

    def chi1(self, number: int) -> Optional[float]:
        res = self.residues.get(number)
        if res is None:
            return None
        third = CHI1_THIRD.get(res.name.upper())
        if third is None:                        # GLY and ALA have no chi1
            return None
        pts = []
        for name in ("N", "CA", "CB", third):
            atom = res.find_atom(name, "*")
            if atom is None:
                return None
            pts.append(np.array([atom.pos.x, atom.pos.y, atom.pos.z]))
        return dihedral(*pts)


def load_chain(path: str | Path, chain_name: Optional[str] = None) -> Optional[Chain]:
    """The named polymer chain (or the longest one) with waters and ligands removed."""
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_waters()
    best: Optional[Chain] = None
    for ch in st[0]:
        residues = {}
        seq = []
        numbers = []
        for res in ch:
            if not is_amino_acid(res):
                continue
            residues[res.seqid.num] = res
            seq.append(gemmi.find_tabulated_residue(res.name).one_letter_code.upper())
            numbers.append(res.seqid.num)
        if not residues:
            continue
        cand = Chain(name=ch.name, sequence="".join(seq), numbers=numbers, residues=residues)
        if chain_name and ch.name == chain_name:
            return cand
        if best is None or len(cand.residues) > len(best.residues):
            best = cand
    return best


def ligand_atoms(path: str | Path, ccd: Optional[str] = None,
                 chain_name: Optional[str] = None) -> list[tuple[str, np.ndarray]]:
    """Heavy atoms of the named ligand, or of the largest non-water heteroatom group."""
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_waters()
    best: list[tuple[str, np.ndarray]] = []
    for ch in st[0]:
        if chain_name and ch.name != chain_name:
            continue
        for res in ch:
            if is_amino_acid(res) or res.is_water():
                continue
            if ccd and res.name.upper() != ccd.upper():
                continue
            atoms = [(a.name, np.array([a.pos.x, a.pos.y, a.pos.z]))
                     for a in res if a.element != gemmi.Element("H")]
            if len(atoms) > len(best):
                best = atoms
    return best


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def dihedral(p0, p1, p2, p3) -> float:
    """Signed dihedral in degrees.

    b0 is p0 - p1, not p1 - p0. The other sign convention gives an answer
    exactly 180 degrees out, which looks entirely plausible and quietly ruins
    every rotamer comparison downstream.
    """
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation and translation taking `mobile` onto `target` (both N x 3)."""
    mc = mobile.mean(axis=0)
    tc = target.mean(axis=0)
    h = (mobile - mc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, tc - rot @ mc


def apply_transform(coords: np.ndarray, rot: np.ndarray, trans: np.ndarray) -> np.ndarray:
    return (rot @ coords.T).T + trans


def rmsd(a: np.ndarray, b: np.ndarray) -> Optional[float]:
    if len(a) == 0 or len(a) != len(b):
        return None
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


# ---------------------------------------------------------------------------
# Sequence mapping
# ---------------------------------------------------------------------------

def align_numbering(model: Chain, reference: Chain) -> dict[int, int]:
    """model residue number -> reference residue number, by sequence alignment.

    Global alignment with BLOSUM62 through biotite. Only aligned pairs of
    identical or substituted residues are kept, so an insertion in either
    construct simply has no partner rather than shifting everything after it.
    """
    from biotite.sequence import ProteinSequence
    from biotite.sequence.align import SubstitutionMatrix, align_optimal

    # Unknown residues are coerced to alanine rather than dropped: removing a
    # character would shift every index after it, and the trace indexes back
    # into `numbers` to recover residue numbers.
    try:
        s1 = ProteinSequence(model.sequence.replace("X", "A"))
        s2 = ProteinSequence(reference.sequence.replace("X", "A"))
    except Exception:
        return {}
    matrix = SubstitutionMatrix.std_protein_matrix()
    alignment = align_optimal(s1, s2, matrix, gap_penalty=(-10, -1), terminal_penalty=False)[0]
    mapping: dict[int, int] = {}
    for i, j in alignment.trace:
        if i < 0 or j < 0:
            continue
        mapping[model.numbers[i]] = reference.numbers[j]
    return mapping


def pocket_residues(reference_path: str | Path, reference: Chain, ccd: Optional[str],
                    radius: float = config.POCKET_RADIUS) -> list[int]:
    """Reference residue numbers with a heavy atom within `radius` of its ligand."""
    lig = ligand_atoms(reference_path, ccd)
    if not lig:
        return []
    lig_xyz = np.array([xyz for _, xyz in lig])
    out = []
    for num, res in reference.residues.items():
        coords = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res
                           if a.element != gemmi.Element("H")])
        if len(coords) == 0:
            continue
        d = np.linalg.norm(coords[:, None, :] - lig_xyz[None, :, :], axis=2)
        if d.min() <= radius:
            out.append(num)
    return sorted(out)


# ---------------------------------------------------------------------------
# The comparison
# ---------------------------------------------------------------------------

def compare(model_path: str | Path, reference_path: str | Path,
            reference_ccd: Optional[str] = None, reference_chain: Optional[str] = None,
            model_chain: Optional[str] = None,
            model_ligand_ccd: Optional[str] = None,
            ligand_smiles: str = "") -> dict[str, Any]:
    """One predicted structure against the crystal, superposed on the pocket."""
    model = load_chain(model_path, model_chain)
    ref = load_chain(reference_path, reference_chain)
    if model is None or ref is None:
        return {"error": "Could not read a protein chain from one of the structures."}

    mapping = align_numbering(model, ref)
    if not mapping:
        return {"error": "The model and the reference share no alignable sequence."}
    inverse = {v: k for k, v in mapping.items()}

    pocket_ref = pocket_residues(reference_path, ref, reference_ccd)
    if not pocket_ref:
        return {"error": "No ligand found in the reference, so there is no pocket to superpose on."}

    pairs = [(inverse[r], r) for r in pocket_ref if r in inverse]
    mob, tgt = [], []
    for m_num, r_num in pairs:
        a, b = model.ca(m_num), ref.ca(r_num)
        if a is not None and b is not None:
            mob.append(a)
            tgt.append(b)
    if len(mob) < 4:
        return {"error": f"Only {len(mob)} pocket Ca atoms could be paired: too few to superpose."}

    rot, trans = kabsch(np.array(mob), np.array(tgt))
    pocket_ca_rmsd = rmsd(apply_transform(np.array(mob), rot, trans), np.array(tgt))

    # Side chains, walking the same pairs and matching atoms by name.
    sc_model, sc_ref = [], []
    chi_hits, chi_total = 0, 0
    displacements: list[dict] = []
    for m_num, r_num in pairs:
        m_atoms = dict(model.heavy_sidechain(m_num))
        r_atoms = dict(ref.heavy_sidechain(r_num))
        shared = [n for n in m_atoms if n in r_atoms]
        if shared:
            m_xyz = apply_transform(np.array([m_atoms[n] for n in shared]), rot, trans)
            r_xyz = np.array([r_atoms[n] for n in shared])
            sc_model.extend(m_xyz)
            sc_ref.extend(r_xyz)
            per_residue = float(np.sqrt(((m_xyz - r_xyz) ** 2).sum(axis=1).mean()))
        else:
            per_residue = None
        chi_m, chi_r = model.chi1(m_num), ref.chi1(r_num)
        chi_delta = None
        if chi_m is not None and chi_r is not None:
            chi_total += 1
            chi_delta = abs((chi_m - chi_r + 180) % 360 - 180)
            if chi_delta <= CHI1_TOLERANCE_DEG:
                chi_hits += 1
        displacements.append({
            "reference_residue": r_num,
            "model_residue": m_num,
            "name": ref.residues[r_num].name,
            "sidechain_rmsd": round(per_residue, 2) if per_residue is not None else None,
            "chi1_model": round(chi_m, 1) if chi_m is not None else None,
            "chi1_reference": round(chi_r, 1) if chi_r is not None else None,
            "chi1_delta": round(chi_delta, 1) if chi_delta is not None else None,
        })

    result: dict[str, Any] = {
        "pocket_residues_reference": pocket_ref,
        "pocket_pairs": len(pairs),
        "pocket_ca_atoms": len(mob),
        "pocket_ca_rmsd": round(pocket_ca_rmsd, 3) if pocket_ca_rmsd is not None else None,
        "pocket_sc_rmsd": round(rmsd(np.array(sc_model), np.array(sc_ref)), 3) if sc_model else None,
        "chi1_agreement": round(chi_hits / chi_total, 3) if chi_total else None,
        "chi1_residues": chi_total,
        "displacements": sorted(
            [d for d in displacements if d["sidechain_rmsd"] is not None],
            key=lambda d: -d["sidechain_rmsd"],
        )[:10],
        "transform": {"rotation": rot.tolist(), "translation": trans.tolist()},
        # Model residue number -> reference residue number, as string keys so
        # the whole result survives a round trip through JSON. The interaction
        # fingerprints need this to compare contacts across two numberings.
        "numbering_map": {str(k): v for k, v in mapping.items()},
        "model_chain": model.name,
        "reference_chain": ref.name,
    }

    # Ligand RMSD, symmetry corrected, in the superposed frame.
    lig_model = ligand_atoms(model_path, model_ligand_ccd)
    lig_ref = ligand_atoms(reference_path, reference_ccd)
    if lig_model and lig_ref:
        result["ligand_rmsd"] = ligand_rmsd(
            model_path, reference_path, rot, trans,
            model_ccd=model_ligand_ccd, reference_ccd=reference_ccd, smiles=ligand_smiles,
        )
        result["ligand_atoms_model"] = len(lig_model)
        result["ligand_atoms_reference"] = len(lig_ref)
    else:
        result["ligand_rmsd"] = None

    result["tm_score"] = tm_score(model, ref, mapping)
    return result


def ligand_rmsd(model_path: str | Path, reference_path: str | Path,
                rot: np.ndarray, trans: np.ndarray, model_ccd: Optional[str] = None,
                reference_ccd: Optional[str] = None, smiles: str = "") -> Optional[float]:
    """Heavy-atom RMSD of the docked ligand to the crystal ligand.

    Symmetry matters: a phenyl ring flipped 180 degrees is the same pose, and
    naive atom-order RMSD calls it a 2 A error. RDKit's CalcRMS enumerates the
    graph automorphisms and takes the best, and unlike GetBestRMS it does not
    re-align the molecules, which would throw away the protein superposition
    this whole measurement rests on.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign

    probe = _ligand_mol(model_path, model_ccd, smiles)
    ref = _ligand_mol(reference_path, reference_ccd, smiles)
    if probe is None or ref is None:
        return None

    conf = probe.GetConformer()
    for i in range(probe.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        x, y, z = apply_transform(np.array([[p.x, p.y, p.z]]), rot, trans)[0]
        conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(float(x), float(y), float(z)))

    try:
        return round(float(rdMolAlign.CalcRMS(probe, ref)), 3)
    except (RuntimeError, ValueError):
        # Different heavy-atom graphs (a different protonation state, a
        # truncated ligand in the crystal): fall back to a centroid distance,
        # clearly labelled by the caller as approximate.
        p = np.array([list(probe.GetConformer().GetAtomPosition(i))
                      for i in range(probe.GetNumAtoms())])
        r = np.array([list(ref.GetConformer().GetAtomPosition(i))
                      for i in range(ref.GetNumAtoms())])
        return round(float(np.linalg.norm(p.mean(axis=0) - r.mean(axis=0))), 3)


def _ligand_mol(path: str | Path, ccd: Optional[str], smiles: str):
    """Extract one ligand as an RDKit molecule with bond orders from the SMILES."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_waters()
    st.remove_hydrogens()
    keep = None
    for ch in st[0]:
        for res in ch:
            if is_amino_acid(res) or res.is_water():
                continue
            if ccd and res.name.upper() != ccd.upper():
                continue
            if keep is None or len(res) > len(keep):
                keep = res
    if keep is None:
        return None

    # Written by gemmi rather than by hand: a PDB line is fixed-width, and an
    # atom name one column out of place parses to a molecule with zero atoms
    # that raises nothing at all until something asks it for a conformer.
    single = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    chain.add_residue(keep.clone())
    model.add_chain(chain)
    single.add_model(model)
    single.setup_entities()
    block = single.make_pdb_string()

    mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=True)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    if smiles:
        template = Chem.MolFromSmiles(smiles)
        if template is not None:
            try:
                mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
            except (ValueError, RuntimeError):
                pass                              # keep the raw connectivity
    return mol


def tm_score(model: Chain, reference: Chain, mapping: dict[int, int]) -> Optional[float]:
    """Whole-chain TM-score, for context rather than for the grade."""
    try:
        from tmtools import tm_align
    except ImportError:
        return None
    m_xyz, m_seq, r_xyz, r_seq = [], [], [], []
    for m_num, r_num in mapping.items():
        a, b = model.ca(m_num), reference.ca(r_num)
        if a is None or b is None:
            continue
        m_xyz.append(a)
        r_xyz.append(b)
        m_seq.append(gemmi.find_tabulated_residue(model.residues[m_num].name).one_letter_code.upper())
        r_seq.append(gemmi.find_tabulated_residue(reference.residues[r_num].name).one_letter_code.upper())
    if len(m_xyz) < 10:
        return None
    try:
        res = tm_align(np.array(m_xyz), np.array(r_xyz), "".join(m_seq), "".join(r_seq))
        return round(float(max(res.tm_norm_chain1, res.tm_norm_chain2)), 3)
    except Exception:
        return None
