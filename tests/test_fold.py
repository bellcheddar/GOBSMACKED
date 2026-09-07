"""The fold stage's choice of method, and what co-folding keeps.

Co-folding is the measured recommendation from the five-structure experiment: it
was the only starting structure whose docked pose came back at rank 1. What
makes it work is asymmetric, and the tests below pin that asymmetry down --
Boltz-2's protein is kept and its ligand is discarded, because on the same
target its pocket was the closest of any starting structure while its own ligand
pose was 4.5 to 4.8 A out.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bundle_template"))

from gobsmacked_run import fold  # noqa: E402


def test_the_default_is_still_esmfold():
    """A campaign written before co-folding existed must not change behaviour."""
    assert fold.method_of({}) == "esmfold"
    assert fold.method_of({"fold": {}}) == "esmfold"
    assert fold.method_of({"fold": {"method": "esmfold"}}) == "esmfold"


def test_co_folding_is_opt_in_and_case_insensitive():
    assert fold.method_of({"fold": {"method": "boltz2"}}) == "boltz2"
    assert fold.method_of({"fold": {"method": "Boltz2"}}) == "boltz2"


def test_an_unknown_method_falls_back_rather_than_raising():
    """A bundle can outlive the runner that reads it. Refusing to fold because
    the campaign named a method this copy has never heard of would strand a run
    that ESMFold could have completed."""
    assert fold.method_of({"fold": {"method": "alphafold9"}}) == "esmfold"
    assert fold.method_of({"fold": {"method": None}}) == "esmfold"


def test_the_cofold_yaml_has_no_template():
    """The affinity stage forces a template because it must score one exact
    pose. Co-folding wants the opposite: the structure module has to build the
    pocket around the ligand, because that pocket is the thing being kept."""
    text = fold.cofold_input("ACDEFGH", "CCO", {"path": None})
    assert "templates:" not in text and "force:" not in text
    assert "sequence: ACDEFGH" in text
    assert "smiles: 'CCO'" in text
    # And no affinity property: this call wants coordinates, not a number.
    assert "affinity" not in text


def test_a_cached_msa_is_handed_to_the_cofold():
    """The MSA is keyed on the sequence and shared with the affinity stage, so
    co-folding must not make the public server compute it a second time."""
    text = fold.cofold_input("ACDEFGH", "CCO", {"path": "/tmp/abc.csv"})
    assert "msa: /tmp/abc.csv" in text


def _write_complex(path):
    """A two-residue protein plus a two-atom ligand, built with gemmi rather
    than hand-written mmCIF. The hand-written version parsed to zero models and
    the test failed against correct code, which is a fixture bug wearing the
    costume of a real one."""
    import gemmi
    st = gemmi.Structure()
    st.name = "cofold"
    model = gemmi.Model("1")
    protein = gemmi.Chain("A")
    for i, (name, x) in enumerate((("ALA", 1.0), ("GLY", 3.0)), start=1):
        res = gemmi.Residue()
        res.name, res.seqid = name, gemmi.SeqId(i, " ")
        for atom_name, dx in (("N", 0.0), ("CA", 1.0)):
            atom = gemmi.Atom()
            atom.name, atom.element = atom_name, gemmi.Element("C" if atom_name == "CA" else "N")
            atom.pos = gemmi.Position(x + dx, 1.0, 1.0)
            atom.occ, atom.b_iso = 1.0, 90.0
            res.add_atom(atom)
        protein.add_residue(res)
    lig_chain = gemmi.Chain("B")
    lig = gemmi.Residue()
    lig.name, lig.seqid = "LIG", gemmi.SeqId(1, " ")
    for atom_name, element, x in (("C1", "C", 8.0), ("O1", "O", 9.0)):
        atom = gemmi.Atom()
        atom.name, atom.element = atom_name, gemmi.Element(element)
        atom.pos = gemmi.Position(x, 1.0, 1.0)
        atom.occ, atom.b_iso = 1.0, 50.0
        lig.add_atom(atom)
    lig_chain.add_residue(lig)
    model.add_chain(protein)
    model.add_chain(lig_chain)
    st.add_model(model)
    st.setup_entities()
    st.make_mmcif_document().write_file(str(path))


def test_split_keeps_the_protein_and_sets_the_ligand_aside(tmp_path):
    """One predicted complex in, a protein-only PDB and a ligand SDF out."""
    gemmi = pytest.importorskip("gemmi")
    cif = tmp_path / "c.cif"
    _write_complex(cif)
    protein, ligand = tmp_path / "model_apo.pdb", tmp_path / "lig.sdf"
    kept, atoms = fold.split_cofold(cif, protein, ligand)

    assert kept == 2 and atoms == 2
    body = protein.read_text(encoding="utf-8")
    assert "ATOM" in body and "TER" in body
    # The discarded pose must not survive into the receptor: docking into a
    # pocket that still contains the ligand would find nowhere to put it.
    assert "LIG" not in body
    assert ligand.read_text(encoding="utf-8").rstrip().endswith("$$$$")
