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


# --- what the CLI says it is about to do -------------------------------------

def _runner():
    """run.py is a script, not a package member, so it is loaded by path."""
    import importlib.util
    path = Path(__file__).resolve().parents[1] / "bundle_template" / "run.py"
    spec = importlib.util.spec_from_file_location("gobsmacked_runner", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["gobsmacked_runner"] = module
    spec.loader.exec_module(module)
    return module


CAMPAIGN = {"md": {"production_ps": 500, "equilibration_ps": 100, "frame_interval_ps": 10},
            "docking": {"mode": "dock"}, "affinity": {"include": True, "n_frames": 5}}


def test_the_plan_names_the_method_it_will_actually_use():
    """Stage 1 said "ESMFold" even when the campaign asked for co-folding."""
    runner = _runner()
    esm = runner.blurb_for("fold", {**CAMPAIGN, "fold": {"method": "esmfold"}})
    boltz = runner.blurb_for("fold", {**CAMPAIGN, "fold": {"method": "boltz2"}})
    assert "ESMFold" in esm
    assert "ESMFold" not in boltz
    assert "co-fold" in boltz and "Boltz-2" in boltz
    # And the other stages are untouched by the campaign.
    assert runner.blurb_for("md", CAMPAIGN) == runner.BLURB["md"]


def test_co_folding_is_not_quoted_at_esmfolds_three_minutes():
    """It quoted ~3 min for a stage that takes about a quarter of an hour."""
    runner = _runner()
    esm = runner.estimate_seconds("fold", {**CAMPAIGN, "fold": {"method": "esmfold"}}, False)
    boltz = runner.estimate_seconds("fold", {**CAMPAIGN, "fold": {"method": "boltz2"}}, False)
    assert esm == pytest.approx(180, abs=1)
    assert boltz >= 600, "co-folding is minutes, not one ESMFold pass"
    assert boltz > 4 * esm


def test_a_supplied_model_still_costs_nothing_either_way():
    """The usual case: the server shipped a structure, so fold is skipped."""
    runner = _runner()
    for method in ("esmfold", "boltz2"):
        campaign = {**CAMPAIGN, "fold": {"method": method}}
        assert runner.estimate_seconds("fold", campaign, skip_fold=True) == 0.0


# --- the frame the box lives in ----------------------------------------------
# The first co-folded run put its ligand 87 A from the site. The box centre had
# been computed at Prepare time from the fetched structure; Boltz-2 returns its
# own origin, 59 A away, and 150 of 9,391 receptor atoms fell inside the box.

def test_the_box_centre_comes_from_the_co_folded_ligand(tmp_path):
    """Not from the pocket residue list, which does not transfer: a campaign can
    carry a pocket in a crystal's numbering while the co-folded model is
    numbered from the sequence, and those name different residues."""
    sdf = tmp_path / "cofold_ligand.sdf"
    sdf.write_text(
        "x\n  GOBSMACKED\n\n  2  0  0  0  0  0  0  0  0  0999 V2000\n"
        "   10.0000    0.0000    0.0000 C   0  0\n"
        "   20.0000   10.0000    4.0000 O   0  0\n"
        "M  END\n$$$$\n", encoding="utf-8")
    assert fold._ligand_centre(sdf) == [15.0, 5.0, 2.0]


def test_a_missing_co_folded_ligand_is_detected_not_assumed(tmp_path):
    """Silently keeping the old centre is the failure being fixed, so the
    absence has to be visible to the caller."""
    assert fold._ligand_centre(tmp_path / "nothing.sdf") is None
    assert fold._distance([0, 0, 0], [3, 4, 0]) == 5.0
    assert fold._distance(None, [1, 2, 3]) is None


# --- co-fold the domain, not the precursor -----------------------------------

def test_only_the_docked_domain_is_co_folded():
    """It handed Boltz all 1,210 residues of EGFR when the campaign asked for
    253, then prep threw away four fifths of the prediction. 54 minutes against
    an estimate of 15."""
    seq = "".join("ACDEFGHIKL"[i % 10] for i in range(1210))
    campaign = {"protein": {"residue_range": [714, 966]}}
    domain, offset = fold._domain_of(campaign, seq, lambda m: None)
    assert len(domain) == 253
    assert offset == 714
    assert domain == seq[713:966]


def test_no_range_means_the_whole_sequence():
    seq = "ACDEFGHIKL"
    for campaign in ({}, {"protein": {}}, {"protein": {"residue_range": None}}):
        assert fold._domain_of(campaign, seq, lambda m: None) == (seq, 1)


def test_a_range_that_does_not_fit_falls_back_loudly():
    """A campaign can be internally inconsistent -- this one carried a pocket in
    one numbering and a range in another. Folding a nonsense slice silently
    would be worse than folding all of it."""
    seq = "ACDEFGHIKL"
    said = []
    got = fold._domain_of({"protein": {"residue_range": [5, 900]}}, seq, said.append)
    assert got == (seq, 1)
    assert said and "does not fit" in said[0]


def test_the_split_renumbers_into_the_campaigns_frame(tmp_path):
    """Boltz numbers what it was given from 1. The campaign calls the same
    residues 714 onward, and prep's trim reads those numbers."""
    pytest.importorskip("gemmi")
    cif = tmp_path / "c.cif"
    _write_complex(cif)
    protein = tmp_path / "m.pdb"
    fold.split_cofold(cif, protein, tmp_path / "l.sdf", first_residue=714)
    numbers = {int(line[22:26]) for line in protein.read_text().splitlines()
               if line.startswith("ATOM")}
    assert numbers == {714, 715}, numbers


def test_the_counts_line_is_not_read_as_an_atom(tmp_path):
    """`  2  0  0  0 ...` parses as the coordinate (2, 0, 0) if the atom block
    is found by looking for lines with three numbers. With two real atoms it
    moved the centre by 5 A."""
    sdf = tmp_path / "l.sdf"
    sdf.write_text(
        "x\n  GOBSMACKED\n\n  2  0  0  0  0  0  0  0  0  0999 V2000\n"
        "   10.0000    0.0000    0.0000 C   0  0\n"
        "   20.0000   10.0000    4.0000 O   0  0\n"
        "M  END\n$$$$\n", encoding="utf-8")
    assert fold._ligand_centre(sdf) == [15.0, 5.0, 2.0]


def test_a_truncated_atom_block_is_refused_rather_than_averaged(tmp_path):
    """A short read would give a centre that looks plausible and is not."""
    sdf = tmp_path / "l.sdf"
    sdf.write_text(
        "x\n  GOBSMACKED\n\n  3  0  0  0  0  0  0  0  0  0999 V2000\n"
        "   10.0000    0.0000    0.0000 C   0  0\n"
        "M  END\n$$$$\n", encoding="utf-8")
    assert fold._ligand_centre(sdf) is None


# --- the archive has to record the box that was actually used ----------------

def test_a_recentred_box_is_written_into_the_archived_campaign():
    """The archived campaign is copied verbatim before the stages run, so a
    correction made during the run never reached it. Co-folding moves the box
    36 A into its own frame; the server then checked "is the ligand inside the
    docking box" against a centre nowhere near the pose, failed it, and capped a
    correct run at 40."""
    runner = _runner()
    before = {"center": [-9.87, 33.08, 14.14], "box": [31.7, 20.7, 23.1]}
    after = {"center": [8.25, 2.23, 8.45], "box": [31.7, 20.7, 23.1]}
    edits = runner._campaign_edits(before, after)
    assert len(edits) == 1
    assert "re-centred by 36" in edits[0]
    assert "the one that ran" in edits[0]


def test_an_unmoved_box_is_not_reported_as_an_edit():
    """Every non-co-folded run goes through this. A note on all of them would be
    noise, and rounding is not a decision."""
    runner = _runner()
    same = {"center": [1.0, 2.0, 3.0]}
    assert runner._campaign_edits(same, same) == []
    assert runner._campaign_edits(same, {"center": [1.0, 2.0, 3.2]}) == []
    assert runner._campaign_edits({}, {}) == []
    assert runner._campaign_edits({"center": None}, {"center": [1, 2, 3]}) == []


def test_the_box_correction_survives_a_resume(tmp_path):
    """The re-centring lived in fold.cofold, and a resume skips fold because
    fold.done exists. A co-folded run that failed in dock and retried came back
    with the campaign's ORIGINAL centre and docked 36 A from its own receptor:
    the retry silently undid the first attempt's fix while looking normal.

    Derived from disk now, so it does not matter which stages are being rerun.
    """
    runner = _runner()
    results = tmp_path / "results"
    (results / "cofold").mkdir(parents=True)
    (results / "cofold" / "cofold_ligand.sdf").write_text(
        "x\n  GOBSMACKED\n\n  2  0  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0\n"
        "   20.0000   10.0000    4.0000 O   0  0\n"
        "M  END\n$$$$\n", encoding="utf-8")

    class Log:
        def __init__(self): self.said = []
        def detail(self, m, **k): self.said.append(m)

    log = Log()
    campaign = {"pocket": {"center": [-9.87, 33.08, 14.14], "box": [30, 20, 23]}}
    runner._apply_cofold_box(campaign, results, log)
    assert campaign["pocket"]["center"] == [10.0, 5.0, 2.0]
    assert campaign["pocket"]["box"] == [30, 20, 23], "only the centre moves"
    assert log.said and "co-folded ligand" in log.said[0]


def test_a_run_that_did_not_co_fold_is_left_alone(tmp_path):
    """Every supplied-structure and ESMFold run goes through this."""
    runner = _runner()
    results = tmp_path / "results"
    results.mkdir()
    campaign = {"pocket": {"center": [1.0, 2.0, 3.0]}}

    class Log:
        def detail(self, m, **k): raise AssertionError("should not have spoken")

    runner._apply_cofold_box(campaign, results, Log())
    assert campaign["pocket"]["center"] == [1.0, 2.0, 3.0]


def test_applying_it_twice_changes_nothing(tmp_path):
    """It runs on every invocation, so it has to be idempotent or a second
    resume would report a move of zero as though it were news."""
    runner = _runner()
    results = tmp_path / "results"
    (results / "cofold").mkdir(parents=True)
    (results / "cofold" / "cofold_ligand.sdf").write_text(
        "x\n  GOBSMACKED\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "    5.0000    5.0000    5.0000 C   0  0\nM  END\n$$$$\n", encoding="utf-8")
    said = []

    class Log:
        def detail(self, m, **k): said.append(m)

    campaign = {"pocket": {"center": [0.0, 0.0, 0.0]}}
    runner._apply_cofold_box(campaign, results, Log())
    runner._apply_cofold_box(campaign, results, Log())
    assert campaign["pocket"]["center"] == [5.0, 5.0, 5.0]
    assert len(said) == 1, "the second pass should be a no-op"
