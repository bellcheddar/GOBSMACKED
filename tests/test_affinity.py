"""The affinity block, on the server side.

The stage runs on the machine with the GPU; this is what the card makes of what
came back. Three states arrive and all three are normal, so all three are tested:
a block that ran, a block that says why it did not, and no block at all.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from app.services import affinity as aff

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bundle_template"))

from gobsmacked_run import affinity as stage  # noqa: E402

RAN = {
    "requested": True, "ran": True, "route": "boltz-template",
    "unit": "affinity_pred_value is log10(IC50) in micromolar",
    "frames": {"scored": ["pre_md", "post_md_0080", "post_md_0090"]},
    "pre_md": {"affinity_pred_value": -1.41, "pic50": 7.41,
               "affinity_probability_binary": 0.79},
    "post_md": {"affinity_pred_value_mean": -1.58, "affinity_pred_value_sd": 0.19,
                "affinity_pred_value_range": [-1.81, -1.30],
                "pic50_mean": 7.58, "pic50_sd": 0.19,
                "probability_binary_mean": 0.84, "probability_binary_sd": 0.04},
    "delta": {"affinity_pred_value": -0.17, "pic50": 0.17, "probability_binary": 0.05},
}


def test_pic50_is_six_minus_the_model_value():
    """The model reports log10(IC50) in micromolar and lower is stronger; pIC50
    is the negative log10 in molar and higher is stronger. They differ by a sign
    as well as an offset, so a mix-up gives a plausible number, not an obvious
    one."""
    assert stage.pic50(-3.0) == 9.0        # 1 nM
    assert stage.pic50(0.0) == 6.0         # 1 uM
    assert stage.pic50(2.0) == 4.0         # 100 uM


def test_a_run_that_happened_is_two_columns_and_a_change():
    card = aff.summarise(RAN)
    assert card["state"] == "ran"
    labels = [r["label"] for r in card["rows"]]
    assert labels == ["pIC50", "Boltz-2 value", "Binder probability"]
    assert card["n_post_frames"] == 2
    assert card["spread"] == 0.51


def test_better_follows_the_direction_of_each_row():
    """Getting this wrong colours an improvement red: the Boltz-2 value is a
    log10(IC50), so a negative change is an improvement, while pIC50 and the
    probability both read the usual way round."""
    rows = {r["label"]: r for r in aff.summarise(RAN)["rows"]}
    assert rows["pIC50"]["better"] is True             # +0.17, higher is stronger
    assert rows["Boltz-2 value"]["better"] is True     # -0.17, lower is stronger
    assert rows["Binder probability"]["better"] is True


def test_a_worse_pose_after_md_is_marked_as_worse():
    block = {**RAN, "delta": {"pic50": -0.4, "affinity_pred_value": 0.4,
                              "probability_binary": -0.1}}
    rows = {r["label"]: r for r in aff.summarise(block)["rows"]}
    assert rows["pIC50"]["better"] is False
    assert rows["Boltz-2 value"]["better"] is False


def test_a_wide_spread_is_called_out_rather_than_averaged_away():
    block = {**RAN, "post_md": {**RAN["post_md"],
                                "affinity_pred_value_range": [-2.6, -0.4]}}
    assert "not settled" in aff.summarise(block)["verdict"]


@pytest.mark.parametrize("block,state", [
    ({}, "absent"),
    ({"requested": False, "ran": False}, "declined"),
    ({"requested": True, "ran": False, "reason": "no network for the MSA"}, "skipped"),
])
def test_every_way_of_having_no_affinity_is_a_panel_not_an_error(block, state):
    card = aff.summarise(block)
    assert card["state"] == state
    assert card["note"]


def test_the_msa_cache_key_is_the_trimmed_sequence():
    """A campaign may trim to a domain, and the alignment for residues 714-966
    of EGFR is not the alignment for the 1,210-residue precursor. Keying on the
    accession would hand the second to a campaign that asked for the first."""
    assert stage.msa_key("MKVLA") == stage.msa_key("MKVLA")
    assert stage.msa_key("MKVLA") != stage.msa_key("MKVLAG")


def test_frames_come_from_the_tail_of_the_run():
    indices = stage.frame_indices(100, {"n_frames": 5, "window_fraction": 0.2})
    assert len(indices) == 5
    assert min(indices) >= 80 and max(indices) <= 99


def test_a_short_run_still_yields_frames():
    assert stage.frame_indices(3, {"n_frames": 5}) == [2]


def test_the_frame_count_does_not_come_from_a_file_written_later():
    """It used to be read from traj/summary.json, which summarise writes AFTER
    this stage runs. The lookup failed, the cluster fell back to the single
    final frame, and a run that asked for five reported one with a spread of
    zero. Nothing errored; the number was simply less than it claimed."""
    import inspect

    # The invariant, not the prose: this takes a count and touches no
    # filesystem, so there is no file whose absence can quietly shrink it.
    # (Checking the source text instead matched the sentence in the docstring
    # describing the bug, which is the second time that trick has failed here.)
    parameters = list(inspect.signature(stage.frame_indices).parameters)
    assert parameters[0] == "n_frames"
    assert stage.frame_indices(50, {"n_frames": 5, "window_fraction": 0.2}) == [40, 42, 44, 46, 48]


def test_the_stage_declines_rather_than_raises_without_a_campaign(tmp_path):
    out = stage.run({"affinity": {"include": True}}, tmp_path, tmp_path, lambda m: None)
    assert "skipped" in out["headline"]
    block = __import__("json").loads(
        (tmp_path / "affinity" / "affinity.json").read_text(encoding="utf-8"))
    assert block["ran"] is False and block["reason"]


def test_opting_out_is_recorded_rather_than_left_blank(tmp_path):
    out = stage.run({"affinity": {"include": False}}, tmp_path, tmp_path, lambda m: None)
    assert out["headline"] == "not requested"
    block = __import__("json").loads(
        (tmp_path / "affinity" / "affinity.json").read_text(encoding="utf-8"))
    assert block["requested"] is False


# --- the template, which is four traps in a row ------------------------------

@pytest.fixture()
def complex_pdb(tmp_path):
    """Two residues of protein and a ligand, as MD writes them."""
    lines = [
        "ATOM      1  N   ALA A 714      11.104   6.134  -6.504  1.00  0.00           N",
        "ATOM      2  CA  ALA A 714      11.639   6.071  -5.147  1.00  0.00           C",
        "ATOM      3  C   ALA A 714      10.755   5.units 0.00  1.00  0.00           C",
        "ATOM      4  N   GLY A 715      13.049   6.011  -5.108  1.00  0.00           N",
        "ATOM      5  CA  GLY A 715      13.855   5.936  -3.903  1.00  0.00           C",
        "HETATM    6  C1  LIG B   1      15.000   5.000  -3.000  1.00  0.00           C",
        "HETATM    7  C2  LIG B   1      16.000   5.500  -3.500  1.00  0.00           C",
        "END",
    ]
    # The third line above is deliberately malformed; drop it rather than ship a
    # fixture that only parses by luck.
    path = tmp_path / "complex.pdb"
    path.write_text("\n".join(lines[:2] + lines[3:]) + "\n", encoding="utf-8")
    return path


def test_the_template_carries_no_ligand(complex_pdb, tmp_path):
    """A ligand chain in a template walks boltz's parse_polymer off the end of
    its residue list, with nothing in the error naming the file or the reason."""
    out = stage.to_cif(complex_pdb, tmp_path / "t.cif")
    text = out.read_text(encoding="utf-8")
    assert "HETATM" not in text
    assert "C1" not in text.split("_atom_site")[-1]


def test_the_template_carries_entity_poly_seq(complex_pdb, tmp_path):
    """gemmi writes the loop only for an entity whose full_sequence is set, and
    setup_entities() does not set it from coordinates. Boltz builds its residue
    list from exactly that loop."""
    out = stage.to_cif(complex_pdb, tmp_path / "t.cif")
    assert "_entity_poly_seq" in out.read_text(encoding="utf-8")


def test_the_sequence_comes_from_the_structure_not_the_campaign(complex_pdb):
    """The campaign carries UniProt's precursor; the pose is a trimmed domain.
    Boltz maps template onto sequence by index, so the two must agree."""
    assert stage.sequence_of(complex_pdb) == "AG"


def test_the_template_is_forced_with_a_threshold(complex_pdb, tmp_path):
    """An unforced template only conditions the trunk. Forcing it is what makes
    the affinity head read the conformation MD produced rather than one the
    model invented, and Boltz refuses `force` without a threshold."""
    cif = stage.to_cif(complex_pdb, tmp_path / "t.cif")
    text = stage.boltz_input("AG", "CCO", cif, {"cached": False}, 2.0)
    assert "force: true" in text
    assert "threshold: 2.0" in text
    assert "templates:" in text and str(cif) in text


def test_potentials_are_requested_because_force_needs_them(tmp_path):
    """Boltz enforces a forced template through a steering potential, so
    without the flag `force: true` is silently inert."""
    cmd = stage.boltz_command(Path("in.yaml"), tmp_path, {}, {"cached": True})
    assert "--use_potentials" in cmd
    # A cached MSA means the server is not asked; an uncached one means it is.
    assert "--use_msa_server" not in cmd
    assert "--use_msa_server" in stage.boltz_command(
        Path("in.yaml"), tmp_path, {}, {"cached": False})


def test_a_cached_msa_is_offered_to_boltz(complex_pdb, tmp_path):
    cif = stage.to_cif(complex_pdb, tmp_path / "t.cif")
    text = stage.boltz_input("AG", "CCO", cif, {"cached": True, "path": "/tmp/m.csv"})
    assert "msa: /tmp/m.csv" in text
    assert "msa:" not in stage.boltz_input("AG", "CCO", cif, {"cached": False})


def test_the_first_pose_pays_for_the_msa_and_the_rest_do_not(tmp_path, monkeypatch):
    """Six poses meant six identical queries to a free public server for the
    same protein, which is what the spec said never to do and what the first
    implementation did anyway."""
    monkeypatch.setattr(stage, "MSA_CACHE", tmp_path / "cache")
    work = tmp_path / "boltz" / "pre_md" / "boltz_results_pre_md" / "msa"
    work.mkdir(parents=True)
    (work / "pre_md_0.csv").write_text("key,sequence\n-1,AG\n1,AC\n", encoding="utf-8")

    before = {"cached": False, "key": "abc123", "path": None, "depth": None}
    after = stage.capture_msa(tmp_path / "boltz" / "pre_md", before, lambda m: None)
    assert after["cached"] is True
    assert after["depth"] == 2
    assert Path(after["path"]).read_text(encoding="utf-8").startswith("key,sequence")

    # And a second call is a no-op rather than a second copy.
    assert stage.capture_msa(tmp_path / "nowhere", after, lambda m: None) == after


def test_a_run_with_no_msa_to_capture_carries_on(tmp_path, monkeypatch):
    monkeypatch.setattr(stage, "MSA_CACHE", tmp_path / "cache")
    before = {"cached": False, "key": "abc123", "path": None, "depth": None}
    assert stage.capture_msa(tmp_path / "empty", before, lambda m: None) == before


def test_msa_depth_counts_both_formats(tmp_path):
    a3m = tmp_path / "m.a3m"
    a3m.write_text(">one\nAG\n>two\nAC\n", encoding="utf-8")
    csv = tmp_path / "m.csv"
    csv.write_text("key,sequence\n-1,AG\n1,AC\n", encoding="utf-8")
    assert stage.depth_of(a3m) == 2
    assert stage.depth_of(csv) == 2


def test_the_panel_describes_the_route_it_actually_ran():
    """On the forced route diffusion still runs, constrained to the pose; on the
    unforced one it runs free. Describing one while reporting the other is a
    lie in the only sentence that says what the number means."""
    forced = aff.summarise({**RAN, "route": "boltz-forced-template",
                            "engine": {"template_threshold_a": 2.0}})
    assert "held within 2 A" in forced["how"]

    unforced = aff.summarise({**RAN, "route": "boltz-template", "engine": {}})
    assert "did not constrain it" in unforced["how"]
    assert "partly" in unforced["how"]
