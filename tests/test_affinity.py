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


def test_frames_come_from_the_tail_of_the_run(tmp_path):
    (tmp_path / "traj").mkdir()
    (tmp_path / "traj" / "summary.json").write_text('{"frames": 100}', encoding="utf-8")
    indices = stage.frame_indices(tmp_path, {"n_frames": 5, "window_fraction": 0.2})
    assert len(indices) == 5
    assert min(indices) >= 80 and max(indices) <= 99


def test_a_short_run_still_yields_frames(tmp_path):
    (tmp_path / "traj").mkdir()
    (tmp_path / "traj" / "summary.json").write_text('{"frames": 3}', encoding="utf-8")
    assert stage.frame_indices(tmp_path, {"n_frames": 5}) == [2]


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
