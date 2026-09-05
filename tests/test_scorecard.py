"""The grades, the composite and the two things that must not drift."""

from __future__ import annotations

import pytest

from app.services import scorecard

PERFECT = {
    "ligand_rmsd": 0.0, "plip_jaccard": 1.0, "pocket_ca_rmsd": 0.0,
    "chi1_agreement": 1.0, "md_drift": 0.0, "rescue": 0.6,
}
ALL_VALID = {"no_clash": True, "bond_lengths": True, "chirality": True, "inside_box": True}


def test_structure_against_itself_scores_100():
    """The only self-evidently correct anchor the score has."""
    card = scorecard.composite(PERFECT, ALL_VALID)
    assert card["score"] == 100.0
    assert card["grade"] == "A"
    assert all(m["grade"] == "A" for m in card["metrics"])


def test_five_angstrom_pose_fails_the_ligand_metric():
    card = scorecard.composite({**PERFECT, "ligand_rmsd": 5.0}, ALL_VALID)
    ligand = next(m for m in card["metrics"] if m["key"] == "ligand_rmsd")
    assert ligand["grade"] == "F"
    assert card["score"] < 100


@pytest.mark.parametrize("value,expected", [
    (0.9, "A"), (1.0, "A"), (1.5, "B"), (2.0, "B"), (2.5, "C"), (3.0, "C"),
    (3.5, "D"), (4.0, "D"), (4.1, "F"),
])
def test_ligand_rmsd_boundaries(value, expected):
    """Thresholds are inclusive at the top of each band, as the About table says."""
    metric = next(m for m in scorecard.METRICS if m.key == "ligand_rmsd")
    assert metric.grade(value) == expected


@pytest.mark.parametrize("value,expected", [
    (0.80, "A"), (0.75, "A"), (0.60, "B"), (0.55, "B"), (0.45, "C"),
    (0.30, "D"), (0.20, "F"),
])
def test_jaccard_boundaries_run_the_other_way(value, expected):
    metric = next(m for m in scorecard.METRICS if m.key == "plip_jaccard")
    assert metric.grade(value) == expected


def test_validity_failure_caps_the_composite():
    card = scorecard.composite(PERFECT, {**ALL_VALID, "no_clash": False})
    assert card["validity"]["pass"] is False
    assert card["validity"]["capped"] is True
    assert card["score"] == scorecard.VALIDITY_FAIL_CAP


def test_a_low_score_is_not_raised_by_the_cap():
    poor = {"ligand_rmsd": 6.0, "plip_jaccard": 0.05, "pocket_ca_rmsd": 4.0,
            "chi1_agreement": 0.1, "md_drift": 4.0, "rescue": -2.0}
    card = scorecard.composite(poor, {**ALL_VALID, "no_clash": False})
    assert card["score"] < scorecard.VALIDITY_FAIL_CAP
    assert card["validity"]["capped"] is False


def test_no_reference_means_no_composite():
    """A run nobody verified must not come back with an A for the two metrics
    that survive without a reference."""
    card = scorecard.composite({"md_drift": 0.3}, ALL_VALID)
    assert card["score"] is None
    assert card["grade"] is None
    assert card["verified"] is False
    assert "Unverified" in card["label"]
    assert "Ligand RMSD" in card["unmeasured"]


def test_missing_metrics_renormalise_rather_than_score_zero():
    without_dynamics = {k: v for k, v in PERFECT.items() if k not in ("md_drift", "rescue")}
    card = scorecard.composite(without_dynamics, ALL_VALID)
    assert card["score"] == 100.0
    assert set(card["unmeasured"]) == {"Drift, last 200 ps", "MD rescue"}


def test_every_metric_carries_a_sentence():
    card = scorecard.composite(PERFECT, ALL_VALID)
    for metric in card["metrics"]:
        assert metric["note"], f"{metric['label']} has no explanation"
        assert len(metric["note"]) > 20


def test_weights_sum_to_one_hundred():
    total = sum(m.weight for m in scorecard.METRICS) + scorecard.VALIDITY_WEIGHT
    assert total == 100
