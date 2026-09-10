"""The Example tab: two campaigns, explained to somebody who just arrived.

They are fixed published results, so the runs are held in the route as data.
That is deliberate and worth guarding: a page describing a finished experiment
must not change because a run was deleted or a thirteenth was started.

The page exists to say that the two campaigns disagree, so the assertions below
pin that disagreement rather than any single number. If a future campaign makes
EGFR and beta-2 adrenergic agree, these tests should fail and the prose should
be rewritten, which is the point of pinning it.
"""

from __future__ import annotations

import re

from app import create_app
from app.routes.example import CAMPAIGNS


def campaign(key):
    return next(c for c in CAMPAIGNS if c["key"] == key)


def test_the_page_renders_and_links_every_run():
    with create_app().test_client() as client:
        response = client.get("/example")
        html = response.get_data(as_text=True)
    assert response.status_code == 200
    for c in CAMPAIGNS:
        for run in c["runs"]:
            assert f"/runs/{run['job']}" in html, f"{c['key']} {run['name']}"


def test_it_is_reachable_from_every_page():
    with create_app().test_client() as client:
        for path in ("/", "/prepare", "/analyze", "/runs", "/about", "/example"):
            assert ">Example</a>" in client.get(path).get_data(as_text=True), path


def test_both_campaigns_are_the_same_experiment():
    """Six runs each, two receptors by three modes, or the comparison is not
    like for like and the whole page is invalid."""
    for c in CAMPAIGNS:
        runs = c["runs"]
        assert len(runs) == 6, c["key"]
        assert sorted(r["mode"] for r in runs) == ["dock", "dock", "flex", "flex",
                                                   "hybrid", "hybrid"], c["key"]
        assert len([r for r in runs if r["receptor"] == "co-folded"]) == 3, c["key"]
        assert len([r for r in runs if r["receptor"] == "AlphaFold"]) == 3, c["key"]


def test_the_kinase_result_is_the_one_that_happened():
    runs = campaign("egfr")["runs"]
    by_job = {r["job"]: r for r in runs}
    alphafold = [r for r in runs if r["receptor"] == "AlphaFold"]
    bs = [r for r in runs if r["grade"] == "B"]
    assert len(bs) == 2
    assert all(r["receptor"] == "co-folded" and r["rmsd"] < 2.0 for r in bs)
    # No AlphaFold run came close, which is this campaign's central claim.
    assert all(r["rmsd"] > 4.0 and r["grade"] == "D" for r in alphafold)
    # Flexible docking was the slowest in both arms, and the worst in both.
    for receptor in ("co-folded", "AlphaFold"):
        arm = [r for r in runs if r["receptor"] == receptor]
        assert max(arm, key=lambda r: r["minutes"])["mode"] == "flex", receptor
        assert max(arm, key=lambda r: r["rmsd"])["mode"] == "flex", receptor
    assert by_job["gs_20260909_irczwdpxzhs7"]["score"] == 87.4


def test_the_gpcr_result_inverts_it():
    """The page's reason for existing: the same six choices, opposite answers."""
    runs = campaign("b2ar")["runs"]
    best = min(runs, key=lambda r: r["rmsd"])
    worst = max(runs, key=lambda r: r["rmsd"])
    # Flexible docking on an AlphaFold receptor is the best run of the twelve --
    # the exact combination that came last on the kinase.
    assert best["mode"] == "flex" and best["receptor"] == "AlphaFold"
    assert best["grade"] == "A" and best["rmsd"] < 1.0
    # Plain docking on a co-folded receptor is the only failure -- the exact
    # combination that came first on the kinase.
    assert worst["mode"] == "dock" and worst["receptor"] == "co-folded"
    assert len([r for r in runs if r["grade"] == "D"]) == 1
    # Every AlphaFold run graded B or better, against none on the kinase.
    assert all(r["grade"] in ("A", "B") for r in runs if r["receptor"] == "AlphaFold")


def test_pocket_accuracy_does_not_predict_pose_accuracy():
    """Quoted on the page as 0.016 A of receptor spread giving 1.13 to 5.69 A of
    answer, over the four runs that docked into a receptor held fixed."""
    rigid = [r for r in campaign("b2ar")["runs"] if r["mode"] in ("dock", "hybrid")]
    assert len(rigid) == 4
    cas = [r["pocket_ca"] for r in rigid]
    rmsds = [r["rmsd"] for r in rigid]
    assert max(cas) - min(cas) <= 0.02, cas
    assert min(rmsds) < 1.2 and max(rmsds) > 5.0, rmsds


def test_the_mode_label_claims_match_the_runs():
    """The page says the GPCR labels agreed six times out of six and the kinase
    labels four out of six. Both are load-bearing, in opposite directions."""
    assert sum(r["mode_match"] for r in campaign("b2ar")["runs"]) == 6
    assert sum(r["mode_match"] for r in campaign("egfr")["runs"]) == 4


def test_it_does_not_overclaim():
    """Two targets, two ligands, twelve runs. The page has to say so, and has to
    keep saying that wall-clock time is not a measurement."""
    with create_app().test_client() as client:
        html = client.get("/example").get_data(as_text=True)
    flat = re.sub(r"\s+", " ", html)
    assert "not a benchmark" in flat
    assert "Two targets and two drugs" in flat
    assert "not as a measurement" in flat
    assert "not a membrane simulation" in flat
