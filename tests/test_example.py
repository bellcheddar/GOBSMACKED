"""The Example tab: one campaign, explained to somebody who just arrived.

It is a fixed published result, so the runs are held in the route as data. That
is deliberate and worth guarding: a page describing a finished experiment must
not change because a run was deleted or a seventh was started.
"""

from __future__ import annotations

from app import create_app
from app.routes.example import RUNS


def test_the_page_renders_and_links_every_run():
    with create_app().test_client() as client:
        response = client.get("/example")
        html = response.get_data(as_text=True)
    assert response.status_code == 200
    for run in RUNS:
        assert f"/runs/{run['job']}" in html, run["name"]


def test_it_is_reachable_from_every_page():
    with create_app().test_client() as client:
        for path in ("/", "/prepare", "/analyze", "/runs", "/about", "/example"):
            assert ">Example</a>" in client.get(path).get_data(as_text=True), path


def test_the_result_it_describes_is_the_one_that_happened():
    """The numbers on the page are the numbers from those runs. Pinned, because
    prose drifts from data silently and a worked example that misreports itself
    is worse than none."""
    by_job = {r["job"]: r for r in RUNS}
    cofolded = [r for r in RUNS if r["receptor"] == "co-folded"]
    alphafold = [r for r in RUNS if r["receptor"] == "AlphaFold"]
    assert len(RUNS) == 6 and len(cofolded) == 3 and len(alphafold) == 3
    # Two B grades, both co-folded, both under 2 A.
    bs = [r for r in RUNS if r["grade"] == "B"]
    assert len(bs) == 2
    assert all(r["receptor"] == "co-folded" and r["rmsd"] < 2.0 for r in bs)
    # No AlphaFold run came close, which is the page's central claim.
    assert all(r["rmsd"] > 4.0 and r["grade"] == "D" for r in alphafold)
    # Flexible docking was the slowest in both arms, which is the other claim.
    for arm in (cofolded, alphafold):
        slowest = max(arm, key=lambda r: r["minutes"])
        assert slowest["mode"] == "flex", arm
    assert by_job["gs_20260909_irczwdpxzhs7"]["score"] == 87.4


def test_it_does_not_overclaim():
    """One target, one ligand, six runs. The page has to say so."""
    import re
    with create_app().test_client() as client:
        html = client.get("/example").get_data(as_text=True)
    # Whitespace-normalised: the assertion is about what the page says, and a
    # line break landing mid-sentence is not a change of meaning.
    flat = re.sub(r"\s+", " ", html)
    assert "not a benchmark" in flat
    assert "One target and one drug" in flat
