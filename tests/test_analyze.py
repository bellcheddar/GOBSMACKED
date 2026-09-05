"""The whole Analyze pipeline, on two archives built from real crystals.

Slow (PLIP three times, PandaMap twice, a trajectory read) and worth it: this is
the only test that proves the parts fit together, and it is what the two
acceptance criteria in the build spec actually say.
"""

from __future__ import annotations

import io
import json

import pytest

from tests.conftest import needs_network, register

pytestmark = needs_network


def upload(client, archive):
    return client.post("/api/upload",
                       data={"file": (io.BytesIO(archive.read_bytes()), "results.tar.gz")},
                       content_type="multipart/form-data")


@pytest.fixture()
def egfr_card(app, client, egfr_archive):
    from app import db

    with app.app_context():
        job_id = register(egfr_archive)
    response = upload(client, egfr_archive)
    assert response.status_code == 200, response.get_json()
    with app.app_context():
        return job_id, json.loads(db.get_job(job_id)["scorecard_json"])


@pytest.fixture()
def beta2_card(app, client, beta2_archive):
    from app import db

    with app.app_context():
        job_id = register(beta2_archive)
    response = upload(client, beta2_archive)
    assert response.status_code == 200, response.get_json()
    with app.app_context():
        return job_id, json.loads(db.get_job(job_id)["scorecard_json"])


def test_egfr_scores_against_1m17_and_is_type_one(egfr_card):
    _, card = egfr_card
    assert card["reference"]["pdb_id"] == "1M17"
    assert card["scorecard"]["verified"] is True
    assert card["scorecard"]["score"] >= 70
    modes = card["modes"]
    assert modes["family"] == "kinase"
    assert modes["predicted"]["label"] == "I"
    assert modes["reference"]["label"] == "I"
    assert modes["predicted"]["dfg"] == "in" and modes["reference"]["dfg"] == "in"
    assert modes["verdict"]["match"] is True


def test_beta2_scores_against_2rh1_and_is_orthosteric_inactive(beta2_card):
    _, card = beta2_card
    assert card["reference"]["pdb_id"] == "2RH1"
    modes = card["modes"]
    assert modes["family"] == "gpcr"
    assert modes["predicted"]["label"] == "orthosteric"
    assert modes["reference"]["label"] == "orthosteric"
    assert modes["predicted"]["state"] == "inactive-like"
    assert modes["reference"]["state"] == "inactive-like"
    assert modes["verdict"]["match"] is True


def test_geometry_is_measured_on_the_pocket(egfr_card):
    _, card = egfr_card
    geometry = card["geometry"]["md_final"]
    assert geometry["pocket_ca_atoms"] >= 20
    assert 0 < geometry["pocket_ca_rmsd"] < 5
    assert geometry["ligand_rmsd"] is not None
    # The whole-chain TM-score is reported alongside, not used for the grade.
    assert 0.5 < geometry["tm_score"] <= 1.0


def test_interaction_fingerprints_are_compared_in_reference_numbering(egfr_card):
    _, card = egfr_card
    assert card["jaccard"]["md_final"] is not None
    assert card["interaction_table"]
    # 1M17 numbers 24 lower than UniProt; without the mapping every contact
    # would look lost and the Jaccard would be zero.
    assert card["jaccard"]["md_final"] > 0.2


def test_dynamics_panels_have_data(egfr_card):
    _, card = egfr_card
    dynamics = card["dynamics"]
    assert dynamics["frames"] == 20
    assert len(dynamics["ligand_rmsd_pose1"]) == 20
    assert dynamics["drift"] is not None
    assert dynamics["contacts"]["residues"]
    assert dynamics["ligand_rmsd_reference"]["series"]


def test_results_page_and_report_render(client, egfr_card):
    job_id, _ = egfr_card
    page = client.get(f"/runs/{job_id}")
    assert page.status_code == 200
    body = page.get_data(as_text=True)
    for section in ("Scorecard", "Complex", "Overlay", "Dynamics", "Mode", "Report"):
        assert section in body
    report = client.get(f"/runs/{job_id}/report")
    assert report.status_code == 200
    assert "GOBSMACK score" in report.get_data(as_text=True) or "Unverified" in report.get_data(as_text=True)


def test_an_archive_for_an_unknown_job_is_refused(client, egfr_archive):
    response = upload(client, egfr_archive)
    assert response.status_code == 404
    assert "no record" in response.get_json()["error"]


def test_upload_rejects_something_that_is_not_an_archive(client):
    response = client.post("/api/upload",
                           data={"file": (io.BytesIO(b"not a tarball"), "notes.txt")},
                           content_type="multipart/form-data")
    assert response.status_code == 400


def test_the_ligand_is_named_per_file(egfr_card):
    """PandaDock names the ligand in the pose complex and OpenMM names it in the
    relaxed one. Filtering PLIP on the wrong name returns an empty site rather
    than an error, which reads as a pose that makes no interactions at all."""
    _, card = egfr_card
    names = card["ligand_resnames"]
    assert set(names) >= {"pose1", "md_final"}
    for state, name in names.items():
        assert name, f"no ligand found in {state}"
    # Both states produced a fingerprint, which is what the per-file name buys.
    assert card["jaccard"]["pose1"] is not None
    assert card["jaccard"]["md_final"] is not None
