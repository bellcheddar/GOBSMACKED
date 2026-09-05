"""Visibility and ownership.

A private run must be absent from listings rather than redacted in them: a
greyed-out row leaks that the run exists and how many there are.
"""

from __future__ import annotations

import pytest

from tests.conftest import register


@pytest.fixture()
def private_job(app, egfr_archive):
    with app.app_context():
        return register(egfr_archive, visibility="private", token="secret-owner-key")


@pytest.fixture()
def public_job(app, beta2_archive):
    with app.app_context():
        return register(beta2_archive, visibility="public", token="other-key")


def test_private_run_is_absent_without_the_token(client, private_job, public_job):
    listed = client.get("/api/runs").get_json()["runs"]
    ids = {run["job_id"] for run in listed}
    assert public_job in ids
    assert private_job not in ids


def test_private_run_appears_with_the_token(client, private_job):
    listed = client.get("/api/runs", headers={"X-Owner-Token": "secret-owner-key"}).get_json()["runs"]
    ids = {run["job_id"] for run in listed}
    assert private_job in ids
    assert next(r for r in listed if r["job_id"] == private_job)["owned"] is True


def test_private_page_shows_only_the_key_prompt(client, private_job):
    response = client.get(f"/runs/{private_job}")
    assert response.status_code == 403
    body = response.get_data(as_text=True)
    assert "This run is private" in body
    # Nothing about the run itself leaks into the prompt page.
    assert "erlotinib" not in body


def test_private_page_opens_with_the_token(client, private_job):
    response = client.get(f"/runs/{private_job}?token=secret-owner-key")
    assert response.status_code == 200
    assert "This run is private" not in response.get_data(as_text=True)


def test_visibility_patch_needs_the_right_token(client, private_job):
    wrong = client.patch(f"/api/runs/{private_job}/visibility",
                         json={"visibility": "public", "token": "not-the-key"})
    assert wrong.status_code == 403

    right = client.patch(f"/api/runs/{private_job}/visibility",
                         json={"visibility": "public", "token": "secret-owner-key"})
    assert right.status_code == 200
    assert right.get_json()["visibility"] == "public"
    assert private_job in {r["job_id"] for r in client.get("/api/runs").get_json()["runs"]}


def test_visibility_value_is_validated(client, private_job):
    response = client.patch(f"/api/runs/{private_job}/visibility",
                            json={"visibility": "unlisted", "token": "secret-owner-key"})
    assert response.status_code == 400


def test_delete_needs_the_job_id_retyped(client, private_job):
    unconfirmed = client.delete(f"/api/runs/{private_job}",
                                json={"token": "secret-owner-key", "confirm": "yes"})
    assert unconfirmed.status_code == 400

    confirmed = client.delete(f"/api/runs/{private_job}",
                              json={"token": "secret-owner-key", "confirm": private_job})
    assert confirmed.status_code == 200
    assert client.get(f"/runs/{private_job}").status_code == 404


def test_downloads_are_guarded_for_private_runs(client, private_job):
    assert client.get(f"/runs/{private_job}/results").status_code == 403
    # 404 rather than 403 once the token is right: the fixture row has no file.
    assert client.get(f"/runs/{private_job}/results?token=secret-owner-key").status_code == 404


def test_owner_token_is_never_stored_in_the_clear(app, private_job):
    from app import db

    with app.app_context():
        row = db.get_job(private_job)
    assert "secret-owner-key" not in str(dict(row))
    assert row["owner_hash"] == db.hash_token("secret-owner-key")


def test_job_ids_are_unguessable():
    from app import db

    ids = {db.new_job_id() for _ in range(200)}
    assert len(ids) == 200
    tail = next(iter(ids)).split("_")[-1]
    assert len(tail) == 12
