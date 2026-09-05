"""Shared fixtures.

Every test runs against a temporary data directory and a temporary database, so
a test run never touches the copy of the app someone is using, and two runs
cannot see each other's rows.
"""

from __future__ import annotations

import os
import shutil
import sys
import tarfile
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _online() -> bool:
    """Is the network up? Several tests fetch a reference structure from RCSB."""
    import socket

    try:
        socket.create_connection(("files.rcsb.org", 443), timeout=5).close()
        return True
    except OSError:
        return False


needs_network = pytest.mark.skipif(not _online(), reason="no network")


@pytest.fixture()
def app(tmp_path, monkeypatch):
    from app import config

    monkeypatch.setattr(config, "DATA_DIR", tmp_path / "data")
    monkeypatch.setattr(config, "RUNS_DIR", tmp_path / "data" / "runs")
    monkeypatch.setattr(config, "STRUCT_CACHE", tmp_path / "data" / "structures")
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "gobsmacked.db")

    from app import create_app

    application = create_app({"TESTING": True})
    yield application


@pytest.fixture()
def client(app):
    return app.test_client()


@pytest.fixture()
def egfr_archive() -> Path:
    path = FIXTURES / "egfr_results.tar.gz"
    if not path.exists():
        pytest.skip("run tests/fixtures/build_fixtures.py first")
    return path


@pytest.fixture()
def beta2_archive() -> Path:
    path = FIXTURES / "beta2_results.tar.gz"
    if not path.exists():
        pytest.skip("run tests/fixtures/build_fixtures.py first")
    return path


def campaign_of(archive: Path) -> dict:
    with tarfile.open(archive) as tar:
        return yaml.safe_load(tar.extractfile("results/campaign.yaml").read())


def register(archive: Path, visibility: str = "public", token: str = "fixture-token") -> str:
    """Insert the jobs row a results upload needs, and return the job ID."""
    from app import db

    campaign = campaign_of(archive)
    job_id = campaign["job_id"]
    db.insert_job(
        job_id=job_id, title=campaign.get("title", ""),
        uniprot=campaign["protein"]["uniprot"], ligand_name=campaign["ligand"]["name"],
        ligand_smiles=campaign["ligand"]["smiles"], family=campaign["protein"]["family"],
        reference_pdb=campaign["reference"]["pdb_id"], status="prepared",
        visibility=visibility, owner_hash=db.hash_token(token),
        campaign_yaml=yaml.safe_dump(campaign, sort_keys=False),
    )
    return job_id
