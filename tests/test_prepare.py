"""Prepare: what the panels post, and what the bundle comes out as."""

from __future__ import annotations

import tarfile

import pytest
import yaml

from app.services import bundle
from tests.conftest import needs_network

EGFR_SMILES = "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"


def test_residue_parsing_defaults_the_chain():
    assert bundle.parse_residues(["A:718", "745", "B:790"]) == [("A", 718), ("A", 745), ("B", 790)]
    assert bundle.parse_residues(["nonsense", ""]) == []


def test_ligand_endpoint_validates_and_draws(client):
    ok = client.post("/api/ligand", json={"smiles": EGFR_SMILES})
    assert ok.status_code == 200
    body = ok.get_json()
    assert body["formula"] == "C22H23N3O4"
    assert body["svg"].startswith("<?xml") or "<svg" in body["svg"]

    bad = client.post("/api/ligand", json={"smiles": "this is not a molecule"})
    assert bad.status_code == 400
    assert "RDKit" in bad.get_json()["error"]


def test_bundle_needs_a_pocket(client):
    response = client.post("/api/bundle", json={
        "protein": {"sequence": "MKV", "chain": "A"},
        "ligand": {"smiles": EGFR_SMILES, "name": "erlotinib"},
        "pocket": {},
    })
    assert response.status_code == 400
    assert "pocket" in response.get_json()["error"].lower()


def test_bundle_is_written_and_is_complete(client, app):
    from app import config, db

    response = client.post("/api/bundle", json={
        "title": "test run",
        "visibility": "private",
        "protein": {"uniprot": "P00533", "sequence": "MRPSGTAGAALLALLAALCPASRA",
                    "chain": "A", "family": "kinase", "source_structure": "afdb",
                    "source_id": "AF-P00533-F1"},
        "ligand": {"name": "erlotinib", "smiles": EGFR_SMILES, "protonation_ph": 7.4},
        "pocket": {"method": "residues", "residues": ["A:745"],
                   "center": [1.0, 2.0, 3.0], "box": [22, 22, 22]},
        "reference": {"pdb_id": "1M17", "ligand_ccd": "AQ4", "chain": "A"},
        "docking": {"mode": "hybrid", "num_poses": 10},
        "md": {"production_ps": 1000, "frame_interval_ps": 10},
    })
    assert response.status_code == 200, response.get_json()
    body = response.get_json()
    job_id = body["job_id"]

    # The owner token is shown once and stored only as a hash.
    with app.app_context():
        row = db.get_job(job_id)
    assert row["visibility"] == "private"
    assert row["owner_hash"] == db.hash_token(body["owner_token"])
    assert body["owner_token"] not in (row["campaign_yaml"] or "")[:0] + str(row["owner_hash"])

    archive = config.RUNS_DIR / job_id / f"run_bundle_{job_id}.tar.gz"
    assert archive.exists()
    with tarfile.open(archive) as tar:
        names = {n.split("/", 1)[-1] for n in tar.getnames()}
        campaign = yaml.safe_load(tar.extractfile(f"run_bundle_{job_id}/campaign.yaml").read())
    for needed in ("run.py", "run.sh", "pixi.toml", "pixi.lock", "campaign.yaml", "README.md",
                   "gobsmacked_run/fold.py", "gobsmacked_run/prep.py",
                   "gobsmacked_run/dock.py", "gobsmacked_run/md.py",
                   "gobsmacked_run/summarise.py", "gobsmacked_run/schema.py"):
        assert needed in names, f"{needed} missing from the bundle"

    assert campaign["job_id"] == job_id
    assert campaign["protein"]["family"] == "kinase"
    assert campaign["pocket"]["center"] == [1.0, 2.0, 3.0]
    # The token rides inside the campaign so uploading results to a private run
    # needs no typing.
    assert campaign["owner_token"] == body["owner_token"]


def test_bundle_download_is_guarded(client, app):
    response = client.post("/api/bundle", json={
        "visibility": "private",
        "protein": {"sequence": "MRPSGTAGAALLALL", "chain": "A"},
        "ligand": {"name": "x", "smiles": "CCO"},
        "pocket": {"residues": ["A:1"], "center": [0, 0, 0], "box": [18, 18, 18]},
    })
    job_id = response.get_json()["job_id"]
    token = response.get_json()["owner_token"]
    assert client.get(f"/runs/{job_id}/bundle").status_code == 403
    assert client.get(f"/runs/{job_id}/bundle?token={token}").status_code == 200


@needs_network
def test_fetch_resolves_an_accession_to_a_model(client):
    response = client.post("/api/fetch", json={"query": "P00533"})
    assert response.status_code == 200
    body = response.get_json()
    assert body["accession"] == "P00533"
    assert body["source_structure"] in ("afdb", "pdb", "esm_atlas")
    assert len(body["sequence"]) == 1210
    # The sequence track counts from 1 and the model may not, so the mapping
    # comes back with the structure.
    assert body["numbering"]["745"]


@needs_network
def test_annotation_routes_a_kinase_to_klifs(client):
    from app.services import fetch as fetch_svc

    entry = fetch_svc.fetch_uniprot("P00533")
    response = client.post("/api/annotate", json={
        "uniprot": "P00533", "sequence": entry["sequence"], "gene": "EGFR",
        "features": entry["features"]})
    body = response.get_json()
    assert body["family"] == "kinase"
    assert body["klifs"]["named_positions"]["gatekeeper"] == 790
    assert any(p["label"] == "gatekeeper" for p in body["positions"])


@needs_network
def test_reference_search_prefers_the_same_ligand(client):
    response = client.post("/api/references", json={"uniprot": "P00533", "smiles": EGFR_SMILES})
    body = response.get_json()
    assert body["default"] == "1M17"
    top = body["entries"][0]
    assert top["best_ligand"]["ccd"] == "AQ4"
    assert top["tanimoto"] == 1.0


@needs_network
def test_reference_site_is_renumbered_onto_the_model(client):
    """1M17 numbers EGFR from the mature protein, 24 lower than UniProt. Pasting
    its site straight into the pocket picker would select real residues that are
    the wrong ones, which is worse than an error."""
    fetched = client.post("/api/fetch", json={"query": "P00533"}).get_json()
    response = client.post("/api/reference_site", json={
        "pdb_id": "1M17", "ligand_ccd": "AQ4",
        "structure_name": fetched["structure_name"], "chain": fetched["chain"]})
    body = response.get_json()
    assert body["renumbered"] is True
    assert body["reference_residues"], "no site found in the crystal"

    def numbers(labels):
        return {int(label.split(":")[-1]) for label in labels}

    crystal = numbers(body["reference_residues"])
    model = numbers(body["residues"])
    # The gatekeeper is Thr766 in 1M17 and Thr790 in UniProt numbering.
    assert 766 in crystal and 790 in model
    offsets = {m - c for m, c in zip(sorted(model), sorted(crystal))}
    assert offsets == {24}, f"expected a uniform +24 shift, got {sorted(offsets)}"


@needs_network
def test_box_is_sized_from_the_ligand_not_the_residue_shell(client):
    """Sizing a docking box from the residues around a ligand makes it five
    times larger than the ligand needs, which samples the true site less
    densely for no benefit."""
    fetched = client.post("/api/fetch", json={"query": "P00533"}).get_json()
    site = client.post("/api/reference_site", json={
        "pdb_id": "1M17", "ligand_ccd": "AQ4",
        "structure_name": fetched["structure_name"], "chain": fetched["chain"]}).get_json()

    shell = client.post("/api/pocket", json={
        "structure_name": fetched["structure_name"], "chain": fetched["chain"],
        "residues": site["residues"]}).get_json()
    sized = client.post("/api/pocket", json={
        "structure_name": fetched["structure_name"], "chain": fetched["chain"],
        "residues": site["residues"],
        "size_from_reference": {"pdb_id": "1M17", "ligand_ccd": "AQ4"}}).get_json()

    def volume(box):
        return box[0] * box[1] * box[2]

    assert sized["center"] == shell["center"], "the centre must still come from the model"
    assert sized["sized_from"] == "1M17 AQ4"
    assert volume(sized["box"]) < volume(shell["box"]) / 3
    # Still large enough for the ligand to rotate in.
    assert min(sized["box"]) >= 18.0


def test_bundle_carries_no_hidden_files(client, app):
    """macOS writes .DS_Store into any directory Finder has opened, and one
    shipped inside the first real bundle."""
    response = client.post("/api/bundle", json={
        "protein": {"sequence": "MRPSGTAGAALLALL", "chain": "A"},
        "ligand": {"name": "x", "smiles": "CCO"},
        "pocket": {"residues": ["A:1"], "center": [0, 0, 0], "box": [18, 18, 18]}})
    job_id = response.get_json()["job_id"]
    from app import config
    with tarfile.open(config.RUNS_DIR / job_id / f"run_bundle_{job_id}.tar.gz") as tar:
        names = [n.split("/", 1)[-1] for n in tar.getnames()]
    hidden = [n for n in names if any(part.startswith(".") for part in n.split("/") if part)]
    assert not hidden, f"hidden files in the bundle: {hidden}"


def test_bundle_environments_do_not_share_a_solve_group():
    """A shared solve-group makes pixi resolve every environment against one
    locked set, so the default environment inherits the GNN extra's
    torch-scatter and fails to install with `No module named 'torch'` while
    building an environment that was meant not to contain torch at all."""
    import tomllib

    from app import config

    manifest = tomllib.loads((config.BUNDLE_TEMPLATE_DIR / "pixi.toml").read_text())
    environments = manifest.get("environments", {})
    assert environments, "the bundle should offer more than one environment"
    for name, spec in environments.items():
        assert "solve-group" not in spec, f"{name} shares a solve group"
    assert not (environments["default"].get("features") or [])


def test_the_run_command_is_one_line_and_self_contained(client, app):
    """Four steps (download, find the file, untar, cd) are four chances to end
    up in the wrong directory, and none of them is interesting."""
    response = client.post("/api/bundle", json={
        "protein": {"sequence": "MRPSGTAGAALLALL", "chain": "A"},
        "ligand": {"name": "x", "smiles": "CCO"},
        "pocket": {"residues": ["A:1"], "center": [0, 0, 0], "box": [18, 18, 18]},
        "md": {"production_ps": 500, "equilibration_ps": 100},
        "docking": {"mode": "dock"}})
    body = response.get_json()
    command = body["command"]
    assert command.count("&&") == 2, "should be fetch, cd, run and nothing else"
    assert command.startswith("curl -fL "), \
        "-f, or an HTTP error is piped into tar and reported as a gzip problem"
    assert "| tar xz" in command
    assert f"cd run_bundle_{body['job_id']}" in command
    assert command.endswith("./run.sh")
    # An absolute URL, or pasting it into a terminal fetches nothing.
    assert "http://localhost/runs/" in command or "https://" in command
    assert body["results_path"].endswith("results/results.tar.gz")
    assert 30 <= body["estimate_minutes"] <= 120, body["estimate_minutes"]


def test_a_private_bundle_command_carries_its_token(client):
    """A private bundle is guarded, so the command has to authenticate."""
    response = client.post("/api/bundle", json={
        "visibility": "private",
        "protein": {"sequence": "MRPSGTAGAALLALL", "chain": "A"},
        "ligand": {"name": "x", "smiles": "CCO"},
        "pocket": {"residues": ["A:1"], "center": [0, 0, 0], "box": [18, 18, 18]}})
    body = response.get_json()
    assert f"token={body['owner_token']}" in body["command"]


def test_a_public_bundle_command_carries_no_token(client):
    response = client.post("/api/bundle", json={
        "visibility": "public",
        "protein": {"sequence": "MRPSGTAGAALLALL", "chain": "A"},
        "ligand": {"name": "x", "smiles": "CCO"},
        "pocket": {"residues": ["A:1"], "center": [0, 0, 0], "box": [18, 18, 18]}})
    body = response.get_json()
    assert "token=" not in body["command"]


def test_the_bundle_carries_its_lock_and_an_executable_bootstrap(client, app):
    """One step means no prerequisites: run.sh installs pixi if the machine has
    not got it, and pixi.lock means no solve, so the environment is the one this
    was tested with rather than whatever resolves today."""
    from app import config

    response = client.post("/api/bundle", json={
        "protein": {"sequence": "MRPSGTAGAALLALL", "chain": "A"},
        "ligand": {"name": "x", "smiles": "CCO"},
        "pocket": {"residues": ["A:1"], "center": [0, 0, 0], "box": [18, 18, 18]}})
    job_id = response.get_json()["job_id"]
    with tarfile.open(config.RUNS_DIR / job_id / f"run_bundle_{job_id}.tar.gz") as tar:
        members = {m.name.split("/", 1)[-1]: m for m in tar.getmembers()}
        lock = tar.extractfile(members["pixi.lock"]).read().decode()
    # A tar that loses the executable bit turns one step back into two.
    assert members["run.sh"].mode & 0o111, "run.sh is not executable inside the archive"
    assert "platforms:" in lock and "osx-arm64" in lock and "linux-64" in lock, \
        "the lock must cover both platforms the bundle claims to support"
    assert len(lock) > 100_000, "that lock looks empty"
