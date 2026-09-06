"""The bundle's dock stage: the command it builds and the checkpoint it needs.

Three bugs stacked here, each hidden by the one in front of it, and every one
of them ended the same way: a run that quietly used the weaker scorer and said
so in a warning that blamed the network.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bundle_template"))

from gobsmacked_run import dock  # noqa: E402

CENTRE = [1.0, 2.0, 3.0]
BOX = [22.0, 22.0, 22.0]


def build(mode: str, model=None, **docking):
    return [str(part) for part in dock.build_command(
        mode, Path("r.pdb"), Path("l.sdf"), CENTRE, BOX,
        {"num_poses": 3, "exhaustiveness": 16, **docking}, Path("out"), model)]


def test_hybrid_is_given_no_flags_it_rejects():
    """`pandadock hybrid` accepts neither --seed nor -e and exits on either,
    before doing any work. Both were being passed, so hybrid could not have run
    even once its checkpoint was found."""
    cmd = build("hybrid", model=Path("model.pt"))
    assert "--seed" not in cmd
    assert "-e" not in cmd and "--exhaustiveness" not in cmd
    assert cmd[:2] == ["pandadock", "hybrid"]
    assert "-m" in cmd and "model.pt" in cmd


def test_dock_keeps_its_seed_and_exhaustiveness():
    """`dock` takes both, and is the mode that stays reproducible."""
    cmd = build("dock")
    assert "--seed" in cmd
    assert cmd[cmd.index("-e") + 1] == "16"
    assert "-m" not in cmd


def test_hybrid_without_a_checkpoint_asks_for_no_model():
    cmd = build("hybrid", model=None)
    assert "-m" not in cmd


def test_the_checkpoint_is_named_as_the_release_publishes_it():
    """The old value was pandadock_gnn.pt, a file that exists nowhere, so a
    correctly downloaded checkpoint would still not have been found."""
    assert dock.MODEL_NAME == "pandadock_gnn_v4.pt"
    assert dock.FALLBACK_URL.endswith("/" + dock.MODEL_NAME)


def test_a_truncated_download_is_never_left_in_the_cache(tmp_path, monkeypatch):
    """An interrupted fetch that landed in the cache would be found by every
    later run and loaded as a model."""
    monkeypatch.setattr(dock, "MODEL_CACHE", tmp_path / "models")
    monkeypatch.setattr(dock, "newest_model_url", lambda log: "https://example.invalid/m.pt")

    def half_a_file(url, dest, log):
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"not a model")

    monkeypatch.setattr(dock, "download", half_a_file)
    assert dock.ensure_gnn_model(tmp_path / "work", lambda message: None) is None
    assert not list((tmp_path / "models").glob("*.pt"))


def test_a_download_that_raises_is_not_fatal(tmp_path, monkeypatch):
    """Falling back to the empirical scorer is a working run; raising here would
    end a campaign three stages in."""
    monkeypatch.setattr(dock, "MODEL_CACHE", tmp_path / "models")
    monkeypatch.setattr(dock, "newest_model_url", lambda log: None)

    def explode(url, dest, log):
        raise OSError("network went away")

    monkeypatch.setattr(dock, "download", explode)
    said = []
    assert dock.ensure_gnn_model(tmp_path / "work", said.append) is None
    assert any("could not be fetched" in line for line in said)


def test_a_cached_checkpoint_is_used_without_a_request(tmp_path, monkeypatch):
    monkeypatch.setattr(dock, "MODEL_CACHE", tmp_path / "models")
    (tmp_path / "models").mkdir(parents=True)
    cached = tmp_path / "models" / dock.MODEL_NAME
    cached.write_bytes(b"x" * (dock.MIN_MODEL_BYTES + 1))

    def never(*args, **kwargs):
        raise AssertionError("the cache should have answered")

    monkeypatch.setattr(dock, "newest_model_url", never)
    monkeypatch.setattr(dock, "download", never)
    assert dock.ensure_gnn_model(tmp_path / "work", lambda message: None) == cached


def test_hybrid_scores_are_read_from_the_file_hybrid_writes(tmp_path):
    """`pandadock hybrid` writes hybrid_results.csv and no JSON at all, with
    columns of its own. Nothing that read `dock`'s output found anything here,
    so a hybrid run came back with three poses and an empty score column, shown
    on the card as "best None kcal/mol"."""
    (tmp_path / "hybrid_results.csv").write_text(
        "rank,gnn_pec50,gnn_energy,vina_energy,activity_prob\n"
        "1,4.993,-6.819,-15.354,0.962\n"
        "2,4.784,-6.533,-15.584,0.932\n", encoding="utf-8")
    rows = dock.write_scores(tmp_path, tmp_path / "scores.csv", lambda m: None)
    assert [r["rank"] for r in rows] == [1, 2]
    # The score column is the GNN energy: in hybrid mode the GNN did the
    # ranking, so the score should be the number the ranking was made on.
    assert rows[0]["score"] == -6.819
    assert rows[0]["gnn_affinity"] == 4.993
    assert isinstance(rows[0]["rank"], int)


def test_pandadocks_own_sdf_tag_names_are_read(tmp_path):
    """score_kcal_per_mol and energy_gnn_pec50 are what it actually writes."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("_Name", "pose1")
    mol.SetProp("rank", "1")
    mol.SetProp("score_kcal_per_mol", "-6.819")
    mol.SetProp("energy_gnn_pec50", "4.993")
    path = tmp_path / "poses.sdf"
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()

    rows = dock.scores_from_sdf(path)
    assert rows[0]["score"] == -6.819
    assert rows[0]["gnn_affinity"] == 4.993
