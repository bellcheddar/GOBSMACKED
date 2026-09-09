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


def test_flex_output_is_looked_for_where_flex_writes_it(tmp_path):
    """`-o` is a directory for dock and hybrid and an output PREFIX for
    pandadock-flex, which appends "_results". A flex run told to write into
    work/docking writes into work/docking_results, so every flex run failed
    after finishing: 45 minutes of completed work in a directory one name away
    from the one being checked."""
    asked = tmp_path / "docking"
    asked.mkdir()
    real = tmp_path / "docking_results"
    real.mkdir()
    (real / "poses.sdf").write_text("", encoding="utf-8")

    assert dock.actual_output_dir(asked, "flex") == real
    # dock and hybrid are unaffected, even with a sibling sitting there.
    assert dock.actual_output_dir(asked, "dock") == asked
    assert dock.actual_output_dir(asked, "hybrid") == asked
    # And when flex did write where it was asked, that is what is used.
    (asked / "poses.sdf").write_text("", encoding="utf-8")
    assert dock.actual_output_dir(asked, "flex") == asked


# --- every mode names its top complex differently ----------------------------

def test_the_top_complex_is_found_whatever_the_mode_called_it(tmp_path):
    """dock and hybrid write complex1.pdb; flex writes
    complexes/flex_complex_1.pdb. Globbing "complex*.pdb" matched the first two
    and missed the third, so a flex run that had succeeded died with "PandaDock
    wrote no complex for the top pose" after forty-five minutes, with all ten
    complexes on disk under a name the glob could not see."""
    from gobsmacked_run.dock import find_top_complex

    plain = tmp_path / "plain"
    plain.mkdir()
    for i in (1, 2, 10):
        (plain / f"complex{i}.pdb").write_text("ATOM\n")
    assert find_top_complex(plain).name == "complex1.pdb"

    flex = tmp_path / "flex" / "complexes"
    flex.mkdir(parents=True)
    for i in range(1, 11):
        (flex / f"flex_complex_{i}.pdb").write_text("ATOM\n")
        (flex / f"flex_ligand_{i}.pdb").write_text("ATOM\n")
    assert find_top_complex(tmp_path / "flex").name == "flex_complex_1.pdb"


def test_a_ligand_only_file_is_never_mistaken_for_a_complex(tmp_path):
    """flex writes flex_ligand_N.pdb beside the complexes. Handing one to MD
    would solvate a ligand with no receptor rather than fail."""
    from gobsmacked_run.dock import find_top_complex
    d = tmp_path / "d" / "complexes"
    d.mkdir(parents=True)
    (d / "flex_ligand_1.pdb").write_text("ATOM\n")
    assert find_top_complex(tmp_path / "d") is None
    (d / "flex_complex_3.pdb").write_text("ATOM\n")
    assert find_top_complex(tmp_path / "d").name == "flex_complex_3.pdb"


def test_ten_does_not_sort_before_one(tmp_path):
    from gobsmacked_run.dock import find_top_complex
    d = tmp_path / "d"
    d.mkdir()
    for i in (10, 2, 1):
        (d / f"complex{i}.pdb").write_text("ATOM\n")
    assert find_top_complex(d).name == "complex1.pdb"


# --- Vinardo decides which pose goes forward ---------------------------------
# Measured across nine pose sets: the engine put a pose within 2 A of the crystal
# first once; Vinardo did it four times; and on the run the engine got right,
# Vinardo picked the same pose.

def test_poses_are_reordered_not_relabelled(tmp_path):
    """"Pose 1" has to mean the same thing in poses.sdf, scores.csv,
    complex_pose1.pdb and on the results page. Leaving the file in engine order
    while carrying a different complex forward would put the overlay and the
    scorecard on different molecules."""
    from gobsmacked_run.dock import reorder_poses
    sdf = tmp_path / "poses.sdf"
    sdf.write_text("".join(f"MOL{i}\nblock\n$$$$\n" for i in (1, 2, 3)), encoding="utf-8")
    reorder_poses(sdf, [3, 1, 2])
    records = [r.strip() for r in sdf.read_text(encoding="utf-8").split("$$$$") if r.strip()]
    assert [r.splitlines()[0] for r in records] == ["MOL3", "MOL1", "MOL2"]


def test_a_count_mismatch_leaves_the_file_untouched(tmp_path):
    """Silently reordering a file whose records do not match the scores would
    scramble which pose is which, which is worse than not reordering."""
    from gobsmacked_run.dock import reorder_poses
    sdf = tmp_path / "poses.sdf"
    original = "MOL1\n$$$$\nMOL2\n$$$$\n"
    sdf.write_text(original, encoding="utf-8")
    reorder_poses(sdf, [3, 2, 1])
    assert sdf.read_text(encoding="utf-8") == original


def test_the_engines_rank_is_kept_as_a_column(tmp_path):
    """Its opinion is being overruled, not deleted."""
    from gobsmacked_run.dock import reorder_scores
    rows = [{"pose_id": f"pose{i}", "score": -10.0 + i, "gnn_affinity": None, "rank": i}
            for i in (1, 2, 3)]
    csv_path = tmp_path / "scores.csv"
    moved = reorder_scores(csv_path, rows, [3, 1, 2])
    assert [r["pose_id"] for r in moved] == ["pose3", "pose1", "pose2"]
    assert [r["rank"] for r in moved] == [1, 2, 3]
    assert [r["engine_rank"] for r in moved] == [3, 1, 2]
    body = csv_path.read_text(encoding="utf-8")
    assert body.splitlines()[0].endswith("engine_rank")


def test_the_complex_follows_the_promoted_pose(tmp_path):
    """PandaDock names its complexes by ITS ranking, so promoting pose 8 means
    carrying complex8.pdb, not the first file on disk."""
    from gobsmacked_run.dock import find_top_complex
    d = tmp_path / "d"
    d.mkdir()
    for i in range(1, 11):
        (d / f"complex{i}.pdb").write_text("ATOM\n")
    assert find_top_complex(d, engine_rank=8).name == "complex8.pdb"
    assert find_top_complex(d).name == "complex1.pdb"
    # A rank with no file falls back rather than failing the run.
    assert find_top_complex(d, engine_rank=99).name == "complex1.pdb"
