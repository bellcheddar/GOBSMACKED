"""The trajectory clip.

Cheap to test on the EGFR fixture, which is a real 20-frame trajectory, so this
covers the actual render rather than a mock: the file is opened afterwards and
asked what it contains.
"""

from __future__ import annotations

import shutil
import subprocess
import tarfile

import pytest

from app.services import movie

pytestmark = pytest.mark.skipif(shutil.which("ffmpeg") is None,
                                reason="ffmpeg is not installed")


@pytest.fixture()
def results(tmp_path, egfr_archive):
    with tarfile.open(egfr_archive) as tar:
        tar.extractall(tmp_path)
    root = tmp_path / "results" if (tmp_path / "results").is_dir() else tmp_path
    return root


def probe(path) -> dict:
    out = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "v:0",
                          "-show_entries", "stream=width,height,nb_frames,pix_fmt",
                          "-of", "default=noprint_wrappers=1", str(path)],
                         capture_output=True, text=True)
    return dict(line.split("=", 1) for line in out.stdout.strip().splitlines())


def test_renders_a_playable_clip(results, tmp_path):
    out = movie.render(results / "traj" / "topology.pdb", results / "traj" / "traj.dcd",
                       tmp_path / "motion.mp4", tmp_path / "motion_poster.webp")
    assert "error" not in out, out
    assert (tmp_path / "motion.mp4").stat().st_size > 5_000
    assert (tmp_path / "motion_poster.webp").stat().st_size > 1_000

    info = probe(tmp_path / "motion.mp4")
    assert (int(info["width"]), int(info["height"])) == (movie.WIDTH, movie.HEIGHT)
    # yuv420p or Safari refuses to decode it at all.
    assert info["pix_fmt"] == "yuv420p"
    # Forwards then backwards: 2n - 2 frames from n frames of trajectory.
    assert int(info["nb_frames"]) == out["frames"] * 2 - 2


def test_finds_the_ligand_and_its_site(results, tmp_path):
    out = movie.render(results / "traj" / "topology.pdb", results / "traj" / "traj.dcd",
                       tmp_path / "motion.mp4", tmp_path / "motion_poster.webp")
    # A pocket of zero residues means the ligand was not found, which is the
    # failure that would silently produce a clip of a protein and nothing else.
    assert out["pocket_residues"] > 4


def test_declines_without_a_trajectory(tmp_path):
    out = movie.render(tmp_path / "nope.pdb", tmp_path / "nope.dcd",
                       tmp_path / "motion.mp4", tmp_path / "poster.webp")
    assert "error" in out
    assert not (tmp_path / "motion.mp4").exists()


def test_declines_past_the_budget(results, tmp_path, monkeypatch):
    monkeypatch.setattr(movie, "BUDGET_ATOM_FRAMES", 10)
    out = movie.render(results / "traj" / "topology.pdb", results / "traj" / "traj.dcd",
                       tmp_path / "motion.mp4", tmp_path / "poster.webp")
    assert "budget" in out["error"]
