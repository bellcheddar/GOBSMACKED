"""The pose overlay: every pose docking made, and which failure it was.

The interesting assertions here are about the verdict, because an earlier
version read RMSD alone and told a run whose poses sat 2.4 A from the crystal
ligand's centre, fully inside the docking box, to "check the box centre". The
box was right. Three numbers are needed to tell the failures apart and the tests
below are one per failure.
"""

from __future__ import annotations

from app.services import poses as poses_svc


def rows(**overrides):
    base = {"rank": 1, "pose_id": "pose1", "score": -7.0, "gnn_energy": -7.0,
            "vina_energy": -15.0, "gnn_pec50": 5.1, "confidence": 0.97,
            "rmsd": 8.2, "centroid_distance": 2.4, "shape_rmsd": 2.5,
            "closest_contact": 3.3}
    return [{**base, **overrides}]


def test_a_correct_pose_is_called_correct():
    verdict = poses_svc._verdict(rows(rmsd=1.2), {"rank": 1, "rmsd": 1.2})
    assert "the search found it and the scoring function agreed" in verdict


def test_a_mis_ranked_crystal_pose_is_a_scoring_failure():
    all_rows = rows() + rows(rank=7, rmsd=1.4, centroid_distance=0.9, shape_rmsd=1.1)
    verdict = poses_svc._verdict(all_rows, all_rows[1])
    assert "scoring failure" in verdict
    assert "ranked it 7" in verdict


def test_the_right_pocket_and_the_wrong_orientation_is_not_a_box_problem():
    """The case that produced the wrong advice: centres 2.4 A apart, every atom
    inside the box, and the old verdict said the search never went there."""
    all_rows = rows()
    verdict = poses_svc._verdict(all_rows, all_rows[0])
    assert "in the right pocket" in verdict
    assert "turned the wrong way round" in verdict
    assert "Check the box centre before anything else" not in verdict


def test_a_genuinely_wrong_pocket_still_says_check_the_box():
    all_rows = rows(rmsd=14.0, centroid_distance=13.5, shape_rmsd=2.4)
    verdict = poses_svc._verdict(all_rows, all_rows[0])
    assert "not in the crystal's pocket at all" in verdict
    assert "Check the box centre" in verdict


def test_no_crystal_means_no_claim():
    all_rows = rows(rmsd=None, centroid_distance=None, shape_rmsd=None)
    assert "No crystal ligand" in poses_svc._verdict(all_rows, None)


def test_the_shape_rmsd_is_free_to_superpose_and_the_main_one_is_not():
    """GetBestRMS superposes before measuring, so it answers "is this the same
    shape" and never "is this in the right place". Using it for the main RMSD
    would report every pose as nearly correct."""
    import inspect

    source = inspect.getsource(poses_svc._rmsd_to)
    assert "GetBestRMS" in source and "best_fit" in source
    # The default path must be the in-place one.
    assert "if best_fit:" in source
