"""Example: one campaign, six runs, and what came out of it.

Written for someone who has just arrived and wants to know whether any of this
works, not how it is built. Results and their meaning only: the About tab has
the pipeline and the thresholds, and each run links out to its own scorecard.

The six runs are held here as data rather than fetched from the database. They
are a fixed published result, and a page describing them should not change
because somebody deleted a run or started a seventh.
"""

from __future__ import annotations

from flask import Blueprint, render_template

bp = Blueprint("example", __name__)

# EGFR kinase domain plus erlotinib, judged against crystal 1M17. One campaign,
# six runs on an M1 Max, varying only where the receptor came from and which
# docking mode searched it.
RUNS = [
    {"name": "Co-folded receptor, hybrid docking", "job": "gs_20260909_zho3oqr4hlem",
     "receptor": "co-folded", "mode": "hybrid", "grade": "B", "score": 84.8,
     "rmsd": 1.66, "minutes": 66},
    {"name": "Co-folded receptor, standard docking", "job": "gs_20260909_irczwdpxzhs7",
     "receptor": "co-folded", "mode": "dock", "grade": "B", "score": 87.4,
     "rmsd": 1.59, "minutes": 61},
    {"name": "Co-folded receptor, flexible docking", "job": "gs_20260909_jq5tpzjs2xf6",
     "receptor": "co-folded", "mode": "flex", "grade": "D", "score": 49.6,
     "rmsd": 8.58, "minutes": 99},
    {"name": "AlphaFold receptor, hybrid docking", "job": "gs_20260909_mhxgoie3trgt",
     "receptor": "AlphaFold", "mode": "hybrid", "grade": "D", "score": 56.2,
     "rmsd": 7.44, "minutes": 63},
    {"name": "AlphaFold receptor, standard docking", "job": "gs_20260909_hggp2krjwh3e",
     "receptor": "AlphaFold", "mode": "dock", "grade": "D", "score": 51.2,
     "rmsd": 4.90, "minutes": 60},
    {"name": "AlphaFold receptor, flexible docking", "job": "gs_20260909_lpqzx23moxjc",
     "receptor": "AlphaFold", "mode": "flex", "grade": "D", "score": 51.2,
     "rmsd": 8.44, "minutes": 98},
]


@bp.route("/example")
def example_page():
    return render_template("example.html", tab="example", runs=RUNS)
