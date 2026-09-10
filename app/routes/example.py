"""Example: two campaigns, six runs each, and what they disagree about.

Written for someone who has just arrived and wants to know whether any of this
works, not how it is built. Results and their meaning only: the About tab has
the pipeline and the thresholds, and each run links out to its own scorecard.

Both campaigns are held here as data rather than fetched from the database.
They are fixed published results, and a page describing them should not change
because somebody deleted a run or started a thirteenth.

Two targets rather than one because the first campaign's conclusions did not
survive the second. EGFR said co-folding decides the outcome and flexible
docking is the worst thing you can do; beta-2 adrenergic says the opposite on
both counts, with the same code and the same settings. What that buys the
reader is the honest version: these are target-specific observations, and the
only way to know which applies to yours is to check against a structure.
"""

from __future__ import annotations

from flask import Blueprint, render_template

bp = Blueprint("example", __name__)

# One campaign per target. Everything within a campaign is held constant except
# where the receptor came from and which docking mode searched it.
CAMPAIGNS = [
    {
        "key": "egfr",
        "family": "kinase",
        "target": "EGFR",
        "target_long": "EGFR kinase domain (P00533, residues 714-966)",
        "ligand": "erlotinib",
        "reference": "1M17",
        "reference_note": "2.6 &#8491;, holding erlotinib itself",
        "verdict": "Co-folding decided it. Nothing else did.",
        "runs": [
            {"name": "Co-folded receptor, standard docking", "job": "gs_20260909_irczwdpxzhs7",
             "receptor": "co-folded", "mode": "dock", "grade": "B", "score": 87.4,
             "rmsd": 1.48, "pocket_ca": 0.649, "minutes": 61,
             "mode_match": False, "predicted": "I 1/2", "crystal_label": "I"},
            {"name": "Co-folded receptor, hybrid docking", "job": "gs_20260909_zho3oqr4hlem",
             "receptor": "co-folded", "mode": "hybrid", "grade": "B", "score": 84.8,
             "rmsd": 1.62, "pocket_ca": 0.817, "minutes": 66,
             "mode_match": True, "predicted": "I", "crystal_label": "I"},
            {"name": "Co-folded receptor, flexible docking", "job": "gs_20260909_jq5tpzjs2xf6",
             "receptor": "co-folded", "mode": "flex", "grade": "D", "score": 49.6,
             "rmsd": 8.58, "pocket_ca": 1.308, "minutes": 99,
             "mode_match": True, "predicted": "I", "crystal_label": "I"},
            {"name": "AlphaFold receptor, standard docking", "job": "gs_20260909_hggp2krjwh3e",
             "receptor": "AlphaFold", "mode": "dock", "grade": "D", "score": 51.2,
             "rmsd": 4.9, "pocket_ca": 0.908, "minutes": 60,
             "mode_match": False, "predicted": "allosteric", "crystal_label": "I"},
            {"name": "AlphaFold receptor, hybrid docking", "job": "gs_20260909_mhxgoie3trgt",
             "receptor": "AlphaFold", "mode": "hybrid", "grade": "D", "score": 56.2,
             "rmsd": 7.44, "pocket_ca": 0.908, "minutes": 63,
             "mode_match": True, "predicted": "I", "crystal_label": "I"},
            {"name": "AlphaFold receptor, flexible docking", "job": "gs_20260909_lpqzx23moxjc",
             "receptor": "AlphaFold", "mode": "flex", "grade": "D", "score": 51.2,
             "rmsd": 8.38, "pocket_ca": 1.61, "minutes": 98,
             "mode_match": True, "predicted": "I", "crystal_label": "I"}
        ],
    },
    {
        "key": "b2ar",
        "family": "GPCR",
        "target": "&beta;2AR",
        "target_long": "Beta-2 adrenergic receptor (P07550, residues 50-326)",
        "ligand": "carazolol",
        "reference": "2RH1",
        "reference_note": "2.4 &#8491;, holding carazolol itself",
        "verdict": "The docking mode decided it. The receptor did not.",
        "runs": [
            {"name": "Co-folded receptor, standard docking", "job": "gs_20260910_2iy7ep3ulygd",
             "receptor": "co-folded", "mode": "dock", "grade": "D", "score": 60.4,
             "rmsd": 5.69, "pocket_ca": 0.307, "minutes": 165,
             "mode_match": True, "predicted": "orthosteric", "crystal_label": "orthosteric"},
            {"name": "Co-folded receptor, hybrid docking", "job": "gs_20260910_j3uwomcgpz2n",
             "receptor": "co-folded", "mode": "hybrid", "grade": "B", "score": 91.6,
             "rmsd": 1.13, "pocket_ca": 0.323, "minutes": 160,
             "mode_match": True, "predicted": "orthosteric", "crystal_label": "orthosteric"},
            {"name": "Co-folded receptor, flexible docking", "job": "gs_20260910_lutzsft3hfh4",
             "receptor": "co-folded", "mode": "flex", "grade": "B", "score": 84.1,
             "rmsd": 1.41, "pocket_ca": 0.767, "minutes": 173,
             "mode_match": True, "predicted": "orthosteric", "crystal_label": "orthosteric"},
            {"name": "AlphaFold receptor, standard docking", "job": "gs_20260910_453twmkb6yns",
             "receptor": "AlphaFold", "mode": "dock", "grade": "B", "score": 83.2,
             "rmsd": 2.84, "pocket_ca": 0.322, "minutes": 175,
             "mode_match": True, "predicted": "orthosteric", "crystal_label": "orthosteric"},
            {"name": "AlphaFold receptor, hybrid docking", "job": "gs_20260910_g43lmwd5cshc",
             "receptor": "AlphaFold", "mode": "hybrid", "grade": "B", "score": 87.4,
             "rmsd": 1.41, "pocket_ca": 0.322, "minutes": 154,
             "mode_match": True, "predicted": "orthosteric", "crystal_label": "orthosteric"},
            {"name": "AlphaFold receptor, flexible docking", "job": "gs_20260910_ffbrd5lvtpsb",
             "receptor": "AlphaFold", "mode": "flex", "grade": "A", "score": 94.2,
             "rmsd": 0.9, "pocket_ca": 1.182, "minutes": 91,
             "mode_match": True, "predicted": "orthosteric", "crystal_label": "orthosteric"}
        ],
    },
]

# The 2x3 grid each campaign's table is also drawn as, so the inversion between
# the two targets is visible without reading twelve rows.
MODES = ["dock", "hybrid", "flex"]
RECEPTORS = ["co-folded", "AlphaFold"]


def grid(campaign: dict) -> list[dict]:
    """Runs arranged receptor x mode, for the comparison grid."""
    by_key = {(r["receptor"], r["mode"]): r for r in campaign["runs"]}
    return [{"receptor": rec, "cells": [by_key.get((rec, m)) for m in MODES]}
            for rec in RECEPTORS]


@bp.route("/example")
def example_page():
    campaigns = [{**c, "grid": grid(c)} for c in CAMPAIGNS]
    return render_template("example.html", tab="example",
                           campaigns=campaigns, modes=MODES)
