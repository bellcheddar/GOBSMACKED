"""About: the pipeline schematic, the software table and the credits.

The software table is generated from `software.yaml` so that it and
THIRD_PARTY.md cannot drift: one file, two renderings.
"""

from __future__ import annotations

from pathlib import Path

import yaml
from flask import Blueprint, render_template

from .. import config

bp = Blueprint("about", __name__)

SOFTWARE_YAML = config.ROOT_DIR / "software.yaml"


def load_software() -> list[dict]:
    if not SOFTWARE_YAML.exists():
        return []
    data = yaml.safe_load(SOFTWARE_YAML.read_text()) or {}
    return data.get("software", [])


@bp.route("/about")
def about_page():
    software = load_software()
    stages: list[str] = []
    for item in software:
        if item.get("stage") and item["stage"] not in stages:
            stages.append(item["stage"])
    grouped = [(s, [i for i in software if i.get("stage") == s]) for s in stages]
    strip = [
        {"name": "Fetch", "text": "UniProt, RCSB, AFDB", "state": "done"},
        {"name": "Annotate", "text": "InterPro, KLIFS, GPCRdb", "state": "done"},
        {"name": "Fold", "text": "ESMFold", "state": "done"},
        {"name": "Dock", "text": "PandaDock", "state": "done"},
        {"name": "MD", "text": "OpenMM, OpenFF", "state": "done"},
        {"name": "Verify", "text": "biotite, tmtools, PLIP", "state": "done"},
        {"name": "Mode", "text": "KLIFS, GPCRdb", "state": "done"},
    ]
    return render_template("about.html", tab="about", stages=strip, grouped=grouped,
                           version=config.VERSION)
