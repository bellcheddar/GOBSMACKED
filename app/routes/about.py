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
        {"name": "Fetch", "state": "ready", "text": "UniProt, RCSB, AlphaFold DB, ESM Atlas"},
        {"name": "Annotate", "state": "ready", "text": "InterPro, KLIFS, GPCRdb"},
        {"name": "Fold", "state": "ready", "text": "ESMFold, in the bundle"},
        {"name": "Dock", "state": "ready", "text": "PandaDock, in the bundle"},
        {"name": "MD", "state": "ready", "text": "OpenMM and OpenFF, in the bundle"},
        {"name": "Verify", "state": "ready", "text": "biotite, tmtools, PLIP"},
        {"name": "Mode", "state": "ready", "text": "KLIFS and GPCRdb labels"},
    ]
    return render_template("about.html", tab="about", stages=strip, grouped=grouped,
                           version=config.VERSION)
