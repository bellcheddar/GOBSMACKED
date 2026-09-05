"""The contract between the bundle and the server.

`REQUIRED` is the same list `app/services/ingest.py` validates against. If a
stage cannot produce one of these files it must fail loudly here rather than
letting the archive be written incomplete: an archive that unpacks and then
fails validation on the server wastes the user's upload as well as their run.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "1.0"

REQUIRED = [
    "manifest.json",
    "campaign.yaml",
    "model_apo.pdb",
    "poses/poses.sdf",
    "poses/scores.csv",
    "complex_pose1.pdb",
    "complex_min.pdb",
    "complex_md_final.pdb",
    "traj/topology.pdb",
    "traj/summary.json",
]

OPTIONAL = ["plddt.json", "traj/traj.dcd", "logs/run.log"]

STAGES = ["fold", "prep", "dock", "md", "summarise"]


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def engine_versions() -> dict[str, str]:
    """Whatever is installed, recorded so a result can be reproduced.

    Every import is guarded: a bundle that skipped folding has no torch, and
    that is not a reason to fail while writing the manifest.
    """
    versions: dict[str, str] = {
        "gobsmacked_run": "1.0.0",
        "python": platform.python_version(),
        "platform": f"{platform.system()} {platform.machine()}",
    }
    for name, module in (("rdkit", "rdkit"), ("openmm", "openmm"), ("mdtraj", "mdtraj"),
                         ("torch", "torch"), ("pdbfixer", "pdbfixer"),
                         ("openmmforcefields", "openmmforcefields"), ("esm", "esm")):
        try:
            mod = __import__(module)
            versions[name] = getattr(mod, "__version__", "installed")
        except Exception:
            continue
    try:
        out = subprocess.run(["pandadock", "--version"], capture_output=True, text=True, timeout=60)
        if out.returncode == 0:
            versions["pandadock"] = out.stdout.strip().splitlines()[-1]
    except Exception:
        pass
    return versions


def campaign_sha256(campaign_path: Path) -> str:
    return hashlib.sha256(campaign_path.read_bytes()).hexdigest()


def write_manifest(results: Path, job_id: str, campaign_path: Path,
                   timings: dict[str, float], warnings: list[str]) -> Path:
    manifest = {
        "schema": SCHEMA_VERSION,
        "job_id": job_id,
        "written": now(),
        "campaign_sha256": campaign_sha256(campaign_path),
        "engine_versions": engine_versions(),
        "timings": {k: round(v, 1) for k, v in timings.items()},
        "warnings": warnings,
        "command": " ".join(sys.argv),
    }
    path = results / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2))
    return path


def check_complete(results: Path) -> list[str]:
    """Required files that are missing. Empty means the archive is valid."""
    return [name for name in REQUIRED if not (results / name).exists()]
