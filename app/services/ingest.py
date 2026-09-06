"""Validating and unpacking a results archive.

The archive comes from a machine this server has never seen, so it is treated as
untrusted input: extraction uses tarfile's `data` filter (no absolute paths, no
symlinks, no device files), the required files are checked by name before
anything is parsed, and a missing file is reported by name rather than as a
stack trace three stages later.
"""

from __future__ import annotations

import hashlib
import json
import tarfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import yaml

SCHEMA_VERSION = "1.0"

# Present in every archive the runner writes, regardless of campaign.
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

# Written when the campaign asked for them; their absence is a fact about the
# run, not a fault in the archive.
OPTIONAL = ["plddt.json", "traj/traj.dcd", "logs/run.log"]


class IngestError(ValueError):
    """An archive cannot be accepted. The message is shown to the user verbatim."""


@dataclass
class Results:
    root: Path                       # the extracted `results/` directory
    manifest: dict = field(default_factory=dict)
    campaign: dict = field(default_factory=dict)
    summary: dict = field(default_factory=dict)
    # Absent from every archive built before stage 5 existed, and from every run
    # that declined it. An empty dict, never a required file: a missing affinity
    # is a missing panel, not an invalid archive.
    affinity: dict = field(default_factory=dict)
    present: list[str] = field(default_factory=list)
    missing_optional: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def path(self, rel: str) -> Optional[Path]:
        p = self.root / rel
        return p if p.exists() else None

    @property
    def job_id(self) -> str:
        return self.manifest.get("job_id") or self.campaign.get("job_id") or ""

    @property
    def owner_token(self) -> str:
        # campaign.yaml is echoed back inside the archive, so uploading results
        # to a private run needs no extra typing.
        return self.campaign.get("owner_token") or ""


def _members(tar: tarfile.TarFile) -> dict[str, tarfile.TarInfo]:
    """Archive members keyed by their path below the single top-level directory."""
    out: dict[str, tarfile.TarInfo] = {}
    for m in tar.getmembers():
        parts = Path(m.name).parts
        if not parts:
            continue
        rel = "/".join(parts[1:]) if parts[0] in ("results", ".") else m.name
        out[rel] = m
    return out


def validate(archive: Path) -> list[str]:
    """Required files missing from the archive. Empty list means it is valid."""
    with tarfile.open(archive, "r:*") as tar:
        names = set(_members(tar))
    return [r for r in REQUIRED if r not in names]


def extract(archive: Path, dest: Path) -> Results:
    """Validate, extract and parse. Raises IngestError with a named file."""
    archive = Path(archive)
    dest = Path(dest)
    try:
        missing = validate(archive)
    except tarfile.TarError as exc:
        raise IngestError(f"That file is not a readable tar.gz archive ({exc}).") from exc
    if missing:
        raise IngestError(
            "The archive is missing " + ", ".join(missing[:4])
            + (f" and {len(missing) - 4} more" if len(missing) > 4 else "")
            + ". Re-run the bundle's summarise stage: it writes results.tar.gz complete or not at all."
        )

    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:*") as tar:
        members = _members(tar)
        for rel, member in members.items():
            if not member.isfile():
                continue
            target = (dest / rel).resolve()
            if not str(target).startswith(str(dest.resolve())):
                raise IngestError(f"The archive tries to write outside its directory ({rel}).")
            target.parent.mkdir(parents=True, exist_ok=True)
            src = tar.extractfile(member)
            if src is None:
                continue
            with open(target, "wb") as fh:
                fh.write(src.read())

    res = Results(root=dest, present=sorted(p for p in members))
    res.manifest = _read_json(dest / "manifest.json", "manifest.json")
    res.summary = _read_json(dest / "traj" / "summary.json", "traj/summary.json")
    res.affinity = _read_json_optional(dest / "affinity" / "affinity.json")
    try:
        res.campaign = yaml.safe_load((dest / "campaign.yaml").read_text()) or {}
    except yaml.YAMLError as exc:
        raise IngestError(f"campaign.yaml in the archive is not valid YAML ({exc}).") from exc

    schema = str(res.manifest.get("schema", ""))
    if schema and schema.split(".")[0] != SCHEMA_VERSION.split(".")[0]:
        raise IngestError(
            f"The archive declares results schema {schema}; this server reads {SCHEMA_VERSION}. "
            "Download a fresh bundle and re-run it."
        )

    manifest_job = res.manifest.get("job_id")
    campaign_job = res.campaign.get("job_id")
    if manifest_job and campaign_job and manifest_job != campaign_job:
        # The manifest wins (it is what `job_id` returns), but a disagreement
        # means the archive was repacked after the run, and the settings shown
        # on the results page may not be the settings that produced it.
        res.warnings.append(
            f"The manifest names job {manifest_job} and the campaign inside it names "
            f"{campaign_job}. This archive has been repacked since it was written."
        )

    res.missing_optional = [o for o in OPTIONAL if not (dest / o).exists()]
    if "traj/traj.dcd" in res.missing_optional:
        res.warnings.append("No trajectory file: the dynamics panels will show the summary only.")
    res.warnings.extend(res.manifest.get("warnings") or [])
    return res


def campaign_matches(res: Results, campaign_yaml: str) -> bool:
    """Does the echoed campaign match the one this server issued?

    A mismatch is not fatal: someone may legitimately edit the MD length before
    running. It is surfaced so the scorecard is not silently attributed to
    settings that were never used.
    """
    stored = hashlib.sha256((campaign_yaml or "").encode("utf-8")).hexdigest()
    claimed = res.manifest.get("campaign_sha256", "")
    return bool(claimed) and claimed == stored


def _read_json_optional(path: Path) -> dict:
    """A block that may not be there, and whose absence is not an error.

    `_read_json` raises on a missing file, which is right for the manifest and
    wrong for anything a run may legitimately not have produced. Using it for
    the affinity block would reject every archive built before stage 5 existed.
    """
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return {}
    return obj if isinstance(obj, dict) else {}


def _read_json(path: Path, label: str) -> dict:
    try:
        obj = json.loads(path.read_text())
    except FileNotFoundError:
        raise IngestError(f"The archive is missing {label}.") from None
    except json.JSONDecodeError as exc:
        raise IngestError(f"{label} in the archive is not valid JSON ({exc}).") from exc
    # json.loads("3") is a valid int, and .get() on it raises AttributeError
    # rather than saying anything useful.
    if not isinstance(obj, dict):
        raise IngestError(f"{label} should be a JSON object, not {type(obj).__name__}.")
    return obj
