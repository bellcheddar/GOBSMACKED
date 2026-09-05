"""Stage 3: dock.

PandaDock is run as a subprocess, in the mode the campaign asked for:

    hybrid  search with the empirical function, rank with the SE(3) GNN
    dock    empirical search and scoring only, no GNN model needed
    flex    induced fit, refining receptor side chains around each pose

The GNN checkpoint is about 82 MB and is downloaded on first use, so `hybrid`
falls back to `dock` when there is no network and no cached model rather than
failing the whole run four stages in.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any, Optional

MODEL_NAME = "pandadock_gnn.pt"


def run(campaign: dict, work: Path, results: Path, log) -> dict[str, Any]:
    docking = campaign.get("docking") or {}
    pocket = campaign.get("pocket") or {}
    mode = docking.get("mode", "hybrid")
    centre = pocket.get("center")
    box = pocket.get("box") or [22, 22, 22]
    if not centre:
        raise RuntimeError("The campaign has no pocket centre, so there is nothing to dock into.")

    receptor = work / "receptor.pdb"
    ligand = work / "ligand.sdf"
    out_dir = work / "docking"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    warnings: list[str] = []

    model = None
    if mode == "hybrid":
        model = ensure_gnn_model(work, log)
        if model is None:
            mode = "dock"
            warnings.append("No PandaDock GNN checkpoint and no network to fetch one: the poses "
                            "were searched and ranked by the empirical function alone.")

    cmd = build_command(mode, receptor, ligand, centre, box, docking, out_dir, model)
    log("dock: " + " ".join(str(c) for c in cmd))
    proc = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    (work / "pandadock.log").write_text((proc.stdout or "") + "\n" + (proc.stderr or ""))
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-3:]
        raise RuntimeError("PandaDock failed: " + " / ".join(tail))

    poses_dir = results / "poses"
    poses_dir.mkdir(parents=True, exist_ok=True)
    poses = out_dir / "poses.sdf"
    if not poses.exists():
        raise RuntimeError(f"PandaDock wrote no poses.sdf in {out_dir}.")
    shutil.copy(poses, poses_dir / "poses.sdf")

    scores = write_scores(out_dir, poses_dir / "scores.csv", log)
    top_complex = find_top_complex(out_dir)
    if top_complex is None:
        raise RuntimeError("PandaDock wrote no complex for the top pose.")
    shutil.copy(top_complex, results / "complex_pose1.pdb")

    log(f"dock: {len(scores)} poses, best score {scores[0]['score'] if scores else '?'}")
    return {"warnings": warnings, "mode": mode, "poses": len(scores),
            "best_score": scores[0]["score"] if scores else None}


def build_command(mode: str, receptor: Path, ligand: Path, centre, box, docking: dict,
                  out_dir: Path, model: Optional[Path]) -> list:
    common = ["-r", receptor, "-l", ligand, "-o", out_dir,
              "--center", centre[0], centre[1], centre[2]]
    if mode == "flex":
        # pandadock-flex takes a radius, not a box: half the largest side is the
        # sphere that contains the campaign's box.
        radius = round(max(box) / 2.0, 1)
        return ["pandadock-flex", *common, "--radius", radius,
                "--initial-poses-to-retain", min(int(docking.get("num_poses", 10)), 5)]
    cmd = ["pandadock", "hybrid" if mode == "hybrid" else "dock", *common,
           "--box", box[0], box[1], box[2],
           "-n", int(docking.get("num_poses", 10)),
           # A fixed seed: two runs of the same campaign should return the same
           # poses, or the scorecard is measuring the sampler's variance.
           "--seed", 20260905]
    exhaustiveness = docking.get("exhaustiveness")
    if exhaustiveness:
        cmd += ["-e", int(exhaustiveness)]
    if mode == "hybrid" and model is not None:
        cmd += ["-m", model]
    return cmd


def ensure_gnn_model(work: Path, log) -> Optional[Path]:
    """The GNN checkpoint, downloading it once if it is not already here."""
    for candidate in (work / "models" / MODEL_NAME,
                      Path.home() / ".pandadock" / MODEL_NAME,
                      Path("models") / MODEL_NAME):
        if candidate.exists():
            return candidate
    log("dock: fetching the PandaDock GNN checkpoint (about 82 MB, once)")
    try:
        proc = subprocess.run(["pandadock", "gnn", "download-model"],
                              capture_output=True, text=True, timeout=1800)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    if proc.returncode != 0:
        return None
    for candidate in (Path("models") / MODEL_NAME, Path.home() / ".pandadock" / MODEL_NAME):
        if candidate.exists():
            return candidate
    return None


def write_scores(out_dir: Path, dest: Path, log) -> list[dict]:
    """pose_id, score, GNN affinity and rank, from whatever PandaDock wrote.

    PandaDock's own JSON is preferred; the SD tags in poses.sdf are the fallback,
    because the JSON's exact filename has changed between versions and the tags
    have not.
    """
    rows: list[dict] = []
    for path in sorted(out_dir.glob("*_poses.json")) + sorted(out_dir.glob("*poses*.json")):
        try:
            data = json.loads(path.read_text())
        except (json.JSONDecodeError, OSError):
            continue
        poses = data.get("poses") if isinstance(data, dict) else data
        if not isinstance(poses, list):
            continue
        for index, pose in enumerate(poses, start=1):
            if not isinstance(pose, dict):
                continue
            rows.append({
                "pose_id": pose.get("pose_id") or pose.get("id") or f"pose{index}",
                "score": _number(pose.get("score", pose.get("energy"))),
                "gnn_affinity": _number(pose.get("gnn_score", pose.get("predicted_affinity",
                                                                       pose.get("pec50")))),
                "rank": pose.get("rank", index),
            })
        if rows:
            break

    if not rows:
        rows = scores_from_sdf(out_dir / "poses.sdf")

    with open(dest, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["pose_id", "score", "gnn_affinity", "rank"])
        writer.writeheader()
        writer.writerows(rows)
    return rows


def scores_from_sdf(path: Path) -> list[dict]:
    if not path.exists():
        return []
    from rdkit import Chem

    rows = []
    for index, mol in enumerate(Chem.SDMolSupplier(str(path), sanitize=False), start=1):
        if mol is None:
            continue
        props = mol.GetPropsAsDict()
        rows.append({
            "pose_id": props.get("pose_id", mol.GetProp("_Name") if mol.HasProp("_Name") else f"pose{index}"),
            "score": _number(props.get("score", props.get("Score", props.get("energy")))),
            "gnn_affinity": _number(props.get("gnn_score", props.get("predicted_affinity"))),
            "rank": props.get("rank", index),
        })
    return rows


def find_top_complex(out_dir: Path) -> Optional[Path]:
    candidates = sorted(out_dir.glob("complex*.pdb"))
    if not candidates:
        candidates = sorted(out_dir.rglob("complex*.pdb"))
    if not candidates:
        return None
    # complex1.pdb is rank 1; sorting lexically would put complex10 before it.
    def rank(path: Path) -> int:
        digits = "".join(c for c in path.stem if c.isdigit())
        return int(digits) if digits else 0
    return sorted(candidates, key=rank)[0]


def _number(value) -> Optional[float]:
    try:
        return round(float(value), 4)
    except (TypeError, ValueError):
        return None
