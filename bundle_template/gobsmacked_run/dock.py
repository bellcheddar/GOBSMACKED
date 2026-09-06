"""Stage 3: dock.

PandaDock is run as a subprocess, in the mode the campaign asked for:

    hybrid  search with the empirical function, rank with the SE(3) GNN
    dock    empirical search and scoring only, no GNN model needed
    flex    induced fit, refining receptor side chains around each pose

The GNN checkpoint is about 82 MB and is fetched on first use into a cache
shared by every bundle on the machine, so `hybrid` falls back to `dock` only
when it genuinely cannot be had, rather than failing the whole run four stages
in.

It is fetched from the GitHub release directly rather than through `pandadock
gnn download-model`, which asks for release v4.0.0. That release carries no
assets at all; the checkpoint is published under v4.1.1. The command therefore
prints "Model not found at release URL" and exits 0, which read here as a
network failure and put every hybrid run onto the empirical scorer with a
warning blaming the user's connection.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Optional

# The asset's real name in the release. The old value here was
# "pandadock_gnn.pt", which is not a file that exists anywhere, so even a
# correctly downloaded checkpoint would not have been found in the cache.
MODEL_NAME = "pandadock_gnn_v4.pt"
RELEASES_API = "https://api.github.com/repos/pritampanda15/PandaDock/releases"
# Where the newest release carrying the asset was, when this was written. Used
# only when the API cannot be reached or answers with nothing usable.
FALLBACK_URL = ("https://github.com/pritampanda15/PandaDock/releases/download/"
                "v4.1.1/pandadock_gnn_v4.pt")
# Shared across bundles: the same 82 MB should not be fetched again for every
# campaign on one machine, and a bundle directory is a temporary thing.
MODEL_CACHE = Path.home() / ".gobsmacked" / "models"
MIN_MODEL_BYTES = 10_000_000


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
            warnings.append("The PandaDock GNN checkpoint could not be fetched, so the poses "
                            "were searched and ranked by the empirical function alone. The run "
                            "is sound; the ranking is the weaker of the two.")

    cmd = build_command(mode, receptor, ligand, centre, box, docking, out_dir, model)
    # num_poses, not n_poses: the campaign's key. Reading the wrong one printed
    # the default of 10 on a run whose command said -n 3, which is the kind of
    # log line that gets believed over the command beside it.
    poses = int(docking.get("num_poses", 10) or 10)
    log(f"dock: {mode} mode, {poses} poses"
        + (", GNN rescoring" if mode == "hybrid"
           else f", exhaustiveness {docking.get('exhaustiveness', 16)}"))
    log("dock: " + " ".join(str(c) for c in cmd))
    captured, returncode = run_with_progress(cmd, log, estimate_seconds(docking))
    (work / "pandadock.log").write_text(captured, encoding="utf-8")
    if returncode != 0:
        tail = captured.strip().splitlines()[-3:]
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

    best = scores[0]["score"] if scores else None
    return {"warnings": warnings, "mode": mode, "poses": len(scores),
            "best_score": best,
            "headline": f"{len(scores)} poses, best {best} kcal/mol" if scores else "no poses"}


def estimate_seconds(docking: dict) -> float:
    """How long PandaDock will take, from the one run that was timed.

    Measured on the EGFR campaign: 937 s for a 29 heavy-atom ligand at
    exhaustiveness 16, of which grid construction was 114 s and the search the
    rest. The search scales with exhaustiveness; the grid does not, it scales
    with the box, which is why it is a constant here rather than a coefficient.
    A prior, not a promise, and the bar says so when it runs past it.
    """
    exhaustiveness = float(docking.get("exhaustiveness", 16) or 16)
    return 114.0 + 823.0 * exhaustiveness / 16.0


def run_with_progress(cmd, log, estimate_s: float):
    """Run PandaDock, showing what it last said and how long it has been going.

    PandaDock prints "Starting docking..." and then nothing at all for a quarter
    of an hour, which is exactly the shape of output that gets a run killed by
    someone who thinks it has hung. Its stdout is drained on a thread so the
    pipe cannot fill and deadlock it, and the newest line it wrote is shown
    beside a clock.
    """
    import threading

    from .console import bar_for

    proc = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, bufsize=1)
    lines: list[str] = []

    def drain():
        for line in proc.stdout:                       # type: ignore[union-attr]
            lines.append(line.rstrip())

    reader = threading.Thread(target=drain, daemon=True)
    reader.start()
    with bar_for(log, "docking", estimate_s=estimate_s) as bar:
        while proc.poll() is None:
            latest = next((line for line in reversed(lines) if line.strip()), "")
            bar.update(note=latest[:48])
            time.sleep(0.2)
    reader.join(timeout=5)
    return "\n".join(lines) + "\n", proc.returncode


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
           "-n", int(docking.get("num_poses", 10))]
    if mode == "hybrid":
        # `pandadock hybrid` accepts neither --seed nor -e: it exits with "No
        # such option" before doing any work. Both were being passed to it, so
        # hybrid could never have run even once the checkpoint was found, and
        # the checkpoint bug was hiding the flag bug behind it.
        #
        # The cost is real and worth stating: hybrid runs are therefore not
        # seeded, so two runs of one campaign can return different poses and any
        # difference between them is the sampler's variance rather than a
        # finding. `dock` remains reproducible.
        if model is not None:
            cmd += ["-m", model]
        return cmd
    # A fixed seed: two runs of the same campaign should return the same poses,
    # or the scorecard is measuring the sampler's variance.
    cmd += ["--seed", 20260905]
    exhaustiveness = docking.get("exhaustiveness")
    if exhaustiveness:
        cmd += ["-e", int(exhaustiveness)]
    return cmd


def ensure_gnn_model(work: Path, log) -> Optional[Path]:
    """The GNN checkpoint, fetched once into a cache shared by every bundle.

    `pandadock gnn download-model` is not used. It requests release v4.0.0,
    which publishes no assets, prints "Model not found at release URL" and
    exits 0, so its failure is invisible to a return-code check. The asset is
    published under v4.1.1, so the release list is asked which release actually
    carries it and the newest one wins.
    """
    for candidate in (MODEL_CACHE / MODEL_NAME,
                      work / "models" / MODEL_NAME,
                      Path.home() / ".pandadock" / MODEL_NAME,
                      Path("models") / MODEL_NAME):
        if candidate.exists() and candidate.stat().st_size > MIN_MODEL_BYTES:
            return candidate

    url = newest_model_url(log) or FALLBACK_URL
    MODEL_CACHE.mkdir(parents=True, exist_ok=True)
    dest = MODEL_CACHE / MODEL_NAME
    # Downloaded beside the destination and moved into place only when whole, so
    # an interrupted fetch cannot leave a truncated checkpoint that every later
    # run then finds in the cache and loads.
    partial = dest.with_suffix(".part")
    try:
        download(url, partial, log)
    except Exception as exc:                       # noqa: BLE001 - the fallback is the point
        partial.unlink(missing_ok=True)
        log(f"dock: the GNN checkpoint could not be fetched: {exc}")
        return None
    # Size read before the unlink, not after: stat() on a file this line has
    # just deleted raises FileNotFoundError, which would turn the graceful
    # fallback to the empirical scorer into a dead run three stages in.
    size = partial.stat().st_size
    if size < MIN_MODEL_BYTES:
        partial.unlink(missing_ok=True)
        log(f"dock: the download was only {size} bytes, which is not a model")
        return None
    partial.replace(dest)
    return dest


def newest_model_url(log) -> Optional[str]:
    """Ask the release list which release actually carries the checkpoint."""
    import json as _json
    import urllib.request

    try:
        with urllib.request.urlopen(RELEASES_API, timeout=30) as response:
            releases = _json.loads(response.read().decode("utf-8"))
    except Exception:                              # noqa: BLE001 - offline is normal here
        return None
    for release in releases if isinstance(releases, list) else []:
        for asset in release.get("assets") or []:
            if asset.get("name") == MODEL_NAME and asset.get("browser_download_url"):
                return asset["browser_download_url"]
    return None


def download(url: str, dest: Path, log) -> None:
    """Fetch with a progress bar, because 82 MB on a slow line looks like a hang."""
    import urllib.request

    from .console import bar_for

    with urllib.request.urlopen(url, timeout=120) as response:
        total = int(response.headers.get("Content-Length") or 0)
        with bar_for(log, "fetching the PandaDock GNN checkpoint",
                     total=total or None) as bar:
            done = 0
            with open(dest, "wb") as fh:
                while True:
                    chunk = response.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
                    done += len(chunk)
                    bar.update(done, note=f"{done / 1e6:.0f} MB")


def write_scores(out_dir: Path, dest: Path, log) -> list[dict]:
    """pose_id, score, GNN affinity and rank, from whatever PandaDock wrote.

    PandaDock's own JSON is preferred; the SD tags in poses.sdf are the fallback,
    because the JSON's exact filename has changed between versions and the tags
    have not.
    """
    rows: list[dict] = []
    for path in sorted(out_dir.glob("*_poses.json")) + sorted(out_dir.glob("*poses*.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
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
        rows = scores_from_hybrid_csv(out_dir / "hybrid_results.csv")
    if not rows:
        rows = scores_from_sdf(out_dir / "poses.sdf")

    with open(dest, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["pose_id", "score", "gnn_affinity", "rank"])
        writer.writeheader()
        writer.writerows(rows)
    return rows


def scores_from_hybrid_csv(path: Path) -> list[dict]:
    """`pandadock hybrid` writes hybrid_results.csv and no JSON at all.

    Its columns are its own: rank, gnn_pec50, gnn_energy, vina_energy,
    activity_prob. Nothing that reads `dock`'s output finds anything here, which
    is how a hybrid run came back with three poses and an empty score column,
    reported on the card as "best None kcal/mol".

    `score` is the GNN energy rather than Vina's, because in hybrid mode the GNN
    is what did the ranking and the score column should be the number the ranking
    was made on.
    """
    if not path.exists():
        return []
    rows = []
    try:
        with open(path, newline="", encoding="utf-8") as fh:
            for index, record in enumerate(csv.DictReader(fh), start=1):
                rows.append({
                    "pose_id": record.get("pose_id") or f"pose{index}",
                    "score": _number(record.get("gnn_energy") or record.get("vina_energy")),
                    "gnn_affinity": _number(record.get("gnn_pec50")),
                    # int, not the float _number returns: this is a position in
                    # a list and it is written straight into the CSV.
                    "rank": int(_number(record.get("rank")) or index),
                })
    except (OSError, csv.Error):
        return []
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
            # PandaDock's own tag names, both modes: score_kcal_per_mol and
            # energy_gnn_pec50 are what it actually writes, and neither was
            # being read.
            "score": _number(props.get("score", props.get("Score", props.get(
                "score_kcal_per_mol", props.get("energy_gnn_energy", props.get("energy")))))),
            "gnn_affinity": _number(props.get("gnn_score", props.get(
                "energy_gnn_pec50", props.get("predicted_affinity")))),
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
