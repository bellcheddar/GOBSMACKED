#!/usr/bin/env python3
"""GOBSMACKED run bundle: fold, prep, dock, MD, summarise.

    pixi run gobsmacked                 # everything, resuming what is done
    pixi run gobsmacked --stage dock    # rerun from dock onward
    pixi run gobsmacked --list          # what would run, and what is already done

Each stage writes a `.done` marker in `work/` when it finishes. A rerun skips
completed stages, which matters because these stages are minutes to hours apart
in cost and an interrupted MD should not mean folding again.

Nothing here contacts the server that wrote the bundle. campaign.yaml goes in,
results.tar.gz comes out.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from gobsmacked_run import dock, fold, md, prep, schema, summarise  # noqa: E402

STAGES = {
    "fold": fold.run,
    "prep": prep.run,
    "dock": dock.run,
    "md": md.run,
    "summarise": summarise.run,
}

# Rough wall clock on one consumer GPU for a 300-residue domain and a drug-sized
# ligand, printed before the run so an hour-long MD is a decision rather than a
# surprise. Folding is the big variable: it is skipped whenever the bundle
# already carries a model, which is the usual case.
ESTIMATE_MIN = {"fold": 3, "prep": 1, "dock": 6, "md": 14, "summarise": 2}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Run a GOBSMACKED campaign.")
    parser.add_argument("--stage", choices=list(STAGES), help="rerun from this stage onward")
    parser.add_argument("--only", choices=list(STAGES), help="run just this stage")
    parser.add_argument("--list", action="store_true", help="show the plan and exit")
    parser.add_argument("--campaign", default="campaign.yaml")
    args = parser.parse_args(argv)

    campaign_path = HERE / args.campaign
    if not campaign_path.exists():
        print(f"No {args.campaign} beside run.py. Unpack the bundle and run from inside it.")
        return 2
    campaign = yaml.safe_load(campaign_path.read_text()) or {}
    job_id = campaign.get("job_id", "unknown")

    work = HERE / "work"
    results = HERE / "results"
    for directory in (work, results, results / "logs", results / "traj", results / "poses"):
        directory.mkdir(parents=True, exist_ok=True)
    log_path = results / "logs" / "run.log"

    def log(message: str) -> None:
        stamp = datetime.now(timezone.utc).strftime("%H:%M:%S")
        line = f"[{stamp}] {message}"
        print(line, flush=True)
        with open(log_path, "a") as fh:
            fh.write(line + "\n")

    plan = build_plan(args, work)
    if args.list:
        for name in STAGES:
            state = "run" if name in plan else ("done" if done_marker(work, name).exists() else "skip")
            print(f"  {name:<11} {state}")
        return 0

    # Folding is the big term and is skipped whenever the bundle carries a
    # model, which is the usual case. Counting it regardless made the printed
    # estimate wrong by three minutes on every ordinary run.
    minutes = sum(ESTIMATE_MIN[name] for name in plan
                  if not (name == "fold" and (HERE / "model_apo.pdb").exists()))
    log(f"GOBSMACKED {job_id}: {len(plan)} stages to run ({', '.join(plan) or 'nothing'}), "
        f"roughly {minutes} minutes on a consumer GPU")
    if not plan:
        log("Everything is already done. Use --stage to rerun from a stage.")
        return 0

    # campaign.yaml is echoed into the results so the server sees exactly what
    # ran, including any edit made here between generating and running.
    (results / "campaign.yaml").write_text(campaign_path.read_text())

    timings: dict[str, float] = {}
    warnings: list[str] = []
    for name in plan:
        started = time.time()
        log(f"--- {name} ---")
        try:
            outcome = STAGES[name](campaign, work, results, log) or {}
        except Exception as exc:
            log(f"{name} failed: {exc}")
            log(traceback.format_exc())
            log(f"Fix the cause and rerun with --stage {name}; earlier stages are already done.")
            return 1
        timings[name] = time.time() - started
        warnings.extend(outcome.get("warnings") or [])
        done_marker(work, name).write_text(json.dumps(
            {"finished": schema.now(), "seconds": round(timings[name], 1), **{
                k: v for k, v in outcome.items() if k != "warnings"}}, indent=2))
        log(f"{name}: {timings[name]:.0f} s")

    # A partial rerun (--stage, --only) would otherwise report only the stages
    # it ran, so the archive from a resumed run would claim docking took no time
    # at all. The .done markers hold what the earlier stages actually cost.
    for name in STAGES:
        if name in timings:
            continue
        marker = done_marker(work, name)
        if marker.exists():
            try:
                timings[name] = float(json.loads(marker.read_text()).get("seconds", 0.0))
            except (json.JSONDecodeError, TypeError, ValueError):
                continue

    schema.write_manifest(results, job_id, campaign_path, timings, warnings)
    missing = schema.check_complete(results)
    if missing:
        log("The run finished but the archive would be incomplete, missing: " + ", ".join(missing))
        log("Rerun the stage that writes them rather than uploading this.")
        return 1

    archive = summarise.pack(results, results / "results.tar.gz", log)
    log(f"Done. Upload {archive.relative_to(HERE)} on the Analyze tab.")
    for warning in warnings:
        log(f"warning: {warning}")
    return 0


def build_plan(args, work: Path) -> list[str]:
    names = list(STAGES)
    if args.only:
        return [args.only]
    if args.stage:
        return names[names.index(args.stage):]
    return [name for name in names if not done_marker(work, name).exists()]


def done_marker(work: Path, name: str) -> Path:
    return work / f"{name}.done"


if __name__ == "__main__":
    sys.exit(main())
