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
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from gobsmacked_run import console as console_mod  # noqa: E402
from gobsmacked_run import dock, fold, md, noise, prep, schema, summarise  # noqa: E402

STAGES = {
    "fold": fold.run,
    "prep": prep.run,
    "dock": dock.run,
    "md": md.run,
    "summarise": summarise.run,
}

# What each stage is for, in one line, printed under its heading so a reader who
# has never seen this pipeline knows what is taking the time.
BLURB = {
    "fold": "ESMFold, unless the bundle already carries a model",
    "prep": "trim, protonate at the campaign pH, and build the ligand conformer",
    "dock": "PandaDock inside the campaign's box",
    "md": "solvate, minimise, equilibrate and run",
    "summarise": "turn the trajectory into numbers and pack the archive",
}

# Every number below was measured, on one run: EGFR plus erlotinib, 253
# residues, a 58,266-atom box, OpenCL on an M1 Max. They are priors, not
# promises, and the bars say so when a stage runs past its estimate.
#
# The flat table these replaced said MD took 14 minutes. It took 54, for half
# the production length the default campaign asks for, so the estimate was out
# by a factor of eight on the one stage where being told the truth matters:
# nobody abandons a 2-minute stage, and everybody eyes a 70-minute one.
MD_SECONDS_PER_NS = 3500.0      # 500 ps of production in 1,752 s, so 24.7 ns/day
MD_SETUP_SECONDS = 450.0        # 290 s solvating, 23 s minimising, the rest imports
SUMMARISE_SECONDS_PER_FRAME = 1.1   # 100 frames in 110 s, pocket volume dominating
FOLD_SECONDS = 180.0            # ESMFold on a ~250-residue chain; the one number not measured here
PREP_SECONDS = 10.0             # measured at 3-5 s, rounded up for a cold RDKit


def estimate_seconds(name: str, campaign: dict, skip_fold: bool) -> float:
    """How long a stage will take on this campaign, not on a typical one.

    MD is the stage worth computing rather than tabulating: it scales directly
    with the production length the campaign asks for, and that is the number the
    person waiting has already chosen. A campaign asking for 1 ns is an hour
    longer than one asking for 500 ps, and a single table entry cannot say so.
    """
    md_cfg = campaign.get("md") or {}
    if name == "fold":
        return 0.0 if skip_fold else FOLD_SECONDS
    if name == "prep":
        return PREP_SECONDS
    if name == "dock":
        return dock.estimate_seconds(campaign.get("docking") or {})
    if name == "md":
        nanoseconds = (float(md_cfg.get("equilibration_ps", 100))
                       + float(md_cfg.get("production_ps", 1000))) / 1000.0
        return MD_SETUP_SECONDS + nanoseconds * MD_SECONDS_PER_NS
    if name == "summarise":
        interval = max(1.0, float(md_cfg.get("frame_interval_ps", 10)))
        frames = float(md_cfg.get("production_ps", 1000)) / interval
        return 20.0 + frames * SUMMARISE_SECONDS_PER_FRAME
    return 60.0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Run a GOBSMACKED campaign.")
    parser.add_argument("--stage", choices=list(STAGES), help="rerun from this stage onward")
    parser.add_argument("--only", choices=list(STAGES), help="run just this stage")
    parser.add_argument("--list", action="store_true", help="show the plan and exit")
    parser.add_argument("--campaign", default="campaign.yaml")
    parser.add_argument("--no-colour", "--no-color", action="store_true", dest="no_colour",
                        help="plain output, no escape sequences")
    args = parser.parse_args(argv)

    campaign_path = HERE / args.campaign
    if not campaign_path.exists():
        print(f"No {args.campaign} beside run.py. Unpack the bundle and run from inside it.")
        return 2
    campaign = yaml.safe_load(campaign_path.read_text(encoding="utf-8")) or {}
    job_id = campaign.get("job_id", "unknown")

    work = HERE / "work"
    results = HERE / "results"
    for directory in (work, results, results / "logs", results / "traj", results / "poses"):
        directory.mkdir(parents=True, exist_ok=True)

    log = console_mod.Console(log_path=results / "logs" / "run.log")
    if args.no_colour:
        log.colour = False
    noise.install(log)

    plan = build_plan(args, work)
    skipped_fold = "fold" in plan and (HERE / "model_apo.pdb").exists()

    def estimate(name: str) -> float:
        return estimate_seconds(name, campaign, skipped_fold)

    log.banner(job_id, campaign.get("title") or "")
    rows = []
    for name in STAGES:
        if name in plan:
            note = ("a model was supplied, nothing to fold" if name == "fold" and skipped_fold
                    else f"{BLURB[name]}  ~{console_mod.human(estimate(name))}")
            rows.append((name, "run", note))
        elif done_marker(work, name).exists():
            rows.append((name, "done", "already done, skipping"))
        else:
            rows.append((name, "skip", "not in this run"))
    minutes = sum(estimate(name) for name in plan) / 60.0
    log.plan(rows, minutes)

    if args.list:
        return 0
    if not plan:
        log.write("  Everything is already done. Use --stage to rerun from a stage.")
        return 0

    # campaign.yaml is echoed into the results so the server sees exactly what
    # ran, including any edit made here between generating and running.
    (results / "campaign.yaml").write_text(campaign_path.read_text(encoding="utf-8"), encoding="utf-8")

    timings: dict[str, float] = {}
    warnings: list[str] = []
    for index, name in enumerate(plan, start=1):
        started = time.time()
        log.stage_start(index, len(plan), name, BLURB[name], estimate(name) or None)
        try:
            outcome = STAGES[name](campaign, work, results, log) or {}
        except Exception as exc:
            log.fail(f"{name} failed: {exc}")
            # The traceback goes to the log file in full and to the terminal as
            # the last three frames. The whole thing on screen buries the one
            # line that says what to do next, and the file has it either way.
            frames = traceback.format_exc().rstrip().splitlines()
            for line in frames[-6:]:
                log.detail(line, logged=line)
            with open(results / "logs" / "run.log", "a",
                      encoding="utf-8", errors="replace") as fh:
                fh.write(traceback.format_exc())
            log.warn(f"Fix the cause and rerun with --stage {name}; "
                     f"earlier stages are already done.")
            return 1
        timings[name] = time.time() - started
        warnings.extend(outcome.get("warnings") or [])
        done_marker(work, name).write_text(json.dumps(
            {"finished": schema.now(), "seconds": round(timings[name], 1), **{
                k: v for k, v in outcome.items() if k != "warnings"}}, indent=2),
            encoding="utf-8")
        log.stage_end(name, timings[name], outcome.get("headline", ""))

    # A partial rerun (--stage, --only) would otherwise report only the stages
    # it ran, so the archive from a resumed run would claim docking took no time
    # at all. The .done markers hold what the earlier stages actually cost.
    for name in STAGES:
        if name in timings:
            continue
        marker = done_marker(work, name)
        if marker.exists():
            try:
                timings[name] = float(json.loads(marker.read_text(encoding="utf-8")).get("seconds", 0.0))
            except (json.JSONDecodeError, TypeError, ValueError):
                continue

    schema.write_manifest(results, job_id, campaign_path, timings, warnings)
    missing = schema.check_complete(results)
    if missing:
        log.fail("The run finished but the archive would be incomplete, missing: "
                 + ", ".join(missing))
        log.warn("Rerun the stage that writes them rather than uploading this.")
        return 1

    # Top level, beside run.py: see summarise.pack for why not inside results/.
    archive = summarise.pack(results, HERE / "results.tar.gz", log)
    ordered = {name: timings[name] for name in STAGES if name in timings}
    log.summary(ordered, warnings, archive.relative_to(HERE), job_id)
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
