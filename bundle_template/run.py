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
import os
import sys
import time
import traceback
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from gobsmacked_run import console as console_mod  # noqa: E402
from gobsmacked_run import affinity, dock, fold, md, noise, prep, schema, summarise  # noqa: E402

STAGES = {
    "fold": fold.run,
    "prep": prep.run,
    "dock": dock.run,
    "md": md.run,
    "affinity": affinity.run,
    "summarise": summarise.run,
}

# What each stage is for, in one line, printed under its heading so a reader who
# has never seen this pipeline knows what is taking the time.
BLURB = {
    # fold's line depends on the campaign, so blurb_for() overrides this one.
    # Left here so every stage still has an entry and the dict reads completely.
    "fold": "ESMFold, unless the bundle already carries a model",
    "prep": "trim, protonate at the campaign pH, and build the ligand conformer",
    "dock": "PandaDock inside the campaign's box",
    "md": "solvate, minimise, equilibrate and run",
    "affinity": "score the docked and the relaxed pose with Boltz-2",
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
# Co-folding is a Boltz-2 pass with the structure module doing the work, not a
# template being steered, so it is nothing like ESMFold's three minutes. Taken
# from the affinity stage's own per-pose timings on the same size of target,
# where a pose took 6 to 17 minutes; this is one pass and the MSA is shared with
# affinity rather than paid for twice. Quoting 3 minutes for it was wrong by
# roughly a factor of five.
COFOLD_SECONDS = 900.0
PREP_SECONDS = 10.0             # measured at 3-5 s, rounded up for a cold RDKit


def blurb_for(name: str, campaign: dict) -> str:
    """What a stage is about to do, for this campaign rather than in general.

    Only fold differs, and it differs in a way the reader needs: co-folding
    builds the pocket around the ligand and then throws the ligand away, which
    is a different thing from folding a sequence and is worth saying before it
    takes a quarter of an hour.
    """
    if name == "fold" and fold.method_of(campaign) == "boltz2":
        return "co-fold the protein with the ligand (Boltz-2), then keep the protein"
    return BLURB[name]


def estimate_seconds(name: str, campaign: dict, skip_fold: bool) -> float:
    """How long a stage will take on this campaign, not on a typical one.

    MD is the stage worth computing rather than tabulating: it scales directly
    with the production length the campaign asks for, and that is the number the
    person waiting has already chosen. A campaign asking for 1 ns is an hour
    longer than one asking for 500 ps, and a single table entry cannot say so.
    """
    md_cfg = campaign.get("md") or {}
    if name == "fold":
        if skip_fold:
            return 0.0
        return COFOLD_SECONDS if fold.method_of(campaign) == "boltz2" else FOLD_SECONDS
    if name == "prep":
        return PREP_SECONDS
    if name == "dock":
        return dock.estimate_seconds(campaign.get("docking") or {})
    if name == "md":
        nanoseconds = (float(md_cfg.get("equilibration_ps", 100))
                       + float(md_cfg.get("production_ps", 1000))) / 1000.0
        return MD_SETUP_SECONDS + nanoseconds * MD_SECONDS_PER_NS
    if name == "affinity":
        if not (campaign.get("affinity") or {}).get("include", True):
            return 0.0
        # The MSA dominates a target that has not been seen before and is free
        # on one that has; the trunk is one pass per pose and the head is cheap.
        cfg = campaign.get("affinity") or {}
        poses = 2 if (cfg.get("frames") or "cluster") == "single" else \
            int(cfg.get("n_frames", 5) or 5) + 1
        return 240.0 + 30.0 * poses
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
                    else f"{blurb_for(name, campaign)}  ~{console_mod.human(estimate(name))}")
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
    # ran, including any edit made here between generating and running. Written
    # now so a run that dies mid-stage still carries it, and written again at the
    # end from the in-memory campaign, because a stage can correct it.
    (results / "campaign.yaml").write_text(campaign_path.read_text(encoding="utf-8"), encoding="utf-8")
    original_pocket = dict((campaign.get("pocket") or {}))

    # Applied here, from disk, on EVERY run rather than inside the fold stage.
    #
    # The re-centring used to live in fold.cofold, which is skipped on a resume
    # because fold.done exists. So a co-folded run that failed in dock and
    # retried came back with the campaign's original centre and docked 36 A from
    # its own receptor: the retry silently undid the fix the first attempt had
    # applied, and looked like a normal run while doing it.
    #
    # The co-folded ligand is already on disk, so the correction can be derived
    # again instead of remembered. Idempotent, and it does not care which stages
    # are being rerun.
    _apply_cofold_box(campaign, results, log)

    timings: dict[str, float] = {}
    warnings: list[str] = []
    for index, name in enumerate(plan, start=1):
        started = time.time()
        log.stage_start(index, len(plan), name, blurb_for(name, campaign), estimate(name) or None)
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

    # BEFORE the manifest and before pack. This block used to sit after
    # summarise.pack, so it rewrote a campaign.yaml that had already been sealed
    # into the archive: two co-folded runs graded F 40.0 on "ligand inside the
    # docking box" with the corrected centre nowhere in the tarball. Anything
    # that edits the results directory has to happen before the directory is
    # packed, and the warning has to exist before the manifest records warnings.
    #
    # Co-folding is the case that exists: it moves the docking box into the
    # frame of the receptor it just built, 36 A on a real run. Prepare computed
    # that box in a different structure's frame, so the server's validity check
    # measured against a centre nowhere near the pose and capped a run whose
    # pose was in exactly the right place.
    edits = _campaign_edits(original_pocket, campaign.get("pocket") or {})
    if edits:
        (results / "campaign.yaml").write_text(
            yaml.safe_dump(campaign, sort_keys=False), encoding="utf-8")
        warnings.extend(edits)
        for edit in edits:
            log.detail(edit)

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


def _apply_cofold_box(campaign: dict, results: Path, log) -> None:
    """Move the docking box into the co-folded receptor's frame, if there is one.

    Does nothing at all when the run did not co-fold, which is every run that
    starts from a supplied or ESMFolded structure.
    """
    ligand = results / "cofold" / "cofold_ligand.sdf"
    if not ligand.exists():
        return
    centre = fold._ligand_centre(ligand)
    if centre is None:
        return
    pocket = campaign.setdefault("pocket", {})
    was = list(pocket.get("center") or [])
    if was and fold._distance(was, centre) is not None and fold._distance(was, centre) < 0.5:
        return                                  # already in this frame
    pocket["center"] = centre
    moved = fold._distance(was, centre)
    log.detail(f"box centre taken from the co-folded ligand: "
               f"{[round(v, 2) for v in centre]}"
               + (f", {moved:.1f} A from the campaign's" if moved is not None else ""))


def _campaign_edits(before: dict, after: dict) -> list[str]:
    """What a stage changed about the pocket, in words, for the archive.

    Only reports a move worth mentioning: a re-centre of a fraction of an
    Angstrom is rounding, not a decision anyone needs to read about.
    """
    a, b = before.get("center"), after.get("center")
    if not (a and b and len(a) == 3 and len(b) == 3):
        return []
    moved = sum((float(x) - float(y)) ** 2 for x, y in zip(a, b)) ** 0.5
    if moved < 0.5:
        return []
    return [f"The docking box was re-centred by {moved:.1f} A during the run, from "
            f"{[round(float(v), 2) for v in a]} to {[round(float(v), 2) for v in b]}, "
            f"to put it in the frame of the receptor this run built. The campaign "
            f"recorded here is the one that ran."]


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
    code = main()
    # os._exit, not sys.exit. sys.exit waits for every non-daemon thread and runs
    # every atexit handler, and this process has imported torch (via pandadock)
    # and OpenMM at module scope whether or not their stages ran. On a real run
    # that wait was 25 minutes: run.py printed its final line, returned 1, and
    # the interpreter then sat there until a watchdog killed it, which reads from
    # outside as a hung stage rather than a finished one.
    #
    # Nothing is lost by leaving early. Every stage writes its own .done marker,
    # the archive is packed and closed, and run.log is opened and closed per
    # write, so the only buffers that could still hold anything are these two.
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.flush()
        except Exception:                             # noqa: BLE001
            pass
    os._exit(code)
