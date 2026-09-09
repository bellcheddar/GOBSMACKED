"""A second opinion on which pose should have been ranked first.

The docking engine searches and ranks with one scoring function. Measured on
five starting structures for EGFR plus erlotinib, the search was reliably better
than the ranking: a near-native pose was in the list every time, and on the
cleanest possible control -- docking a ligand back into the crystal it came from
-- the 1.51 A pose was ranked EIGHTH of ten, behind an 8.77 A decoy. Rescoring
those same ten coordinates with Vinardo put it first.

So this stage re-ranks what docking already produced. It samples nothing, moves
nothing and changes no output the rest of the pipeline reads: pose 1 remains
pose 1 and the scorecard still grades it. What it adds is a second column, and a
sentence when the two functions disagree about which pose deserved the top slot.

It now decides which pose goes forward, and the reason it took a while to get
there is worth recording. The first version reported the disagreement and did
not act on it, on the grounds that Vinardo lost on the one receptor that
worked. That comparison was wrong: it came from a re-DOCKING sweep, not from
re-ranking the shipped poses. Measured properly, across nine pose sets whose
distance to the crystal is known:

    engine   1/9 within 2 A, mean 7.02 A
    Vinardo  4/9 within 2 A, mean 4.03 A

and on the run the engine got right, Vinardo picks the same pose. There was no
case where promoting it cost anything, which removed the only argument against.

Consensus was tested at the same time and dropped: a three-way Borda count and a
"Vinardo unless the other two overrule" rule both scored 4/9 and a mean of 4.03,
identical to Vinardo alone to two decimals. The machinery would have bought
nothing.

One target and one ligand, so `docking.rank_by: engine` turns it off.
"""

from __future__ import annotations

import csv
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Optional

from .console import bar_for

# Vinardo rather than Vina: on the fixed pose sets it had the lowest mean top-1
# RMSD of the seven functions tried (3.50 A, against an oracle of 3.09 A for the
# best pose present), and it separated the right pose from the runner-up by
# 0.41 kcal/mol where Vina managed 0.10.
FUNCTION = "vinardo"
TIMEOUT_S = 600


def available() -> bool:
    return shutil.which("pixi") is not None


def run(receptor: Path, poses: Path, dest: Path, log) -> dict[str, Any]:
    """Score every pose in `poses` against `receptor`; never raise.

    A second opinion is a decoration. If it cannot be produced the run carries
    on with the ranking it already has, because the alternative is losing an
    hour of docking and MD to a nicety.
    """
    if not (receptor.exists() and poses.exists()):
        return {"ran": False, "reason": "no receptor or poses to score"}
    if not available():
        return {"ran": False, "reason": "pixi is not on PATH"}

    cmd = ["pixi", "run", "-e", "rescore", "smina",
           "-r", str(receptor), "-l", str(poses),
           "--score_only", "--scoring", FUNCTION]
    try:
        with bar_for(log, f"re-ranking the poses with {FUNCTION}"):
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  encoding="utf-8", errors="replace", timeout=TIMEOUT_S)
    except (subprocess.TimeoutExpired, OSError) as exc:
        return {"ran": False, "reason": f"{type(exc).__name__}"}
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-2:]
        return {"ran": False, "reason": " / ".join(t.strip() for t in tail) or "non-zero exit"}

    scores = [float(m) for m in re.findall(r"^Affinity:\s+(-?\d+\.?\d*)", proc.stdout, re.M)]
    if not scores:
        return {"ran": False, "reason": "no affinities were parsed from the output"}

    order = sorted(range(len(scores)), key=lambda i: scores[i])
    block = {
        "ran": True,
        "function": FUNCTION,
        "scores": [round(v, 3) for v in scores],
        # 1-based, to match the pose numbering everything else uses.
        "ranking": [i + 1 for i in order],
        "top_pose": order[0] + 1,
        "agrees": order[0] == 0,
    }
    write(dest, block)
    if block["agrees"]:
        log(f"rescore: {FUNCTION} agrees, pose 1 is its best of {len(scores)}")
    else:
        log(f"rescore: {FUNCTION} would rank pose {block['top_pose']} first, not pose 1")
    return block


def write(dest: Path, block: dict) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["pose_id", f"{block.get('function', FUNCTION)}_kcal_per_mol", "rank"])
        scores = block.get("scores") or []
        ranks = {pose: i + 1 for i, pose in enumerate(block.get("ranking") or [])}
        for index, value in enumerate(scores, start=1):
            writer.writerow([f"pose{index}", value, ranks.get(index, "")])
