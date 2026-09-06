"""Reading the bundle's affinity block for the scorecard.

All the computing happened on the machine with the GPU. This turns what came
back into the two rows the card shows, and into the sentences beside them.

It is deliberately tolerant. Three states arrive here and all three are normal:
a block that ran, a block that says why it did not, and no block at all, which
is every archive built before stage 5 existed. None of them is an error, and
none of them changes a single graded metric: the composite measures agreement
with an experimental structure, and a predicted affinity has no counterpart in
the run to agree with.
"""

from __future__ import annotations

from typing import Any, Optional


def summarise(block: dict) -> dict[str, Any]:
    """The panel's content: state, two columns, the change, and the caveats."""
    if not block:
        return {"state": "absent",
                "note": "This archive was built before the affinity stage existed, "
                        "so there is nothing to show."}
    if not block.get("requested", True):
        return {"state": "declined",
                "note": "This run did not ask for an affinity prediction."}
    if not block.get("ran"):
        return {"state": "skipped",
                "note": "The affinity stage could not run: "
                        f"{block.get('reason') or 'no reason was recorded'}."}

    pre = block.get("pre_md") or {}
    post = block.get("post_md") or {}
    delta = block.get("delta") or {}
    rows = [
        row("pIC50", pre.get("pic50"), post.get("pic50_mean"), post.get("pic50_sd"),
            delta.get("pic50"), higher_is_better=True),
        row("Boltz-2 value", pre.get("affinity_pred_value"),
            post.get("affinity_pred_value_mean"), post.get("affinity_pred_value_sd"),
            delta.get("affinity_pred_value"), higher_is_better=False),
        row("Binder probability", pre.get("affinity_probability_binary"),
            post.get("probability_binary_mean"), post.get("probability_binary_sd"),
            delta.get("probability_binary"), higher_is_better=True),
    ]
    frames = (block.get("frames") or {}).get("scored") or []
    return {
        "state": "ran",
        "rows": [r for r in rows if r is not None],
        "n_post_frames": max(0, len(frames) - 1),
        "spread": spread(post),
        "msa": block.get("msa") or {},
        "seconds": block.get("seconds"),
        "route": block.get("route"),
        "how": how(block),
        "unit": block.get("unit", ""),
        "verdict": verdict(pre, post, delta),
    }


def how(block: dict) -> str:
    """What was actually done, in the reader's words rather than the route's.

    The wording has to follow the route or it lies. On the forced route Boltz's
    diffusion still runs, constrained to stay within a threshold of the pose;
    on the unforced one it runs free and the pose only conditions the trunk.
    Those are different measurements and the panel should not describe one while
    reporting the other.
    """
    engine = block.get("engine") or {}
    threshold = engine.get("template_threshold_a")
    if block.get("route") == "boltz-forced-template" and threshold:
        return (f"Boltz-2's affinity head. The complex was held within {threshold:g} A of the "
                f"structure this pipeline produced while Boltz-2 rebuilt it, so the head reads "
                f"the docked and relaxed conformations rather than one of its own.")
    return ("Boltz-2's affinity head. The structure this pipeline produced conditioned the "
            "prediction but did not constrain it, so the coordinates scored are partly "
            "Boltz-2's own.")


def row(label: str, before, after, sd, change, higher_is_better: bool) -> Optional[dict]:
    if before is None and after is None:
        return None
    return {
        "label": label,
        "before": before,
        "after": after,
        "sd": sd or 0.0,
        "delta": change,
        # Which direction is good differs per row: the Boltz-2 value is a
        # log10(IC50), so lower is stronger, while pIC50 and the probability
        # both read the usual way round. Getting this wrong colours an
        # improvement red.
        "better": None if change in (None, 0) else (
            (change > 0) if higher_is_better else (change < 0)),
    }


def spread(post: dict) -> Optional[float]:
    """The full range across the scored frames, which the mean hides."""
    values = post.get("affinity_pred_value_range")
    if not values or len(values) != 2:
        return None
    return round(abs(values[1] - values[0]), 3)


def verdict(pre: dict, post: dict, delta: dict) -> str:
    """One sentence, and never one that claims the pose is right."""
    if not post:
        return "Only the docked pose was scored, so there is nothing to compare it against."
    mean = post.get("pic50_mean")
    if mean is None:
        return "The affinity head ran but returned no value."
    text = f"Predicted pIC50 {mean:.2f} after MD"
    change = (delta or {}).get("pic50")
    if change is not None:
        direction = ("stronger" if change > 0 else "weaker" if change < 0 else "unchanged")
        text += f", {abs(change):.2f} log units {direction} than the docked pose"
    width = spread(post)
    if width and width >= 1.0:
        text += (f". The scored frames disagree by {width:.1f} log units, so the pose is not "
                 f"settled and the mean is doing more work than it should")
    return text + "."
