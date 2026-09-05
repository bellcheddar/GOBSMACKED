"""Runs: the table of every prepared run, and the per-run page.

Visibility is the whole design here. A private run is absent from listings
rather than greyed out, because a greyed row leaks that it exists and how many
there are. Ownership is proved by a 32-character token issued at Prepare time
and stored only as a sha256, sent by the browser in `X-Owner-Token` (it is kept
in localStorage, so the owner never types it after the first time).
"""

from __future__ import annotations

import json
from pathlib import Path

from flask import (Blueprint, abort, jsonify, redirect, render_template,
                   request, send_file, url_for)

from .. import config, db

bp = Blueprint("runs", __name__)

STAGE_ORDER = ["Fetch", "Annotate", "Fold", "Dock", "MD", "Verify", "Mode"]


def owner_token_from_request() -> str | None:
    return (request.headers.get("X-Owner-Token")
            or request.args.get("token")
            or (request.get_json(silent=True) or {}).get("token"))


def stage_states(row) -> list[str]:
    """The seven-cell ministrip for one run, in the rail's own vocabulary."""
    status = row["status"]
    if status == "prepared":
        return ["ready", "ready", "pending", "pending", "pending", "pending", "pending"]
    if status == "results_uploaded":
        return ["ready"] * 5 + ["pending", "pending"]
    if status == "failed":
        return ["ready", "ready", "failed", "failed", "failed", "failed", "failed"]
    verify = "ready" if row["reference_pdb"] else "optional"
    mode = "ready" if row["mode_match"] else ("optional" if row["mode_predicted"] else "pending")
    return ["ready"] * 5 + [verify, mode]


@bp.route("/runs")
def runs_page():
    counts = db.status_counts()
    total = sum(counts.values())
    stages = [
        {"name": "Prepared", "state": "ready" if counts["prepared"] else "pending",
         "text": f"{counts['prepared']} bundles written, not yet run"},
        {"name": "Uploaded", "state": "ready" if counts["results_uploaded"] else "pending",
         "text": f"{counts['results_uploaded']} archives in, not yet scored"},
        {"name": "Analysed", "state": "ready" if counts["analysed"] else "pending",
         "text": f"{counts['analysed']} runs scored"},
        {"name": "Failed", "state": "failed" if counts["failed"] else "pending",
         "text": (f"{counts['failed']} runs could not be analysed" if counts["failed"]
                  else "none")},
        {"name": "Kinase", "state": "optional", "text": "KLIFS pocket numbering"},
        {"name": "GPCR", "state": "optional", "text": "GPCRdb generic numbering"},
        {"name": "Other", "state": "optional", "text": "UniProt features only"},
    ]
    return render_template("runs.html", tab="runs", stages=stages)


@bp.get("/api/runs")
def api_runs():
    rows = db.list_jobs(owner_token=owner_token_from_request())
    for r in rows:
        r["stages"] = stage_states(r)
        r["url"] = url_for("runs.run_page", job_id=r["job_id"])
    return jsonify({"runs": rows, "counts": db.status_counts()})


@bp.route("/runs/<job_id>")
def run_page(job_id: str):
    row = db.get_job(job_id)
    if row is None:
        abort(404)
    token = owner_token_from_request()
    owned = db.token_matches(row, token)
    if row["visibility"] == "private" and not owned:
        # The URL is unguessable, but the page still asks: a link forwarded to
        # someone else must not hand over the run.
        return render_template("run_private.html", tab="runs", job_id=job_id), 403

    scorecard = json.loads(row["scorecard_json"]) if row["scorecard_json"] else None
    from .analyze import render_results
    return render_results(row, scorecard, owned=owned)


@bp.patch("/api/runs/<job_id>/visibility")
def api_visibility(job_id: str):
    row = db.get_job(job_id)
    if row is None:
        return jsonify({"error": "No such run."}), 404
    payload = request.get_json(silent=True) or {}
    if not db.token_matches(row, payload.get("token") or owner_token_from_request()):
        return jsonify({"error": "That owner key does not match this run."}), 403
    wanted = payload.get("visibility")
    if wanted not in ("public", "private"):
        return jsonify({"error": "visibility must be public or private."}), 400
    db.update_job(job_id, visibility=wanted)
    return jsonify({"job_id": job_id, "visibility": wanted})


@bp.patch("/api/runs/<job_id>/title")
def api_title(job_id: str):
    row = db.get_job(job_id)
    if row is None:
        return jsonify({"error": "No such run."}), 404
    payload = request.get_json(silent=True) or {}
    if not db.token_matches(row, payload.get("token") or owner_token_from_request()):
        return jsonify({"error": "That owner key does not match this run."}), 403
    db.update_job(job_id, title=(payload.get("title") or "").strip()[:120])
    return jsonify({"job_id": job_id, "title": payload.get("title")})


@bp.delete("/api/runs/<job_id>")
def api_delete(job_id: str):
    row = db.get_job(job_id)
    if row is None:
        return jsonify({"error": "No such run."}), 404
    payload = request.get_json(silent=True) or {}
    if not db.token_matches(row, payload.get("token") or owner_token_from_request()):
        return jsonify({"error": "That owner key does not match this run."}), 403
    if payload.get("confirm") != job_id:
        return jsonify({"error": "Retype the job ID to confirm deletion."}), 400
    import shutil
    run_dir = config.RUNS_DIR / job_id
    if run_dir.exists():
        shutil.rmtree(run_dir, ignore_errors=True)
    db.delete_job(job_id)
    return jsonify({"deleted": job_id})


@bp.get("/runs/<job_id>/report")
def report(job_id: str):
    """One printable page. `?download=1` sends it as a file instead."""
    row = _guarded(job_id)
    if not row["scorecard_json"]:
        abort(404)
    card = json.loads(row["scorecard_json"])
    html = render_template("report.html", tab="runs", row=row, card=card,
                           job_id=job_id, generated=db.now_iso())
    if request.args.get("download"):
        from flask import Response
        return Response(html, mimetype="text/html", headers={
            "Content-Disposition": f'attachment; filename="gobsmacked_{job_id}.html"'})
    return html


# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------

def _guarded(job_id: str):
    row = db.get_job(job_id)
    if row is None:
        abort(404)
    if row["visibility"] == "private" and not db.token_matches(row, owner_token_from_request()):
        abort(403)
    return row


@bp.get("/runs/<job_id>/bundle")
def download_bundle(job_id: str):
    _guarded(job_id)
    path = config.RUNS_DIR / job_id / f"run_bundle_{job_id}.tar.gz"
    if not path.exists():
        abort(404)
    return send_file(path, as_attachment=True, download_name=path.name)


@bp.get("/runs/<job_id>/results")
def download_results(job_id: str):
    _guarded(job_id)
    path = config.RUNS_DIR / job_id / "results.tar.gz"
    if not path.exists():
        abort(404)
    return send_file(path, as_attachment=True, download_name=f"results_{job_id}.tar.gz")


@bp.get("/runs/<job_id>/file/<path:name>")
def run_file(job_id: str, name: str):
    """Serve one extracted results file (a PDB for Mol*, a plot's JSON, a PNG)."""
    _guarded(job_id)
    base = (config.RUNS_DIR / job_id / "results").resolve()
    target = (base / name).resolve()
    if not str(target).startswith(str(base)) or not target.exists():
        abort(404)
    mimetypes = {".pdb": "chemical/x-pdb", ".cif": "chemical/x-cif",
                 ".json": "application/json", ".png": "image/png",
                 ".sdf": "chemical/x-mdl-sdfile", ".csv": "text/csv",
                 ".dcd": "application/octet-stream", ".log": "text/plain"}
    return send_file(target, mimetype=mimetypes.get(target.suffix, "application/octet-stream"))
