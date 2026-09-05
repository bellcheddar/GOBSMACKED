"""Application factory.

Four blueprints, one SQLite file, no background workers: everything the droplet
does is inside a request, and everything that would need a GPU happens in the
bundle on the user's own machine.
"""

from __future__ import annotations

import os
from pathlib import Path

from flask import Flask, render_template
from werkzeug.middleware.proxy_fix import ProxyFix

from . import config, db


def create_app(test_config: dict | None = None) -> Flask:
    app = Flask(__name__, instance_relative_config=False)
    app.config.update(
        SECRET_KEY=os.environ.get("SECRET_KEY", "gobsmacked-dev-key"),
        MAX_CONTENT_LENGTH=config.MAX_UPLOAD_BYTES,
        JSON_SORT_KEYS=False,
        VERSION=config.VERSION,
    )
    if test_config:
        app.config.update(test_config)

    # nginx terminates TLS and proxies to gunicorn over http, so without this
    # every url_for(..., _external=True) says http and the wrong host. That
    # matters now that the app hands out a curl command people paste.
    app.wsgi_app = ProxyFix(app.wsgi_app, x_for=1, x_proto=1, x_host=1)

    config.ensure_dirs()
    db.init_db()

    from .routes import prepare, analyze, runs, about
    app.register_blueprint(prepare.bp)
    app.register_blueprint(analyze.bp)
    app.register_blueprint(runs.bp)
    app.register_blueprint(about.bp)

    @app.route("/")
    def index():
        return prepare.prepare_page()

    @app.route("/healthz")
    def healthz():
        return {"ok": True, "version": config.VERSION}

    @app.errorhandler(404)
    def not_found(_e):
        return render_template("error.html", message="No such page.",
                               fix="Check the address, or start again from Prepare."), 404

    @app.errorhandler(413)
    def too_large(_e):
        limit_mb = config.MAX_UPLOAD_BYTES // (1024 * 1024)
        return render_template(
            "error.html",
            message=f"That upload is larger than the {limit_mb} MB limit.",
            fix="Re-run the bundle's summarise stage with a longer frame interval, "
                "or trim the trajectory before packing it.",
        ), 413

    @app.errorhandler(500)
    def server_error(_e):
        return render_template("error.html", message="Something failed on the server.",
                               fix="Try again; if it repeats, the run log is in journalctl."), 500

    # Two template helpers the results page leans on. They live here rather
    # than in the template so the arc geometry and the LED rule are testable.
    @app.template_global()
    def dial_arc(score: float) -> str:
        """SVG path for the score dial's filled arc.

        The track runs from (20,110) to (180,110) on a radius-80 semicircle
        centred at (100,110). A score of 100 would put the end point exactly on
        the start point, which SVG draws as nothing at all, so the fill stops a
        hair short of the full sweep.
        """
        import math

        fraction = max(0.0, min(float(score or 0), 99.9)) / 100.0
        theta = math.radians(180.0 - 180.0 * fraction)
        x = 100 + 80 * math.cos(theta)
        y = 110 - 80 * math.sin(theta)
        return f"M20 110 A80 80 0 0 1 {x:.1f} {y:.1f}"

    @app.template_global()
    def led(predicted, reference) -> str:
        """Switch-list LED class: green when the prediction matches the crystal,
        amber when they differ, grey when neither has the feature."""
        if predicted is None and reference is None:
            return "off"
        if reference is None:
            return "" if predicted else "off"
        if predicted == reference:
            return "" if predicted else "off"
        return "amber"

    # Static assets are cached hard by nginx, so every reference carries a
    # ?v=<mtime> stamp. Without it a returning visitor keeps the old CSS after
    # a deploy and the change is invisible.
    @app.context_processor
    def asset_helper():
        def asset(path: str) -> str:
            full = Path(app.static_folder or "") / path
            stamp = int(full.stat().st_mtime) if full.exists() else 0
            return f"/static/{path}?v={stamp}"
        return {"asset": asset, "version": config.VERSION}

    # Jinja's default HTML template caching is fine, but the response needs an
    # explicit Cache-Control: Flask sends none for a rendered template, and a
    # browser's heuristic caching then pins the old ?v= asset URLs.
    @app.after_request
    def no_html_cache(resp):
        if resp.mimetype == "text/html":
            resp.headers.setdefault("Cache-Control", "no-cache")
        return resp

    return app
