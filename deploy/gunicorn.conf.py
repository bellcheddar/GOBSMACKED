"""Gunicorn config for the GOBSMACKED web service.

Two workers, not three: the droplet has 3.8 GB and an analysis holds a
trajectory, three PLIP subprocesses in sequence and a PandaMap render. Two
concurrent analyses is the most this box should attempt.

The timeout is long because the upload endpoint does the whole analysis inline
(ingest, superpose, PLIP three times, PandaMap, modes, dynamics) rather than
handing it to a queue. Measured at 20 to 40 s for a 20-frame archive.
"""
import os

bind = os.environ.get("BIND_ADDR", "127.0.0.1:8009")
workers = int(os.environ.get("WEB_WORKERS", "2"))
worker_class = "sync"
timeout = int(os.environ.get("WEB_TIMEOUT", "300"))
graceful_timeout = 60
keepalive = 5
# Log to stdout/stderr so journald captures it (journalctl -u gobsmacked-web).
accesslog = "-"
errorlog = "-"
loglevel = os.environ.get("LOG_LEVEL", "info")
proc_name = "gobsmacked-web"
