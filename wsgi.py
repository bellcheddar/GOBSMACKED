"""WSGI entry point. `gunicorn wsgi:app` in the systemd unit."""

from app import create_app

app = create_app()
