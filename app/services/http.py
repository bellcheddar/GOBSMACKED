"""One retrying HTTP session for every external service.

UniProt, InterPro, KLIFS, GPCRdb, RCSB, AlphaFold DB and the ESM Atlas are all
public and all occasionally return a transient 5xx or a rate-limit. A single
session gives connection reuse; tenacity gives polite backoff. Adapted from
AlphaFraud, where the same three-line shape has run weekly for months.
"""

from __future__ import annotations

import os
import threading
from pathlib import Path
from typing import Any, Optional

import requests
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from .. import config

_SESSION = requests.Session()
_SESSION.headers.update({"User-Agent": config.USER_AGENT})

_RETRYABLE_STATUS = {429, 500, 502, 503, 504}


class RetryableStatus(Exception):
    """Raised for a retryable HTTP status so tenacity re-attempts the request."""


_retry = retry(
    reraise=True,
    stop=stop_after_attempt(config.HTTP_MAX_RETRIES),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    retry=retry_if_exception_type((requests.ConnectionError, requests.Timeout, RetryableStatus)),
)


@_retry
def get(url: str, **kwargs) -> requests.Response:
    kwargs.setdefault("timeout", config.HTTP_TIMEOUT)
    resp = _SESSION.get(url, **kwargs)
    if resp.status_code in _RETRYABLE_STATUS:
        raise RetryableStatus(f"{resp.status_code} for GET {url}")
    return resp


@_retry
def post(url: str, **kwargs) -> requests.Response:
    kwargs.setdefault("timeout", config.HTTP_TIMEOUT)
    resp = _SESSION.post(url, **kwargs)
    if resp.status_code in _RETRYABLE_STATUS:
        raise RetryableStatus(f"{resp.status_code} for POST {url}")
    return resp


def get_json(url: str, **kwargs) -> Optional[Any]:
    """GET returning parsed JSON, or None on 404 (a miss, not a failure)."""
    resp = get(url, **kwargs)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    return resp.json()


def post_json(url: str, payload: dict[str, Any], **kwargs) -> Optional[Any]:
    resp = post(url, json=payload, **kwargs)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    return resp.json()


def download(url: str, dest: Path, skip_if_exists: bool = True) -> Optional[Path]:
    """Stream a URL to `dest`. Returns dest, or None on 404."""
    dest = Path(dest)
    if skip_if_exists and dest.exists() and dest.stat().st_size > 0:
        return dest
    resp = get(url, stream=True)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Unique temp name per writer so two workers fetching the same structure
    # never clobber each other's partial write.
    tmp = dest.with_suffix(f"{dest.suffix}.{os.getpid()}.{threading.get_ident()}.part")
    try:
        with open(tmp, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 16):
                fh.write(chunk)
        tmp.replace(dest)
    finally:
        if tmp.exists():
            tmp.unlink(missing_ok=True)
    return dest
