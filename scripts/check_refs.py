#!/usr/bin/env python3
"""Check every DOI in software.yaml against the DOI registry.

Deliberately not a request to the publisher's page: Oxford, ACS, Science and RSC
all answer a scripted HEAD with 403, so a check against the article URL reports
fifteen failures for fourteen perfectly good DOIs and one real one. Crossref's
API answers machine requests, returns the registered title, and is the actual
authority on whether a DOI exists.

The title is printed so a DOI that resolves to the wrong paper is visible too,
which the status code alone would never catch. `make check-refs` runs this.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import requests
import yaml

ROOT = Path(__file__).resolve().parent.parent
TIMEOUT = 30
# Crossref asks for a contact address in the User-Agent and gives politer
# service in return.
HEADERS = {"User-Agent": "GOBSMACKED reference checker (mailto:marc@marcdeller.com)"}
CROSSREF = "https://api.crossref.org/works/{doi}"
HANDLE = "https://doi.org/api/handles/{doi}"


def check(doi: str) -> tuple[bool, str]:
    try:
        resp = requests.get(CROSSREF.format(doi=doi), timeout=TIMEOUT, headers=HEADERS)
        if resp.status_code == 200:
            message = resp.json().get("message", {})
            title = (message.get("title") or ["(no title)"])[0]
            year = ((message.get("issued") or {}).get("date-parts") or [[None]])[0][0]
            return True, f"{title} ({year})"
        if resp.status_code == 404:
            # Not every DOI is registered with Crossref (DataCite, mEDRA and
            # others exist), so fall back to the handle system, which knows
            # every DOI regardless of registration agency.
            handle = requests.get(HANDLE.format(doi=doi), timeout=TIMEOUT, headers=HEADERS)
            if handle.status_code == 200 and handle.json().get("responseCode") == 1:
                return True, "registered, not in Crossref"
            return False, "not registered"
        return False, f"Crossref returned {resp.status_code}"
    except requests.RequestException as exc:
        return False, str(exc)


def main() -> int:
    data = yaml.safe_load((ROOT / "software.yaml").read_text()) or {}
    failures = []
    checked = 0
    for item in data.get("software", []):
        doi = item.get("doi")
        if not doi:
            continue
        checked += 1
        ok, detail = check(doi)
        print(f"{'ok  ' if ok else 'FAIL'}  {item['tool']:<28} {doi:<32} {detail[:80]}")
        if not ok:
            failures.append((item["tool"], doi, detail))
        time.sleep(0.2)
    print(f"\n{checked} DOIs checked, {len(failures)} failed")
    for tool, doi, detail in failures:
        print(f"  {tool}: {doi} ({detail})")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
