#!/usr/bin/env python3
"""Regenerate THIRD_PARTY.md from software.yaml.

The About page reads the same file, so the page and the file cannot disagree.
Run `make third-party` after editing software.yaml.
"""

from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
HEADER = """# Third-party software

GOBSMACKED itself is MIT licensed. This file lists everything it uses, what it
uses it for, and under what licence.

**PLIP is GPL-2.0.** It is invoked as a subprocess and never imported, so no GPL
code is linked into this application, and it is not vendored into the run
bundle. Everything else here is permissively licensed or public data.

This file is generated from `software.yaml` by `make third-party`. Edit that file,
not this one.

"""


def main() -> int:
    data = yaml.safe_load((ROOT / "software.yaml").read_text()) or {}
    items = data.get("software", [])
    stages: list[str] = []
    for item in items:
        if item.get("stage") and item["stage"] not in stages:
            stages.append(item["stage"])

    lines = [HEADER]
    for stage in stages:
        lines.append(f"## {stage}\n")
        lines.append("| Tool | Role | Licence | Repository | Reference |")
        lines.append("|---|---|---|---|---|")
        for item in [i for i in items if i.get("stage") == stage]:
            doi = item.get("doi")
            reference = item.get("reference", "")
            if doi:
                reference = f"[{reference or doi}](https://doi.org/{doi})"
            repo = item.get("repo", "")
            repo_cell = f"[{repo.replace('https://', '')}]({repo})" if repo else ""
            lines.append(f"| {item['tool']} | {item['role']} | {item['licence']} | "
                         f"{repo_cell} | {reference} |")
        lines.append("")

    (ROOT / "THIRD_PARTY.md").write_text("\n".join(lines) + "\n")
    print(f"THIRD_PARTY.md: {len(items)} entries across {len(stages)} stages")
    return 0


if __name__ == "__main__":
    sys.exit(main())
