#!/usr/bin/env python3
"""Nightly prune: drop old results archives, keep everything that stays useful.

A public run older than the retention window loses its `results.tar.gz` and its
extracted trajectory, which are the large files. The scorecard, the report, the
final PDBs and the Runs row all stay, so the run remains readable: what goes is
the ability to download the raw archive and redraw the trajectory panels.

Private runs are never pruned. Their owner is the only person who can see them,
and deleting their data on a schedule they never agreed to would be wrong.
"""

from __future__ import annotations

import shutil
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from app import config, db  # noqa: E402

# Kept when the archive goes: everything the results page reads directly.
KEEP = {"complex_md_final.pdb", "complex_min.pdb", "complex_pose1.pdb", "model_apo.pdb",
        "reference.pdb", "manifest.json", "campaign.yaml", "pandamap_md_final.png",
        "pandamap_reference.png"}


def main() -> int:
    cutoff = datetime.now(timezone.utc) - timedelta(days=config.ARCHIVE_RETENTION_DAYS)
    freed = 0
    pruned = 0
    with db.cursor() as conn:
        rows = conn.execute(
            "SELECT job_id, created, visibility FROM jobs "
            "WHERE visibility = 'public' AND status = 'analysed'"
        ).fetchall()

    for row in rows:
        try:
            created = datetime.fromisoformat(row["created"])
        except ValueError:
            continue
        if created.tzinfo is None:
            created = created.replace(tzinfo=timezone.utc)
        if created > cutoff:
            continue
        run_dir = config.RUNS_DIR / row["job_id"]
        archive = run_dir / "results.tar.gz"
        if archive.exists():
            freed += archive.stat().st_size
            archive.unlink()
            pruned += 1
        traj = run_dir / "results" / "traj"
        if traj.exists():
            freed += sum(p.stat().st_size for p in traj.rglob("*") if p.is_file())
            shutil.rmtree(traj, ignore_errors=True)
        results = run_dir / "results"
        if results.exists():
            for path in results.iterdir():
                if path.is_file() and path.name not in KEEP:
                    freed += path.stat().st_size
                    path.unlink()
        db.update_job(row["job_id"], results_path=None)

    print(f"pruned {pruned} runs older than {config.ARCHIVE_RETENTION_DAYS} days, "
          f"freed {freed / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
