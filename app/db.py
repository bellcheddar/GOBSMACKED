"""SQLite schema and helpers.

One file, WAL mode, no ORM. Three tables: `jobs` (one row per Prepare, filled in
by Analyze), `annotation_cache` and `reference_cache` (both keyed by source and
identifier, both with a TTL).

Ownership: a run carries an `owner_token`, a 32-character secret issued at
Prepare time and stored only as its sha256. Possession of the token proves
ownership; the server never holds anything that could reconstruct it.
"""

from __future__ import annotations

import hashlib
import json
import secrets
import sqlite3
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from typing import Any, Iterable, Optional

from . import config

SCHEMA = """
CREATE TABLE IF NOT EXISTS jobs (
    job_id          TEXT PRIMARY KEY,
    created         TEXT NOT NULL,
    updated         TEXT NOT NULL,
    title           TEXT,
    uniprot         TEXT,
    protein_name    TEXT,
    ligand_name     TEXT,
    ligand_smiles   TEXT,
    family          TEXT,
    reference_pdb   TEXT,
    status          TEXT NOT NULL,      -- prepared | results_uploaded | analysed | failed
    visibility      TEXT NOT NULL,      -- public | private
    owner_hash      TEXT NOT NULL,      -- sha256 of the owner token
    campaign_yaml   TEXT,
    results_path    TEXT,
    scorecard_json  TEXT,
    gobsmack_score  REAL,
    grade           TEXT,
    mode_predicted  TEXT,
    mode_reference  TEXT,
    mode_match      INTEGER,            -- 1 match, 0 differ, NULL unverified
    error           TEXT
);
CREATE INDEX IF NOT EXISTS jobs_created ON jobs(created DESC);
CREATE INDEX IF NOT EXISTS jobs_visibility ON jobs(visibility, created DESC);

CREATE TABLE IF NOT EXISTS annotation_cache (
    key       TEXT PRIMARY KEY,         -- "<source>:<identifier>"
    source    TEXT NOT NULL,
    fetched   REAL NOT NULL,            -- unix seconds
    payload   TEXT NOT NULL             -- JSON
);

CREATE TABLE IF NOT EXISTS reference_cache (
    uniprot   TEXT PRIMARY KEY,
    fetched   REAL NOT NULL,
    payload   TEXT NOT NULL
);
"""

# Columns the Runs table and the run header read. Kept as one list so the table
# view, the JSON API and the row-to-dict helper cannot drift apart.
RUN_COLUMNS = (
    "job_id", "created", "updated", "title", "uniprot", "protein_name",
    "ligand_name", "ligand_smiles", "family", "reference_pdb", "status",
    "visibility", "gobsmack_score", "grade", "mode_predicted",
    "mode_reference", "mode_match", "error",
)


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def connect() -> sqlite3.Connection:
    config.DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(config.DB_PATH, timeout=30)
    conn.row_factory = sqlite3.Row
    # WAL lets the nightly prune write while a request reads.
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    return conn


@contextmanager
def cursor():
    conn = connect()
    try:
        yield conn
        conn.commit()
    finally:
        conn.close()


def init_db() -> None:
    with cursor() as conn:
        conn.executescript(SCHEMA)


# ---------------------------------------------------------------------------
# Ownership
# ---------------------------------------------------------------------------

def new_owner_token() -> str:
    """A 32-character owner key. Shown once, stored only as a hash."""
    return secrets.token_hex(16)


def hash_token(token: str) -> str:
    return hashlib.sha256(token.strip().encode("utf-8")).hexdigest()


def token_matches(row: sqlite3.Row | dict, token: Optional[str]) -> bool:
    """Constant-time comparison of a supplied token against a row's hash."""
    if not token:
        return False
    stored = row["owner_hash"] if isinstance(row, sqlite3.Row) else row.get("owner_hash")
    if not stored:
        return False
    return secrets.compare_digest(hash_token(token), stored)


def new_job_id() -> str:
    """`gs_<yyyymmdd>_<12 base32 chars>`.

    The date is for the human reading a directory listing; the 12 random
    characters (60 bits) are what makes a private run's URL unguessable, which
    the Runs tab relies on.
    """
    alphabet = "abcdefghijklmnopqrstuvwxyz234567"
    tail = "".join(secrets.choice(alphabet) for _ in range(12))
    return f"gs_{datetime.now(timezone.utc):%Y%m%d}_{tail}"


# ---------------------------------------------------------------------------
# Jobs
# ---------------------------------------------------------------------------

def insert_job(**fields: Any) -> None:
    fields.setdefault("created", now_iso())
    fields.setdefault("updated", now_iso())
    cols = ", ".join(fields)
    marks = ", ".join("?" for _ in fields)
    with cursor() as conn:
        conn.execute(f"INSERT INTO jobs ({cols}) VALUES ({marks})", tuple(fields.values()))


def update_job(job_id: str, **fields: Any) -> None:
    if not fields:
        return
    fields["updated"] = now_iso()
    sets = ", ".join(f"{k} = ?" for k in fields)
    with cursor() as conn:
        conn.execute(f"UPDATE jobs SET {sets} WHERE job_id = ?", (*fields.values(), job_id))


def get_job(job_id: str) -> Optional[sqlite3.Row]:
    with cursor() as conn:
        return conn.execute("SELECT * FROM jobs WHERE job_id = ?", (job_id,)).fetchone()


def delete_job(job_id: str) -> None:
    with cursor() as conn:
        conn.execute("DELETE FROM jobs WHERE job_id = ?", (job_id,))


def list_jobs(owner_token: Optional[str] = None, limit: int = 500) -> list[dict]:
    """Public runs, plus private runs whose owner hash matches `owner_token`.

    Private runs the caller does not own are absent rather than redacted: a
    greyed-out row would leak that the run exists, and the count of them.
    """
    sql = f"SELECT {', '.join(RUN_COLUMNS)}, owner_hash FROM jobs"
    params: list[Any] = []
    if owner_token:
        sql += " WHERE visibility = 'public' OR owner_hash = ?"
        params.append(hash_token(owner_token))
    else:
        sql += " WHERE visibility = 'public'"
    sql += " ORDER BY created DESC LIMIT ?"
    params.append(limit)
    with cursor() as conn:
        rows = conn.execute(sql, params).fetchall()
    out = []
    for r in rows:
        d = {k: r[k] for k in RUN_COLUMNS}
        d["owned"] = bool(owner_token) and r["owner_hash"] == hash_token(owner_token)
        out.append(d)
    return out


def status_counts() -> dict[str, int]:
    """Aggregate counts for the Runs page stage strip (public runs only)."""
    with cursor() as conn:
        rows = conn.execute(
            "SELECT status, COUNT(*) AS n FROM jobs WHERE visibility = 'public' GROUP BY status"
        ).fetchall()
    counts = {"prepared": 0, "results_uploaded": 0, "analysed": 0, "failed": 0}
    for r in rows:
        counts[r["status"]] = r["n"]
    return counts


# ---------------------------------------------------------------------------
# Caches
# ---------------------------------------------------------------------------

def cache_get(source: str, identifier: str, ttl_days: int | None = None) -> Optional[Any]:
    ttl = (ttl_days if ttl_days is not None else config.CACHE_TTL_DAYS) * 86400
    key = f"{source}:{identifier}"
    with cursor() as conn:
        row = conn.execute(
            "SELECT fetched, payload FROM annotation_cache WHERE key = ?", (key,)
        ).fetchone()
    if not row or (time.time() - row["fetched"]) > ttl:
        return None
    try:
        return json.loads(row["payload"])
    except json.JSONDecodeError:
        return None


def cache_put(source: str, identifier: str, payload: Any) -> None:
    key = f"{source}:{identifier}"
    with cursor() as conn:
        conn.execute(
            "INSERT INTO annotation_cache (key, source, fetched, payload) VALUES (?,?,?,?) "
            "ON CONFLICT(key) DO UPDATE SET fetched = excluded.fetched, payload = excluded.payload",
            (key, source, time.time(), json.dumps(payload)),
        )


def reference_cache_get(uniprot: str, ttl_days: int | None = None) -> Optional[Any]:
    ttl = (ttl_days if ttl_days is not None else config.CACHE_TTL_DAYS) * 86400
    with cursor() as conn:
        row = conn.execute(
            "SELECT fetched, payload FROM reference_cache WHERE uniprot = ?", (uniprot,)
        ).fetchone()
    if not row or (time.time() - row["fetched"]) > ttl:
        return None
    try:
        return json.loads(row["payload"])
    except json.JSONDecodeError:
        return None


def reference_cache_put(uniprot: str, payload: Any) -> None:
    with cursor() as conn:
        conn.execute(
            "INSERT INTO reference_cache (uniprot, fetched, payload) VALUES (?,?,?) "
            "ON CONFLICT(uniprot) DO UPDATE SET fetched = excluded.fetched, payload = excluded.payload",
            (uniprot, time.time(), json.dumps(payload)),
        )
