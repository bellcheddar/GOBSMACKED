"""Central configuration.

Paths anchor to the repository root (the parent of this package) so the app runs
identically from a laptop checkout and from /opt/gobsmacked on the droplet.
Secrets and the deploy target come from a gitignored .env; everything else has a
default here.
"""

from __future__ import annotations

import os
from pathlib import Path

try:
    from dotenv import load_dotenv
except ImportError:  # dotenv is optional (bare python3 running `init`)
    def load_dotenv(*_args, **_kwargs):  # type: ignore
        return False

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PACKAGE_DIR = Path(__file__).resolve().parent
ROOT_DIR = PACKAGE_DIR.parent
DATA_DIR = Path(os.environ.get("GOBSMACKED_DATA", ROOT_DIR / "data"))
RUNS_DIR = DATA_DIR / "runs"                 # data/runs/<job_id>/
STRUCT_CACHE = DATA_DIR / "structures"       # fetched PDB/AFDB/ESM models
DB_PATH = Path(os.environ.get("GOBSMACKED_DB", ROOT_DIR / "gobsmacked.db"))
BUNDLE_TEMPLATE_DIR = ROOT_DIR / "bundle_template"
STATIC_DIR = PACKAGE_DIR / "static"
VENDOR_DIR = STATIC_DIR / "vendor"

load_dotenv(ROOT_DIR / ".env")

# ---------------------------------------------------------------------------
# External API endpoints
# ---------------------------------------------------------------------------
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
AFDB_PREDICTION_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
ESM_ATLAS_FOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
RCSB_FILE_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
INTERPRO_URL = (
    "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{accession}/"
    "?page_size=100"
)
KLIFS_BASE = "https://klifs.net/api"
GPCRDB_BASE = "https://gpcrdb.org/services"

HTTP_TIMEOUT = 60
HTTP_MAX_RETRIES = 5
USER_AGENT = "GOBSMACKED/1.0 (+https://gobsmacked.mdeller.com)"

# Annotation and reference caches are keyed `source:identifier` and expire after
# this many days. KLIFS and GPCRdb change on a release cadence measured in
# months, so a month of staleness costs nothing and saves every page load a
# round trip to three services.
CACHE_TTL_DAYS = 30

# ---------------------------------------------------------------------------
# Pocket geometry
# ---------------------------------------------------------------------------
# Box side = residue extent + this padding, floored at BOX_MIN_SIDE. Both in A.
BOX_PADDING = 8.0
BOX_MIN_SIDE = 18.0
# Residues within this distance of the reference ligand define "the pocket" for
# superposition. 8 A is the usual binding-site definition and keeps fold errors
# far from the site out of the ligand RMSD.
POCKET_RADIUS = 8.0

# ---------------------------------------------------------------------------
# Uploads and retention
# ---------------------------------------------------------------------------
MAX_UPLOAD_BYTES = int(os.environ.get("MAX_UPLOAD_BYTES", 512 * 1024 * 1024))
# Nightly prune drops results.tar.gz for public runs older than this many days.
# The scorecard, report and final PDBs stay, so the Runs row stays useful.
ARCHIVE_RETENTION_DAYS = int(os.environ.get("ARCHIVE_RETENTION_DAYS", "90"))

# ---------------------------------------------------------------------------
# Deploy / serving (from .env)
# ---------------------------------------------------------------------------
DROPLET_SSH = os.environ.get("DROPLET_SSH", "")
DROPLET_PATH = os.environ.get("DROPLET_PATH", "/opt/gobsmacked")
SERVER_NAME = os.environ.get("SERVER_NAME", "gobsmacked.mdeller.com")
BIND_ADDR = os.environ.get("BIND_ADDR", "127.0.0.1:8009")

VERSION = "1.0.0"


def ensure_dirs() -> None:
    """Create the runtime directory tree if missing. Safe to call repeatedly."""
    for d in (DATA_DIR, RUNS_DIR, STRUCT_CACHE):
        d.mkdir(parents=True, exist_ok=True)
