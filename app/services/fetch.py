"""Input resolution: what did the user paste, and what structure do we start from?

Panel 1 of Prepare accepts a UniProt accession, a raw sequence, a PDB ID or an
uploaded coordinate file, and has to answer two questions: what is the canonical
sequence, and which structure does the bundle start from. The structure priority
is fixed and reported to the user rather than being silently applied:

    user upload  >  named PDB entry  >  AlphaFold DB  >  ESM Atlas  >  fold in the bundle

"Fold in the bundle" is the fallback, not a failure: ESMFold runs on the user's
GPU box, which is the whole point of the two-stage design.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import gemmi

from .. import config, db
from . import http

# A UniProt accession: the two shapes the ID scheme actually allows.
UNIPROT_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")
PDB_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
AA_LETTERS = set("ACDEFGHIKLMNPQRSTVWYXBZUO")

# The ESM Atlas fold endpoint refuses long sequences and is slow well before its
# own limit. Past this we go straight to folding in the bundle, where there is a
# GPU and no request timeout.
ESM_ATLAS_MAX_LEN = 400


@dataclass
class ProteinInput:
    """Everything Panel 1 resolved, ready for the campaign and the UI."""
    kind: str                            # uniprot | sequence | pdb | upload
    accession: Optional[str] = None
    protein_name: str = ""
    gene: str = ""
    organism: str = ""
    sequence: str = ""
    features: list[dict] = field(default_factory=list)
    source_structure: str = "fold"       # afdb | esm_atlas | pdb | user_pdb | fold
    source_id: Optional[str] = None
    chain: str = "A"
    structure_path: Optional[str] = None  # server-side path to the fetched model
    mean_plddt: Optional[float] = None
    resolution: Optional[float] = None
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        d = dict(self.__dict__)
        # The absolute path is a server detail: the browser gets a URL instead.
        d.pop("structure_path", None)
        d["has_structure"] = self.structure_path is not None
        return d


def classify_input(text: str) -> str:
    """uniprot | pdb | sequence | unknown, from what the user typed."""
    t = text.strip().replace(" ", "").replace("\n", "")
    if not t:
        return "unknown"
    upper = t.upper()
    if UNIPROT_RE.match(upper):
        return "uniprot"
    if PDB_RE.match(upper):
        return "pdb"
    letters = set(upper)
    if len(upper) >= 20 and letters <= AA_LETTERS:
        return "sequence"
    return "unknown"


def clean_sequence(text: str) -> str:
    """FASTA or free text to a bare one-letter sequence."""
    lines = [ln for ln in text.splitlines() if not ln.startswith(">")]
    return re.sub(r"[^A-Za-z]", "", "".join(lines)).upper()


# ---------------------------------------------------------------------------
# UniProt
# ---------------------------------------------------------------------------

# Feature types worth drawing on the sequence track. UniProt emits dozens;
# these are the ones that say something about a binding site.
TRACK_FEATURES = {
    "Binding site", "Active site", "Site", "Metal binding", "Transmembrane",
    "Domain", "Region", "Motif", "Modified residue", "Disulfide bond",
    "Mutagenesis", "Glycosylation",
}


def fetch_uniprot(accession: str) -> Optional[dict]:
    """Canonical sequence, name and features for an accession. None if unknown."""
    accession = accession.strip().upper()
    cached = db.cache_get("uniprot", accession)
    if cached is not None:
        return cached
    data = http.get_json(config.UNIPROT_ENTRY_URL.format(accession=accession))
    if not data:
        return None
    seq = (data.get("sequence") or {}).get("value", "")
    desc = data.get("proteinDescription") or {}
    rec = desc.get("recommendedName") or (desc.get("submissionNames") or [{}])[0]
    name = ((rec or {}).get("fullName") or {}).get("value", "")
    genes = data.get("genes") or []
    gene = ((genes[0].get("geneName") or {}).get("value", "")) if genes else ""
    organism = (data.get("organism") or {}).get("scientificName", "")

    features = []
    for f in data.get("features", []):
        ftype = f.get("type", "")
        if ftype not in TRACK_FEATURES:
            continue
        loc = f.get("location") or {}
        start = ((loc.get("start") or {}).get("value"))
        end = ((loc.get("end") or {}).get("value"))
        if start is None or end is None:
            continue
        features.append({
            "type": ftype,
            "start": int(start),
            "end": int(end),
            "description": f.get("description", ""),
        })

    payload = {
        "accession": accession,
        "name": name,
        "gene": gene,
        "organism": organism,
        "sequence": seq,
        "length": len(seq),
        "features": features,
    }
    db.cache_put("uniprot", accession, payload)
    return payload


def uniprot_for_sequence(sequence: str) -> Optional[str]:
    """Best-effort accession for a pasted sequence, by exact-length UniProt search.

    Only an exact sequence match counts. A near match would put the wrong
    annotation track under the user's residues, which is worse than none.
    """
    if len(sequence) < 30:
        return None
    try:
        data = http.get_json(
            config.UNIPROT_SEARCH_URL,
            params={"query": f"length:[{len(sequence)} TO {len(sequence)}]",
                    "fields": "accession,sequence", "size": 25, "format": "json"},
        )
    except Exception:
        return None
    for item in (data or {}).get("results", []):
        if (item.get("sequence") or {}).get("value", "") == sequence:
            return item.get("primaryAccession")
    return None


# ---------------------------------------------------------------------------
# Structures
# ---------------------------------------------------------------------------

def fetch_afdb(accession: str) -> Optional[dict]:
    """The AlphaFold DB model for an accession, downloaded to the cache."""
    data = http.get_json(config.AFDB_PREDICTION_URL.format(accession=accession))
    if not data:
        return None
    # Fragment F1 covers everything under 2700 residues, which is every target
    # this app is realistically pointed at; take the first if several.
    entry = data[0]
    url = entry.get("pdbUrl")
    if not url:
        return None
    entry_id = entry.get("entryId", f"AF-{accession}-F1")
    dest = config.STRUCT_CACHE / f"{entry_id}.pdb"
    path = http.download(url, dest)
    if path is None:
        return None
    return {
        "source_structure": "afdb",
        "source_id": entry_id,
        "path": str(path),
        "mean_plddt": _as_float(entry.get("globalMetricValue")),
        "sequence": entry.get("uniprotSequence") or entry.get("sequence") or "",
    }


def fetch_esm_atlas(sequence: str, job_hint: str = "esm") -> Optional[dict]:
    """Fold a short sequence with the ESM Atlas API.

    The public endpoint returns a PDB string with pLDDT in the B-factor column.
    Anything it refuses (length, load, downtime) falls through to folding in the
    bundle, so a failure here is not an error.
    """
    if len(sequence) > ESM_ATLAS_MAX_LEN:
        return None
    try:
        resp = http.post(config.ESM_ATLAS_FOLD_URL, data=sequence.encode("utf-8"))
        if resp.status_code != 200 or "ATOM" not in resp.text:
            return None
    except Exception:
        return None
    dest = config.STRUCT_CACHE / f"esm_{job_hint}.pdb"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(resp.text)
    return {
        "source_structure": "esm_atlas",
        "source_id": f"ESMAtlas:{job_hint}",
        "path": str(dest),
        "mean_plddt": mean_bfactor(dest),
    }


def fetch_pdb(pdb_id: str) -> Optional[dict]:
    """Download an RCSB entry (mmCIF) and read its resolution."""
    pdb_id = pdb_id.strip().upper()
    if not PDB_RE.match(pdb_id):
        return None
    dest = config.STRUCT_CACHE / f"{pdb_id}.cif"
    path = http.download(config.RCSB_FILE_URL.format(pdb_id=pdb_id), dest)
    if path is None:
        return None
    resolution = None
    try:
        meta = http.get_json(config.RCSB_ENTRY_URL.format(pdb_id=pdb_id))
        res = ((meta or {}).get("rcsb_entry_info") or {}).get("resolution_combined")
        if res:
            resolution = float(res[0])
    except Exception:
        pass
    return {
        "source_structure": "pdb",
        "source_id": pdb_id,
        "path": str(path),
        "resolution": resolution,
    }


def mean_bfactor(path: Path) -> Optional[float]:
    """Mean CA B-factor, which for a predicted model is the mean pLDDT."""
    try:
        st = gemmi.read_structure(str(path))
    except Exception:
        return None
    vals = [a.b_iso for model in st for ch in model for res in ch
            for a in res if a.name == "CA"]
    return round(sum(vals) / len(vals), 1) if vals else None


def structure_chains(path: str | Path) -> list[dict]:
    """Chain id, one-letter sequence and residue span for each polymer chain."""
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    out = []
    for chain in st[0]:
        poly = chain.get_polymer()
        if len(poly) < 10:
            continue
        seq = gemmi.one_letter_code(poly.extract_sequence())
        residues = [r.seqid.num for r in poly]
        out.append({
            "chain": chain.name,
            "sequence": seq.upper(),
            "length": len(seq),
            "first": min(residues) if residues else None,
            "last": max(residues) if residues else None,
        })
    return out


def structure_ligands(path: str | Path, min_atoms: int = 6) -> list[dict]:
    """Non-water heteroatom groups, largest first: the entry's candidate ligands."""
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_waters()
    seen: dict[str, dict] = {}
    for chain in st[0]:
        for res in chain:
            if res.het_flag != "H" or res.is_water():
                continue
            if len(res) < min_atoms:
                continue
            key = res.name
            entry = seen.setdefault(key, {"ccd": res.name, "copies": 0, "atoms": len(res),
                                          "chain": chain.name, "seqid": res.seqid.num})
            entry["copies"] += 1
    return sorted(seen.values(), key=lambda d: -d["atoms"])


def _as_float(v) -> Optional[float]:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# The priority ladder
# ---------------------------------------------------------------------------

def resolve_protein(text: str, pdb_id: str = "", uploaded: Optional[Path] = None,
                    job_hint: str = "job") -> ProteinInput:
    """Run the whole ladder and report which rung was used."""
    text = (text or "").strip()
    kind = "upload" if uploaded else classify_input(text)
    result = ProteinInput(kind=kind)

    # 1. Identity and sequence.
    accession = None
    if kind == "uniprot":
        accession = text.upper()
    elif kind == "sequence":
        result.sequence = clean_sequence(text)
        accession = uniprot_for_sequence(result.sequence)
        if accession:
            result.warnings.append(f"Sequence matches UniProt {accession} exactly; using its annotation.")

    if accession:
        entry = fetch_uniprot(accession)
        if entry:
            result.accession = accession
            result.protein_name = entry["name"]
            result.gene = entry["gene"]
            result.organism = entry["organism"]
            result.features = entry["features"]
            if not result.sequence:
                result.sequence = entry["sequence"]
        else:
            result.warnings.append(f"UniProt has no entry {accession}.")

    # 2. Structure, in priority order.
    if uploaded is not None:
        chains = structure_chains(uploaded)
        result.source_structure = "user_pdb"
        result.source_id = Path(uploaded).name
        result.structure_path = str(uploaded)
        result.mean_plddt = mean_bfactor(Path(uploaded))
        if chains:
            result.chain = chains[0]["chain"]
            if not result.sequence:
                result.sequence = chains[0]["sequence"]
        return result

    wanted_pdb = pdb_id.strip().upper() or (text.upper() if kind == "pdb" else "")
    if wanted_pdb:
        got = fetch_pdb(wanted_pdb)
        if got:
            result.source_structure = "pdb"
            result.source_id = got["source_id"]
            result.structure_path = got["path"]
            result.resolution = got.get("resolution")
            chains = structure_chains(got["path"])
            if chains:
                result.chain = chains[0]["chain"]
                if not result.sequence:
                    result.sequence = chains[0]["sequence"]
            return result
        result.warnings.append(f"RCSB has no entry {wanted_pdb}; falling back to a model.")

    if result.accession:
        got = fetch_afdb(result.accession)
        if got:
            result.source_structure = "afdb"
            result.source_id = got["source_id"]
            result.structure_path = got["path"]
            result.mean_plddt = got.get("mean_plddt")
            if not result.sequence:
                result.sequence = got.get("sequence", "")
            return result

    if result.sequence:
        got = fetch_esm_atlas(result.sequence, job_hint=job_hint)
        if got:
            result.source_structure = "esm_atlas"
            result.source_id = got["source_id"]
            result.structure_path = got["path"]
            result.mean_plddt = got.get("mean_plddt")
            return result
        result.source_structure = "fold"
        result.warnings.append(
            "No experimental or database model found: the bundle will fold this "
            "sequence with ESMFold on your GPU."
        )
        return result

    result.warnings.append("Nothing recognisable in that input: paste a UniProt accession, a PDB ID or a sequence.")
    return result
