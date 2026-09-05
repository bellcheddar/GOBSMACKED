"""Reference structure search: which crystal should the prediction be judged against?

RCSB is asked for every entry whose polymer entity maps to the target's UniProt
accession. Each entry's bound ligands are pulled in one GraphQL round trip, and
each ligand is compared to the docked ligand by Morgan fingerprint Tanimoto. The
default choice is the most similar ligand at better than 2.5 A, because a
reference with a different chemotype is a weak test of a docking pose.

"None" is a first-class answer: verification is switched off, the scorecard shows
an unverified banner, and everything else still works.
"""

from __future__ import annotations

from typing import Any, Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import DataStructs
from rdkit.Chem import rdFingerprintGenerator

from .. import config, db
from . import http

# RDKit is loud about sanitisation on ligand SMILES pulled from the PDB chemical
# component dictionary, and those messages are not actionable here.
RDLogger.DisableLog("rdApp.*")

_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

# Components present in most crystals that are never the ligand of interest.
IGNORED_CCD = {
    "HOH", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "MPD", "TRS", "ACT",
    "CL", "NA", "MG", "CA", "ZN", "MN", "K", "IOD", "BR", "NO3", "FMT", "DMS",
    "EPE", "IMD", "CIT", "TLA", "BME", "NH4", "CO3", "AZI", "FLC", "1PE", "P6G",
}

_ENTRY_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    rcsb_entry_info { resolution_combined }
    polymer_entities {
      rcsb_polymer_entity_container_identifiers { auth_asym_ids }
      entity_poly { rcsb_sample_sequence_length }
    }
    nonpolymer_entities {
      rcsb_nonpolymer_entity_container_identifiers { auth_asym_ids }
      nonpolymer_comp {
        chem_comp { id name formula_weight }
        rcsb_chem_comp_descriptor { SMILES_stereo SMILES }
      }
    }
  }
}
"""


def search_by_uniprot(accession: str, limit: int = 100) -> list[str]:
    """PDB IDs whose polymer entity references this UniProt accession."""
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers."
                             "reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": accession.upper(),
            },
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit},
                            "results_content_type": ["experimental"]},
    }
    resp = http.post(config.RCSB_SEARCH_URL, json=query)
    # RCSB answers a search with no hits with 204 No Content, not an empty
    # result set, and .json() on that raises rather than returning nothing.
    if resp.status_code == 204:
        return []
    resp.raise_for_status()
    return [r["identifier"] for r in resp.json().get("result_set", [])]


def entry_details(pdb_ids: list[str]) -> list[dict]:
    """Title, resolution, chains and bound ligands for each entry, in batches."""
    out: list[dict] = []
    for i in range(0, len(pdb_ids), 40):
        batch = pdb_ids[i:i + 40]
        try:
            data = http.post_json(config.RCSB_GRAPHQL_URL,
                                  {"query": _ENTRY_QUERY, "variables": {"ids": batch}})
        except Exception:
            continue
        for entry in ((data or {}).get("data") or {}).get("entries") or []:
            if not entry:
                continue
            res = (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or []
            chains: list[str] = []
            for pe in entry.get("polymer_entities") or []:
                ids = (pe.get("rcsb_polymer_entity_container_identifiers") or {}).get("auth_asym_ids") or []
                chains.extend(ids)
            ligands = []
            for ne in entry.get("nonpolymer_entities") or []:
                comp = (ne.get("nonpolymer_comp") or {})
                chem = comp.get("chem_comp") or {}
                desc = comp.get("rcsb_chem_comp_descriptor") or {}
                ccd = chem.get("id", "")
                if not ccd or ccd in IGNORED_CCD:
                    continue
                ligands.append({
                    "ccd": ccd,
                    "name": chem.get("name", ""),
                    "mw": chem.get("formula_weight"),
                    "smiles": desc.get("SMILES_stereo") or desc.get("SMILES") or "",
                    "chains": (ne.get("rcsb_nonpolymer_entity_container_identifiers") or {}).get("auth_asym_ids") or [],
                })
            out.append({
                "pdb_id": entry.get("rcsb_id", ""),
                "title": (entry.get("struct") or {}).get("title", ""),
                "resolution": float(res[0]) if res else None,
                "chains": sorted(set(chains)),
                "ligands": ligands,
            })
    return out


def fingerprint(smiles: str):
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    return _MORGAN.GetFingerprint(mol) if mol else None


def tanimoto(fp_a, fp_b) -> Optional[float]:
    if fp_a is None or fp_b is None:
        return None
    return round(DataStructs.TanimotoSimilarity(fp_a, fp_b), 3)


def rank_references(accession: str, smiles: str, refresh: bool = False) -> list[dict]:
    """Candidate reference entries, most chemically similar ligand first.

    The RCSB half of this (search plus GraphQL) is cached per accession, since
    it does not depend on the ligand. The Tanimoto half is recomputed every time
    because it does.
    """
    entries = None if refresh else db.reference_cache_get(accession)
    if entries is None:
        ids = search_by_uniprot(accession)
        entries = entry_details(ids) if ids else []
        db.reference_cache_put(accession, entries)

    query_fp = fingerprint(smiles)
    ranked: list[dict] = []
    for entry in entries:
        best = None
        for lig in entry["ligands"]:
            sim = tanimoto(query_fp, fingerprint(lig["smiles"]))
            row = dict(lig, tanimoto=sim)
            if best is None or (sim or -1) > (best.get("tanimoto") or -1):
                best = row
        ranked.append({
            "pdb_id": entry["pdb_id"],
            "title": entry["title"],
            "resolution": entry["resolution"],
            "chains": entry["chains"],
            "n_ligands": len(entry["ligands"]),
            "best_ligand": best,
            "tanimoto": (best or {}).get("tanimoto"),
        })

    # Holo entries first, then by similarity, then by resolution. An apo entry
    # (no ligand at all) still belongs in the list: it is the optional apo
    # reference the overlay can draw.
    ranked.sort(key=lambda e: (
        e["best_ligand"] is None,
        -(e["tanimoto"] if e["tanimoto"] is not None else -1),
        e["resolution"] if e["resolution"] is not None else 99,
    ))
    return ranked


def default_reference(ranked: list[dict], max_resolution: float = 2.5) -> Optional[dict]:
    """The entry to preselect: same ligand if one exists, else best similarity
    at better than `max_resolution`.

    Resolution is the tie-break, not the filter, when the reference contains the
    very ligand being docked. EGFR and erlotinib make the case: 1M17 holds
    erlotinib itself at 2.6 A, while the best sub-2.5 A entry holds gefitinib at
    Tanimoto 0.41. Judging an erlotinib pose against a gefitinib crystal because
    the crystal is 0.9 A sharper would measure the wrong thing.
    """
    holo = [e for e in ranked if e["best_ligand"] and e["tanimoto"] is not None]
    if not holo:
        return None
    same = [e for e in holo if e["tanimoto"] >= 0.99]
    if same:
        return sorted(same, key=lambda e: e["resolution"] if e["resolution"] is not None else 99)[0]
    good = [e for e in holo if (e["resolution"] or 99) < max_resolution]
    return (good or holo)[0]
