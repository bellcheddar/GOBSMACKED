"""Annotation and the family router.

InterPro supplies Pfam domains; the Pfam accessions decide which specialist
service is consulted:

    PF00069 / PF07714              -> kinase  -> KLIFS
    PF00001 / PF00002 / PF00003 / PF10324 -> gpcr -> GPCRdb
    anything else                  -> other   -> UniProt features only

Both specialist clients share the same cache pattern as everything else here, so
a repeat visit to a target costs no external calls at all.
"""

from __future__ import annotations

import re
from typing import Any, Optional

from .. import config, db
from . import http

KINASE_PFAM = {"PF00069", "PF07714"}
GPCR_PFAM = {"PF00001", "PF00002", "PF00003", "PF10324"}

# KLIFS positions this app names. The full 85-residue pocket is fetched and
# mapped, but these five carry the DFG / alphaC logic and the hinge.
KLIFS_BETA3_LYS = 17
KLIFS_ALPHAC_GLU = 24
KLIFS_GATEKEEPER = 45
KLIFS_HINGE = (46, 47, 48)
KLIFS_DFG_ASP = 81
KLIFS_DFG_PHE = 82
KLIFS_DFG_GLY = 83

# The eight KLIFS subpockets the Mode view reports occupancy for, as spans of
# the 85-position pocket numbering.
KLIFS_REGIONS = {
    "beta1": (1, 3),
    "g.loop": (4, 9),
    "beta2": (10, 13),
    "beta3": (14, 19),
    "alphaC": (20, 30),
    "back.loop": (31, 37),
    "beta4": (38, 41),
    "beta5": (42, 44),
    "gatekeeper": (45, 45),
    "hinge": (46, 48),
    "linker": (49, 52),
    "alphaD": (53, 59),
    "alphaE": (60, 64),
    "beta6": (65, 67),
    "cat.loop": (68, 75),
    "beta7": (76, 78),
    "beta8": (79, 80),
    "xDFG": (81, 84),
    "a.loop": (85, 85),
}

# GPCRdb generic numbers that the microswitch panel measures. Ballesteros-
# Weinstein, so they are family-independent by construction.
GPCR_SWITCHES = {
    "TM3_ionic_lock": "3.50",
    "TM6_ionic_lock": "6.30",
    "toggle": "6.48",
    "PIF_P": "5.50",
    "PIF_I": "3.40",
    "PIF_F": "6.44",
    "NPxxY_N": "7.49",
    "NPxxY_Y": "7.53",
    "sodium_D": "2.50",
}


# ---------------------------------------------------------------------------
# Pfam via InterPro
# ---------------------------------------------------------------------------

def pfam_domains(accession: str) -> list[dict]:
    """Pfam entries hit by this accession, with their residue spans."""
    cached = db.cache_get("pfam", accession)
    if cached is not None:
        return cached
    try:
        data = http.get_json(config.INTERPRO_URL.format(accession=accession))
    except Exception:
        return []
    out: list[dict] = []
    for item in (data or {}).get("results", []):
        meta = item.get("metadata", {})
        acc = meta.get("accession", "")
        for prot in item.get("proteins", []):
            for loc in prot.get("entry_protein_locations", []) or []:
                frags = loc.get("fragments") or []
                if not frags:
                    continue
                out.append({
                    "pfam": acc,
                    "name": meta.get("name", ""),
                    "start": int(frags[0]["start"]),
                    "end": int(frags[-1]["end"]),
                })
    out.sort(key=lambda d: d["start"])
    db.cache_put("pfam", accession, out)
    return out


def family_from_pfam(domains: list[dict]) -> str:
    accs = {d["pfam"] for d in domains}
    if accs & KINASE_PFAM:
        return "kinase"
    if accs & GPCR_PFAM:
        return "gpcr"
    return "other"


# ---------------------------------------------------------------------------
# KLIFS
# ---------------------------------------------------------------------------

def klifs_annotation(accession: str, sequence: str, gene: str = "") -> Optional[dict]:
    """KLIFS kinase identity, the 85-residue pocket mapped onto `sequence`, and
    the DFG / alphaC state of every KLIFS structure for the kinase.

    The pocket is returned by KLIFS as an 85-character string with '-' for
    positions absent in that kinase. It carries no residue numbers, so the
    mapping is done here by locating each ungapped run in the target sequence.
    A run that is absent or ambiguous is left unmapped rather than guessed: a
    wrong number here would put the gatekeeper on the wrong residue and every
    downstream label with it.
    """
    cached = db.cache_get("klifs", accession)
    if cached is None:
        cached = _klifs_fetch(accession, gene)
        if cached is None:
            return None
        db.cache_put("klifs", accession, cached)
    if not cached:
        return None
    out = dict(cached)
    out["pocket_map"] = map_pocket_to_sequence(cached.get("pocket", ""), sequence)
    out["named_positions"] = {
        "beta3_lysine": out["pocket_map"].get(str(KLIFS_BETA3_LYS)),
        "alphaC_glutamate": out["pocket_map"].get(str(KLIFS_ALPHAC_GLU)),
        "gatekeeper": out["pocket_map"].get(str(KLIFS_GATEKEEPER)),
        "hinge": [out["pocket_map"].get(str(p)) for p in KLIFS_HINGE],
        "dfg_asp": out["pocket_map"].get(str(KLIFS_DFG_ASP)),
        "dfg_phe": out["pocket_map"].get(str(KLIFS_DFG_PHE)),
        "dfg_gly": out["pocket_map"].get(str(KLIFS_DFG_GLY)),
    }
    return out


def _klifs_fetch(accession: str, gene: str) -> Optional[dict]:
    """Kinase identity and structure list from klifs.net. {} means 'not a kinase there'."""
    kinase_id = None
    info: dict[str, Any] = {}
    for name in filter(None, [gene, accession]):
        try:
            hits = http.get_json(f"{config.KLIFS_BASE}/kinase_ID",
                                 params={"kinase_name": name, "species": "HUMAN"})
        except Exception:
            hits = None
        if not hits:
            continue
        # Prefer the hit whose UniProt accession matches; KLIFS name search is
        # fuzzy and will happily return a paralogue for a gene symbol.
        chosen = next((h for h in hits if h.get("uniprot") == accession), hits[0])
        kinase_id = chosen.get("kinase_ID")
        info = chosen
        break
    if kinase_id is None:
        return {}

    structures = []
    try:
        rows = http.get_json(f"{config.KLIFS_BASE}/structures_list",
                             params={"kinase_ID": kinase_id}) or []
        for s in rows:
            structures.append({
                "pdb": (s.get("pdb") or "").upper(),
                "chain": s.get("chain", ""),
                "alt": s.get("alt", ""),
                "ligand": s.get("ligand", ""),
                "dfg": s.get("DFG", ""),
                "ac_helix": s.get("ac_helix", ""),
                "resolution": s.get("resolution"),
                "pocket": s.get("pocket", ""),
            })
    except Exception:
        pass

    return {
        "kinase_id": kinase_id,
        "name": info.get("name", ""),
        "full_name": info.get("full_name", ""),
        "family": info.get("family", ""),
        "group": info.get("group", ""),
        "uniprot": info.get("uniprot", ""),
        "pocket": info.get("pocket", "") or (structures[0]["pocket"] if structures else ""),
        "structures": structures,
    }


def map_pocket_to_sequence(pocket: str, sequence: str) -> dict[str, int]:
    """KLIFS position (1-85, as a string key) -> residue number in `sequence`.

    The 85-character pocket string is a concatenation of the KLIFS regions, and
    only the regions are contiguous in the sequence: EGFR's beta2 ends at Lys725
    and its beta3 begins at Val742, with sixteen residues in between that the
    pocket does not contain. So the mapping is done region by region, longest
    run first: try to place regions i..j as one substring, shorten j until a
    match is found, then continue from where that match ended. Matching left to
    right with a moving cursor keeps the result monotonic.

    A '-' in the pocket marks a position this kinase does not have. Its
    neighbours are still adjacent in the sequence, so gaps are stripped before
    searching and skipped when the residue numbers are handed out.

    A region that cannot be placed is left out rather than guessed. A wrong
    number here would put the gatekeeper on the wrong residue, and every label
    downstream with it.
    """
    if not pocket or not sequence:
        return {}
    sequence = sequence.upper()
    pocket = pocket.upper()
    bounds = sorted(KLIFS_REGIONS.values())          # [(1,3), (4,9), ...]
    mapping: dict[str, int] = {}
    cursor = 0
    i = 0
    while i < len(bounds):
        placed = False
        for j in range(len(bounds), i, -1):
            lo = bounds[i][0]
            hi = bounds[j - 1][1]
            raw = pocket[lo - 1:hi]
            seg = raw.replace("-", "")
            if len(seg) < 3:
                continue
            at = sequence.find(seg, cursor)
            if at < 0:
                continue
            offset = 0
            for k, ch in enumerate(raw):
                if ch == "-":
                    continue
                mapping[str(lo + k)] = at + offset + 1     # both 1-based
                offset += 1
            cursor = at + len(seg)
            i = j
            placed = True
            break
        if not placed:
            i += 1
    _fill_gapped_regions(mapping, pocket, sequence, bounds)
    return mapping


def _fill_gapped_regions(mapping: dict[str, int], pocket: str, sequence: str,
                         bounds: list[tuple[int, int]]) -> None:
    """Place regions the exact search could not, inside the window their
    neighbours leave open.

    KLIFS's pocket is an alignment, so a region can differ from the sequence by
    an insertion: EGFR's back loop reads VDPHVCR in the pocket and VDNPHVCR in
    UniProt, one residue apart, and an exact substring search finds nothing.
    Here the search is bounded on both sides by already-placed positions, and
    only matching blocks of three residues or more are used, so a partial
    region contributes what it can and nothing is invented.
    """
    import difflib

    for lo, hi in bounds:
        if any(str(p) in mapping for p in range(lo, hi + 1)):
            continue
        raw = pocket[lo - 1:hi]
        seg = raw.replace("-", "")
        if len(seg) < 3:
            continue
        before = [mapping[str(p)] for p in range(1, lo) if str(p) in mapping]
        after = [mapping[str(p)] for p in range(hi + 1, 86) if str(p) in mapping]
        start = max(before) if before else 0          # 1-based residue number
        end = min(after) - 1 if after else len(sequence)
        if end - start < len(seg):
            continue
        window = sequence[start:end]                  # sequence[start] is residue start+1
        matcher = difflib.SequenceMatcher(None, seg, window, autojunk=False)
        for a, b, size in matcher.get_matching_blocks():
            if size < 3:
                continue
            # `a` indexes the ungapped segment; walk `raw` to recover the
            # KLIFS position each matched residue belongs to.
            ungapped = 0
            for k, ch in enumerate(raw):
                if ch == "-":
                    continue
                if a <= ungapped < a + size:
                    mapping[str(lo + k)] = start + b + (ungapped - a) + 1
                ungapped += 1


def klifs_region_of(position: int) -> str:
    for name, (lo, hi) in KLIFS_REGIONS.items():
        if lo <= position <= hi:
            return name
    return ""


# ---------------------------------------------------------------------------
# GPCRdb
# ---------------------------------------------------------------------------

def gpcrdb_annotation(accession: str) -> Optional[dict]:
    """Generic (Ballesteros-Weinstein) numbering and segment assignment.

    GPCRdb keys everything on its own entry name (`adrb2_human`), so the
    accession is resolved first. The residue endpoint returns one row per
    residue with `display_generic_number`, which is the number every GPCR paper
    quotes and the only stable way to compare two receptors' pockets.
    """
    cached = db.cache_get("gpcrdb", accession)
    if cached is not None:
        return cached or None

    try:
        prot = http.get_json(f"{config.GPCRDB_BASE}/protein/accession/{accession}/")
    except Exception:
        prot = None
    if not prot or not prot.get("entry_name"):
        db.cache_put("gpcrdb", accession, {})
        return None
    entry_name = prot["entry_name"]

    try:
        residues = http.get_json(f"{config.GPCRDB_BASE}/residues/extended/{entry_name}/") or []
    except Exception:
        residues = []

    generic: dict[str, int] = {}     # "3.50" -> residue number
    segments: list[dict] = []
    current: Optional[dict] = None
    for r in residues:
        num = r.get("sequence_number")
        disp = (r.get("display_generic_number") or "").split("x")[0]
        seg = r.get("protein_segment") or ""
        if disp and num:
            generic.setdefault(disp, int(num))
        if seg:
            if current and current["segment"] == seg and num == current["end"] + 1:
                current["end"] = int(num)
            else:
                current = {"segment": seg, "start": int(num), "end": int(num)}
                segments.append(current)

    payload = {
        "entry_name": entry_name,
        "receptor_family": (prot.get("family") or {}).get("name", "") if isinstance(prot.get("family"), dict) else prot.get("family", ""),
        "receptor_class": prot.get("receptor_class", ""),
        "species": prot.get("species", ""),
        "generic": generic,
        "segments": segments,
        "switches": {name: generic.get(gn) for name, gn in GPCR_SWITCHES.items()},
    }
    db.cache_put("gpcrdb", accession, payload)
    return payload


# ---------------------------------------------------------------------------
# The router
# ---------------------------------------------------------------------------

def annotate(accession: Optional[str], sequence: str, gene: str = "",
             features: Optional[list[dict]] = None) -> dict:
    """Everything Panel 2 renders: family, domains, and the family-specific track."""
    domains = pfam_domains(accession) if accession else []
    family = family_from_pfam(domains)
    out: dict[str, Any] = {
        "family": family,
        "pfam": domains,
        "features": features or [],
        "klifs": None,
        "gpcrdb": None,
        "positions": [],       # phos ticks on the sequence track
        "notes": [],
    }

    if family == "kinase" and accession:
        klifs = klifs_annotation(accession, sequence, gene)
        if klifs:
            out["klifs"] = klifs
            named = klifs["named_positions"]
            for label, value in (("beta3 Lys", named["beta3_lysine"]),
                                 ("alphaC Glu", named["alphaC_glutamate"]),
                                 ("gatekeeper", named["gatekeeper"]),
                                 ("DFG Asp", named["dfg_asp"]),
                                 ("DFG Phe", named["dfg_phe"])):
                if value:
                    out["positions"].append({"residue": value, "label": label})
            for i, h in enumerate(named["hinge"]):
                if h:
                    out["positions"].append({"residue": h, "label": f"hinge {i + 1}"})
            mapped = len(klifs.get("pocket_map", {}))
            out["notes"].append(f"KLIFS pocket: {mapped} of 85 positions mapped onto this sequence.")
            if mapped < 60:
                out["notes"].append(
                    "Fewer than 60 KLIFS positions mapped: the sequence may be a "
                    "construct or an isoform. Check the DFG and gatekeeper residues before docking."
                )
        else:
            out["notes"].append("Pfam says kinase, but KLIFS has no entry for this accession.")

    elif family == "gpcr" and accession:
        gp = gpcrdb_annotation(accession)
        if gp:
            out["gpcrdb"] = gp
            for name, num in gp["switches"].items():
                if num:
                    out["positions"].append({"residue": num, "label": name.replace("_", " ")})
            out["notes"].append(
                f"GPCRdb {gp['entry_name']}: {len(gp['generic'])} residues carry generic numbering."
            )
        else:
            out["notes"].append("Pfam says GPCR, but GPCRdb has no entry for this accession.")

    else:
        binding = [f for f in (features or []) if f["type"] in ("Binding site", "Active site", "Site")]
        for f in binding:
            out["positions"].append({"residue": f["start"], "label": f["type"].lower()})
        if not binding:
            out["notes"].append("No Pfam family this app specialises in, and UniProt lists no binding sites: pick the pocket by hand.")

    return out
