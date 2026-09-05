"""Interaction fingerprints: PLIP for the contacts, PandaMap for the picture.

PLIP is GPL-2.0. It is run as a **subprocess** and never imported, so this
MIT-licensed app links against nothing of it: what crosses the boundary is an
XML file. That is also why it is not vendored into the run bundle.

PandaMap is MIT, so it is imported directly, and gives the 2D interaction
diagram and its own empirical binding-energy estimate (reported as PandaMap's
number, with PandaMap's caveat, not as a GOBSMACKED score).

The comparison that matters is between fingerprints: the set of
(interaction type, residue) pairs, translated into the reference structure's
numbering so a predicted complex and a crystal can be compared residue by
residue.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Optional

# PLIP's XML groups interactions by tag. These are the ones with a residue on
# the protein side, which is all a fingerprint needs.
PLIP_GROUPS = {
    "hydrophobic_interactions": "hydrophobic",
    "hydrogen_bonds": "hbond",
    "water_bridges": "water bridge",
    "salt_bridges": "salt bridge",
    "pi_stacks": "pi stacking",
    "pi_cation_interactions": "pi cation",
    "halogen_bonds": "halogen bond",
    "metal_complexes": "metal",
}

PLIP_TIMEOUT = 180


class PlipUnavailable(RuntimeError):
    """PLIP is not installed here. Analysis continues without a fingerprint."""


def run_plip(structure: str | Path, ligand_ccd: Optional[str] = None,
             timeout: int = PLIP_TIMEOUT) -> dict[str, Any]:
    """PLIP's interactions for one complex.

    PLIP 3 reads PDB, not mmCIF, and answers a .cif with "no valid file format
    provided", so callers convert first (`as_pdb` below).
    """
    structure = Path(structure)
    if not structure.exists():
        return {"error": f"{structure.name} is not in the results archive."}

    with tempfile.TemporaryDirectory() as tmp:
        cmd = [sys.executable, "-m", "plip.plipcmd", "-f", str(structure),
               "-x", "-o", tmp, "--name", "report"]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        except FileNotFoundError as exc:
            raise PlipUnavailable("PLIP is not installed on this server.") from exc
        except subprocess.TimeoutExpired:
            return {"error": f"PLIP timed out after {timeout} s on {structure.name}."}
        report = Path(tmp) / "report.xml"
        if not report.exists():
            tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-1:] or ["no output"]
            return {"error": f"PLIP produced no report for {structure.name}: {tail[0]}"}
        return parse_plip_xml(report.read_text(), ligand_ccd)


def parse_plip_xml(xml_text: str, ligand_ccd: Optional[str] = None) -> dict[str, Any]:
    """PLIP's XML into a flat list of interactions for one binding site."""
    root = ET.fromstring(xml_text)
    sites = []
    for site in root.iter("bindingsite"):
        ident = site.find("identifiers")
        hetid = (ident.findtext("hetid") or "").strip() if ident is not None else ""
        if ligand_ccd and hetid.upper() != ligand_ccd.upper():
            continue
        rows: list[dict] = []
        interactions = site.find("interactions")
        if interactions is None:
            continue
        for group in interactions:
            kind = PLIP_GROUPS.get(group.tag)
            if kind is None:
                continue
            for item in group:
                resnr = item.findtext("resnr")
                if not resnr:
                    continue
                rows.append({
                    "type": kind,
                    "resnr": int(resnr),
                    "restype": (item.findtext("restype") or "").strip(),
                    "chain": (item.findtext("reschain") or "").strip(),
                    "distance": _as_float(item.findtext("dist")
                                          or item.findtext("dist_h-a")
                                          or item.findtext("centdist")),
                })
        sites.append({"ligand": hetid, "n": len(rows), "interactions": rows})

    if not sites:
        return {"ligand": ligand_ccd or "", "interactions": [], "counts": {}, "n": 0}
    # The site with the most interactions is the one being scored; a crystal
    # often carries buffer components PLIP profiles too.
    best = max(sites, key=lambda s: s["n"])
    counts: dict[str, int] = {}
    for row in best["interactions"]:
        counts[row["type"]] = counts.get(row["type"], 0) + 1
    return {"ligand": best["ligand"], "interactions": best["interactions"],
            "counts": counts, "n": best["n"]}


def fingerprint(result: dict, mapping: Optional[dict[int, int]] = None) -> set[tuple[str, int]]:
    """(interaction type, residue number) pairs, in reference numbering.

    `mapping` translates this structure's residue numbers into the reference's.
    Residues with no counterpart are dropped rather than compared by raw number:
    1M17 numbers EGFR from the mature protein and UniProt from the precursor, 24
    apart, and comparing raw numbers would make every contact look lost.
    """
    out: set[tuple[str, int]] = set()
    for row in result.get("interactions", []):
        num = row["resnr"]
        if mapping is not None:
            if num not in mapping:
                continue
            num = mapping[num]
        out.add((row["type"], num))
    return out


def jaccard(a: set, b: set) -> Optional[float]:
    if not a and not b:
        return None
    union = a | b
    return round(len(a & b) / len(union), 3) if union else None


def comparison_table(predicted: set[tuple[str, int]], reference: set[tuple[str, int]],
                     names: Optional[dict[int, str]] = None) -> list[dict]:
    """One row per residue, saying which interactions each structure makes."""
    names = names or {}
    residues = sorted({num for _, num in predicted | reference})
    rows = []
    for num in residues:
        p = sorted(t for t, n in predicted if n == num)
        r = sorted(t for t, n in reference if n == num)
        rows.append({
            "residue": num,
            "name": names.get(num, ""),
            "predicted": p,
            "reference": r,
            "shared": sorted(set(p) & set(r)),
            "status": "match" if set(p) & set(r) else ("predicted only" if p else "missed"),
        })
    return rows


# ---------------------------------------------------------------------------
# PandaMap
# ---------------------------------------------------------------------------

def pandamap(structure: str | Path, out_png: str | Path,
             ligand_ccd: Optional[str] = None) -> dict[str, Any]:
    """The 2D interaction map and PandaMap's own empirical dG."""
    try:
        from pandamap import HybridProtLigMapper
    except ImportError:
        return {"error": "PandaMap is not installed on this server."}
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    try:
        mapper = HybridProtLigMapper(str(structure), ligand_resname=ligand_ccd)
        # DSSP is an external binary the droplet does not carry; PandaMap's own
        # solvent-accessibility fallback is used instead.
        mapper.run_analysis(output_file=str(out_png), use_dssp=False)
        affinity = mapper.estimate_binding_affinity()
    except Exception as exc:
        return {"error": f"PandaMap could not read {Path(structure).name}: {exc}"}
    return {
        "png": out_png.name,
        "dg": _as_float((affinity or {}).get("dG_estimated")) if isinstance(affinity, dict) else None,
        "interpretation": (affinity or {}).get("interpretation", "") if isinstance(affinity, dict) else "",
        "note": (affinity or {}).get("note", "") if isinstance(affinity, dict) else "",
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def as_pdb(source: str | Path, dest: str | Path, chain: Optional[str] = None,
           keep_ligand: bool = True) -> Path:
    """Write `source` as a PDB PLIP will accept, optionally one chain only."""
    import gemmi

    st = gemmi.read_structure(str(source))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_waters()
    if chain:
        for model in st:
            for name in [c.name for c in model if c.name != chain]:
                model.remove_chain(name)
    if not keep_ligand:
        st.remove_ligands_and_waters()
    st.setup_entities()
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(st.make_pdb_string())
    return dest


def _as_float(v) -> Optional[float]:
    try:
        return round(float(v), 3)
    except (TypeError, ValueError):
        return None
