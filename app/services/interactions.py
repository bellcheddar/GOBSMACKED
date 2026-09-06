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

import contextlib

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
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  encoding="utf-8", errors="replace", timeout=timeout)
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
                    # PLIP reports both ends of every interaction. Keeping them
                    # is what lets the 3D view draw the same interactions the
                    # table lists, rather than a second opinion computed by the
                    # viewer that would quietly disagree with it.
                    "ligcoo": _coords(item.find("ligcoo")),
                    "protcoo": _coords(item.find("protcoo")),
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

# The instrument-panel tokens, for the 2D map.
PANEL = "#243044"
PANEL_2 = "#2c3a50"
LINE = "#3a4a63"
PHOS = "#5de1e6"
MUTED = "#9fb0c7"

# PandaMap's residue chips are drawn with facecolor='white' hard-coded, and the
# label inside each one carries a white outline, so no matplotlib setting can
# darken them. They are remapped afterwards, in the image, by luminance: the
# near-greyscale light pixels become panel colours and everything with a hue
# survives untouched. That only works if the text is NOT near-greyscale, which
# is why the labels are drawn in phosphor rather than in white.
DARK_STYLE = {
    "figure.facecolor": PANEL, "axes.facecolor": PANEL, "savefig.facecolor": PANEL,
    "axes.edgecolor": LINE, "axes.titlecolor": MUTED,
    "legend.facecolor": PANEL_2, "legend.edgecolor": LINE, "legend.labelcolor": MUTED,
    "text.color": PHOS,
}
# Guard rails on the remap: if PandaMap's palette ever changes, a mask that
# catches almost nothing or almost everything means the assumption has broken,
# and the untouched figure is better than a mangled one.
REMAP_MIN, REMAP_MAX = 0.01, 0.40
PANDAMAP_DPI = 150


def pandamap(structure: str | Path, out_png: str | Path,
             ligand_ccd: Optional[str] = None) -> dict[str, Any]:
    """The 2D interaction map, themed to the page, and PandaMap's empirical dG."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from pandamap import HybridProtLigMapper
    except ImportError:
        return {"error": "PandaMap is not installed on this server."}
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    try:
        with plt.rc_context(DARK_STYLE), _pandamap_layout(plt) as captured:
            mapper = HybridProtLigMapper(str(structure), ligand_resname=ligand_ccd)
            # DSSP is an external binary the droplet does not carry; PandaMap's
            # own solvent-accessibility fallback is used instead.
            mapper.run_analysis(output_file=str(out_png), use_dssp=False, dpi=PANDAMAP_DPI)
            affinity = mapper.estimate_binding_affinity()
    except Exception as exc:
        return {"error": f"PandaMap could not read {Path(structure).name}: {exc}"}

    remapped = _darken_light_greys(out_png)
    from PIL import Image

    try:
        with Image.open(out_png) as image:
            size = list(image.size)
    except Exception:                            # noqa: BLE001 - the map still stands
        size = None
    return {
        "png": out_png.name,
        "themed": remapped,
        "size": size,
        "hotspots": captured.get("spots") or [],
        "dg": _as_float((affinity or {}).get("dG_estimated")) if isinstance(affinity, dict) else None,
        "interpretation": (affinity or {}).get("interpretation", "") if isinstance(affinity, dict) else "",
        "note": (affinity or {}).get("note", "") if isinstance(affinity, dict) else "",
    }


# What each thing on the map means, for the tooltips. Keyed on the labels
# PandaMap writes, so a label it stops writing simply loses its tooltip rather
# than showing the wrong one.
INTERACTION_HELP = {
    "Hydrogen Bond": "A hydrogen bond: donor and acceptor within about 3.5 A and roughly "
                     "in line. The interaction most likely to be conserved between a "
                     "prediction and the crystal.",
    "Carbon-Pi": "A C-H group pointing into the face of an aromatic ring. Weak "
                 "individually, common in kinase pockets.",
    "Donor-Pi": "A hydrogen-bond donor pointing into an aromatic face rather than at a "
                "lone pair.",
    "Amide-Pi": "A backbone or side-chain amide stacked against an aromatic ring.",
    "Hydrophobic": "Non-polar contact: two carbon-rich surfaces close enough to exclude "
                   "water. Most of the buried surface of a drug-sized ligand.",
    "Ionic": "A salt bridge: full charges of opposite sign within about 4 A.",
    "Alkyl-Pi": "An aliphatic side chain packed against an aromatic ring face.",
    "Attractive Charge": "Charged groups of opposite sign close enough to attract without "
                         "meeting the geometry for a salt bridge.",
    "Pi-Cation": "A positive charge, usually a lysine or arginine, over an aromatic face.",
    "Unidirectional interaction": "The arrow points from donor to acceptor.",
    "Bidirectional interaction": "Both partners donate and accept.",
    "Solvent accessible": "The shaded halo marks residues still exposed to water rather "
                          "than buried against the ligand.",
}
AMINO_ACIDS = {
    "ALA": "alanine", "ARG": "arginine", "ASN": "asparagine", "ASP": "aspartate",
    "CYS": "cysteine", "GLN": "glutamine", "GLU": "glutamate", "GLY": "glycine",
    "HIS": "histidine", "ILE": "isoleucine", "LEU": "leucine", "LYS": "lysine",
    "MET": "methionine", "PHE": "phenylalanine", "PRO": "proline", "SER": "serine",
    "THR": "threonine", "TRP": "tryptophan", "TYR": "tyrosine", "VAL": "valine",
}


def _tooltip_for(label: str) -> Optional[str]:
    """The sentence a label deserves, or nothing if it explains itself."""
    text = label.strip()
    if text in INTERACTION_HELP:
        return f"{text}. {INTERACTION_HELP[text]}"
    parts = text.replace("\n", " ").split()
    if len(parts) == 2 and parts[0].upper() in AMINO_ACIDS and parts[1].isdigit():
        name = AMINO_ACIDS[parts[0].upper()]
        return (f"{parts[0].title()}{parts[1]}: {name} lining the site, numbered as the "
                f"crystal numbers it.")
    return None


def _hotspots(fig, dpi: float, pad: float = 0.1) -> list[dict]:
    """Where every label ended up, as fractions of the image that gets saved.

    Fractions rather than pixels because the page scales the picture to the
    panel, and the overlay has to scale with it. Everything is computed in
    inches: a Text artist's window extent is in display units at the FIGURE's
    dpi, which is not the dpi savefig is about to use, and mixing the two puts
    every hotspot in the wrong place by the ratio between them.
    """
    renderer = fig.canvas.get_renderer()
    tight = fig.get_tightbbox(renderer)          # inches
    left, bottom = tight.x0 - pad, tight.y0 - pad
    width, height = tight.width + 2 * pad, tight.height + 2 * pad
    if width <= 0 or height <= 0:
        return []

    spots = []
    for text in fig.findobj(match=lambda a: a.__class__.__name__ == "Text"):
        label = (text.get_text() or "").strip()
        if not label:
            continue
        help_text = _tooltip_for(label)
        if not help_text:
            continue
        try:
            box = text.get_window_extent(renderer)
        except Exception:                        # noqa: BLE001 - an artist with no extent
            continue
        x0, x1 = box.x0 / fig.dpi, box.x1 / fig.dpi
        y0, y1 = box.y0 / fig.dpi, box.y1 / fig.dpi
        spots.append({
            "label": label,
            "help": help_text,
            # Top-down for CSS, bottom-up in matplotlib.
            "x": round((x0 - left) / width, 5),
            "y": round(1 - (y1 - bottom) / height, 5),
            "w": round((x1 - x0) / width, 5),
            "h": round((y1 - y0) / height, 5),
        })
    return spots


@contextlib.contextmanager
def _pandamap_layout(plt):
    """Three changes to how PandaMap draws, made while it draws.

    None of them is reachable through its API: `visualize` takes a `title` but
    calls `plt.title` in both branches of the `if`, so there is no way to ask for
    none, and neither the legend's placement nor the figure size is exposed at
    all. `run_analysis` forwards neither. So the three calls are wrapped for the
    duration and restored afterwards, which is narrower than vendoring the
    function and re-writing it.

    * The title goes. "Protein-Ligand Interactions: superposed_md_final" is the
      name of a working file, and the panel it sits in is already headed
      Interaction map.
    * The legend moves out of the axes to the right. Inside, at `upper right`,
      it sits on top of the diagram and covers whichever residues happen to be
      there; outside, the diagram gets the whole square and the legend gets a
      column of its own.
    * The figure grows and turns landscape, because a square figure with a
      column beside it is no longer square.
    """
    original_title, original_subplots = plt.title, plt.subplots
    from matplotlib.axes import Axes

    original_legend = Axes.legend
    original_savefig = plt.savefig
    captured: dict = {"spots": [], "size": None}

    def no_title(*args, **kwargs):
        return None

    def wider(*args, **kwargs):
        # Smaller canvas, same point sizes: PandaMap lays the diagram out in
        # data coordinates, so shrinking the figure does not shrink the drawing,
        # it makes every label a larger fraction of the image. That is the only
        # lever on the residue label size, which PandaMap hard-codes in its own
        # ax.text calls.
        kwargs["figsize"] = (13, 9.0)
        return original_subplots(*args, **kwargs)

    def legend_at_the_side(self, *args, **kwargs):
        # PandaMap's first handle carries the same label as the legend's title,
        # so the heading "Interacting structural groups" is immediately followed
        # by an entry saying "Interacting structural groups". The duplicate goes;
        # the swatch it carried explains itself, since every residue box on the
        # diagram is drawn that way.
        title = kwargs.get("title")
        handles = kwargs.get("handles") or (args[0] if args else None)
        if handles and title:
            kept = [h for h in handles
                    if getattr(h, "get_label", lambda: "")() != title]
            if kept and len(kept) < len(handles):
                kwargs["handles"] = kept
                args = ()
        # Bigger type, tighter rows. Spacing was carrying the column's height
        # before, which is a waste of the space a larger font could have used:
        # thirteen entries at 19 point with 1.1 spacing stand as tall as the
        # diagram and can actually be read.
        kwargs.update(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False,
                      fontsize=19, title_fontsize=21, labelspacing=1.1,
                      handletextpad=0.9, borderaxespad=0.0)
        return original_legend(self, *args, **kwargs)

    def savefig_and_measure(*args, **kwargs):
        # The last moment the figure exists: visualize() saves and immediately
        # closes it, so anything the overlay needs has to be read here.
        try:
            figure = plt.gcf()
            figure.canvas.draw()
            captured["spots"] = _hotspots(figure, kwargs.get("dpi") or figure.dpi)
        except Exception:                        # noqa: BLE001 - tooltips are not the map
            captured["spots"] = []
        return original_savefig(*args, **kwargs)

    plt.title, plt.subplots = no_title, wider
    Axes.legend, plt.savefig = legend_at_the_side, savefig_and_measure
    try:
        yield captured
    finally:
        plt.title, plt.subplots = original_title, original_subplots
        Axes.legend, plt.savefig = original_legend, original_savefig


def _darken_light_greys(path: Path) -> bool:
    """Turn the white chips and outlines into panel colours, in place.

    Only pixels that are both light and close to greyscale are touched, so every
    coloured element (the interaction lines, the markers, the phosphor labels)
    is left exactly as PandaMap drew it. Antialiased edges land smoothly because
    the replacement is a ramp rather than a single colour.
    """
    try:
        import numpy as np
        from PIL import Image
    except ImportError:
        return False
    try:
        image = Image.open(path).convert("RGB")
        a = np.asarray(image).astype(np.float32) / 255.0
        mx, mn = a.max(axis=2), a.min(axis=2)
        saturation = np.where(mx > 0, (mx - mn) / np.maximum(mx, 1e-6), 0.0)
        mask = (saturation < 0.15) & (mx > 0.72)
        fraction = float(mask.mean())
        if not (REMAP_MIN <= fraction <= REMAP_MAX):
            return False
        low = np.array([0x3a, 0x4a, 0x63], dtype=np.float32) / 255.0
        high = np.array([0x2c, 0x3a, 0x50], dtype=np.float32) / 255.0
        t = np.clip((mx - 0.72) / 0.28, 0.0, 1.0)[..., None]
        a[mask] = (low * (1.0 - t) + high * t)[mask]
        Image.fromarray((a * 255).astype(np.uint8)).save(path, optimize=True)
        return True
    except Exception:
        return False


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


def _coords(node) -> Optional[list[float]]:
    if node is None:
        return None
    try:
        return [float(node.findtext(axis)) for axis in ("x", "y", "z")]
    except (TypeError, ValueError):
        return None


def _as_float(v) -> Optional[float]:
    try:
        return round(float(v), 3)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Drawing the interactions in 3D
# ---------------------------------------------------------------------------

# One colour per interaction type, matching the app's tokens. These are the
# colours the viewer uses; the legend on the page is generated from this map so
# the two cannot disagree.
INTERACTION_COLOURS = {
    "hbond": 0x5de1e6,
    "hydrophobic": 0x9fb0c7,
    "water bridge": 0x7ec8e2,
    "salt bridge": 0xffb454,
    "pi stacking": 0xc39cff,
    "pi cation": 0xff8a5c,
    "halogen bond": 0x7ee2a8,
    "metal": 0xff5c5c,
}

# Dashes per interaction. Mol*'s viewer build exports no shape builder, so the
# lines are drawn as structures: each dash is a two-atom fragment with a CONECT
# record, and a run of them along the interaction vector reads as a dashed line.
DASHES = 4
DASH_FRACTION = 0.55        # of each segment that is drawn rather than gap


def write_interaction_lines(result: dict, out_dir: Path, prefix: str) -> list[dict]:
    """One small PDB per interaction type present, for the viewer to load.

    Returns a list of {type, file, colour, count} so the page can build both the
    toggle and the legend from the same data.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    by_type: dict[str, list[dict]] = {}
    for row in result.get("interactions", []):
        if not row.get("ligcoo") or not row.get("protcoo"):
            continue
        by_type.setdefault(row["type"], []).append(row)

    written = []
    for kind, rows in sorted(by_type.items()):
        lines = []
        serial = 1
        conect = []
        for row in rows:
            for start, end in _dash_segments(row["protcoo"], row["ligcoo"]):
                for point in (start, end):
                    lines.append(
                        f"HETATM{serial:5d}  C   DSH X{len(conect) + 1:4d}    "
                        f"{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}  1.00  0.00           C"
                    )
                    serial += 1
                conect.append((serial - 2, serial - 1))
        if not conect:
            continue
        for a, b in conect:
            lines.append(f"CONECT{a:5d}{b:5d}")
        lines.append("END")
        name = f"{prefix}_{kind.replace(' ', '_')}.pdb"
        (out_dir / name).write_text("\n".join(lines) + "\n")
        written.append({"type": kind, "file": name,
                        "colour": INTERACTION_COLOURS.get(kind, 0x9fb0c7),
                        "count": len(rows)})
    return written


def _dash_segments(start, end, dashes: int = DASHES):
    """Split a line into `dashes` drawn pieces with gaps between them."""
    ax, ay, az = start
    bx, by, bz = end
    for i in range(dashes):
        t0 = i / dashes
        t1 = t0 + DASH_FRACTION / dashes
        yield (
            [ax + (bx - ax) * t0, ay + (by - ay) * t0, az + (bz - az) * t0],
            [ax + (bx - ax) * t1, ay + (by - ay) * t1, az + (bz - az) * t1],
        )


def write_pocket_sticks(structure: str | Path, residues, dest: Path,
                        chain: Optional[str] = None) -> Optional[Path]:
    """Just the pocket residues, as their own file.

    The Complex view wants them drawn as sticks over the cartoon, and Mol*'s
    viewer build has no selection language to carve them out of a loaded
    structure. Writing them as a second small structure is the same route the
    interaction lines take, and it costs a few kilobytes.
    """
    import gemmi

    wanted = set()
    for item in residues or []:
        tail = str(item).split(":")[-1]
        try:
            wanted.add(int(tail))
        except ValueError:
            continue
    if not wanted:
        return None

    st = gemmi.read_structure(str(structure))
    st.setup_entities()
    st.remove_waters()
    st.remove_hydrogens()
    out = gemmi.Structure()
    model = gemmi.Model("1")
    for ch in st[0]:
        if chain and ch.name != chain:
            continue
        keep = [res.clone() for res in ch if res.seqid.num in wanted]
        if not keep:
            continue
        new_chain = gemmi.Chain(ch.name)
        for res in keep:
            new_chain.add_residue(res)
        model.add_chain(new_chain)
    if not len(model):
        return None
    out.add_model(model)
    out.setup_entities()
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(out.make_pdb_string())
    return dest
