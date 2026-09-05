"""Writing the run bundle.

A bundle is a self-contained directory that a user unpacks on a machine with a
GPU: the campaign, the starting structure if one was found, a pixi environment
and the five-stage runner. Nothing in it phones home, and nothing in it needs
this server again until the results come back.

The campaign file is the contract between the two halves of the app, so it is
written from one place and validated by `ingest.py` on the way back.
"""

from __future__ import annotations

import hashlib
import io
import re
import shutil
import tarfile
from pathlib import Path
from typing import Any, Iterable, Optional

import gemmi
import yaml

from .. import config

CAMPAIGN_VERSION = "1.0"

RESIDUE_RE = re.compile(r"^(?:([A-Za-z0-9]+):)?(-?\d+)$")


def _is_amino_acid(res: gemmi.Residue) -> bool:
    """gemmi keeps this on the tabulated component, not on the residue."""
    info = gemmi.find_tabulated_residue(res.name)
    return bool(info and info.is_amino_acid())


# ---------------------------------------------------------------------------
# Pocket geometry
# ---------------------------------------------------------------------------

def parse_residues(residues: Iterable[str], default_chain: str = "A") -> list[tuple[str, int]]:
    """["A:718", "745"] -> [("A", 718), ("A", 745)]. Unparseable entries are dropped."""
    out: list[tuple[str, int]] = []
    for item in residues:
        m = RESIDUE_RE.match(str(item).strip())
        if not m:
            continue
        chain, num = m.groups()
        out.append((chain or default_chain, int(num)))
    return out


def pocket_box(structure_path: str | Path, residues: Iterable[str],
               default_chain: str = "A") -> Optional[dict]:
    """Centroid and box side lengths for a residue selection.

    The box is the extent of the selected residues' heavy atoms plus
    `BOX_PADDING` on each side, floored at `BOX_MIN_SIDE`: a two-residue
    selection would otherwise produce a box too small for the ligand to rotate
    in, and docking would report failure that was really a setup mistake.
    """
    wanted = set(parse_residues(residues, default_chain))
    if not wanted:
        return None
    st = gemmi.read_structure(str(structure_path))
    st.setup_entities()
    st.remove_waters()
    coords: list[tuple[float, float, float]] = []
    found: set[tuple[str, int]] = set()
    for chain in st[0]:
        for res in chain:
            key = (chain.name, res.seqid.num)
            if key not in wanted:
                continue
            found.add(key)
            for atom in res:
                if atom.element == gemmi.Element("H"):
                    continue
                coords.append((atom.pos.x, atom.pos.y, atom.pos.z))
    if not coords:
        return None
    xs, ys, zs = zip(*coords)
    centre = [round(sum(v) / len(v), 3) for v in (xs, ys, zs)]
    box = [
        round(max(config.BOX_MIN_SIDE, (max(v) - min(v)) + 2 * config.BOX_PADDING), 1)
        for v in (xs, ys, zs)
    ]
    return {
        "center": centre,
        "box": box,
        "n_residues": len(found),
        "missing": sorted(f"{c}:{n}" for c, n in (wanted - found)),
    }


def residues_near_ligand(structure_path: str | Path, ccd: str,
                         radius: float = config.POCKET_RADIUS) -> list[str]:
    """Polymer residues with a heavy atom within `radius` of a named ligand.

    This is the "use the reference ligand's site" button on Panel 3, and the
    same definition superpose.py uses for the pocket, so the residues the user
    docks into are the residues the verification is measured on.
    """
    st = gemmi.read_structure(str(structure_path))
    st.setup_entities()
    st.remove_waters()
    ligand_atoms = [
        atom.pos for chain in st[0] for res in chain
        if res.name.upper() == ccd.upper() for atom in res
        if atom.element != gemmi.Element("H")
    ]
    if not ligand_atoms:
        return []
    hits: set[tuple[str, int]] = set()
    for chain in st[0]:
        for res in chain:
            if res.name.upper() == ccd.upper() or res.is_water():
                continue
            if not _is_amino_acid(res):
                continue
            for atom in res:
                if atom.element == gemmi.Element("H"):
                    continue
                if any(atom.pos.dist(p) <= radius for p in ligand_atoms):
                    hits.add((chain.name, res.seqid.num))
                    break
    return [f"{c}:{n}" for c, n in sorted(hits, key=lambda t: (t[0], t[1]))]


# ---------------------------------------------------------------------------
# The campaign file
# ---------------------------------------------------------------------------

def build_campaign(job_id: str, protein: dict, ligand: dict, pocket: dict,
                   reference: dict, docking: dict, md: dict,
                   owner_token: str = "", title: str = "") -> dict:
    """Assemble campaign.yaml's content. Key order is the order in the file."""
    return {
        "gobsmacked_version": CAMPAIGN_VERSION,
        "job_id": job_id,
        "title": title or "",
        "owner_token": owner_token,
        "protein": {
            "uniprot": protein.get("uniprot"),
            "sequence": protein.get("sequence", ""),
            "source_structure": protein.get("source_structure", "fold"),
            "source_id": protein.get("source_id"),
            "chain": protein.get("chain", "A"),
            "residue_range": protein.get("residue_range"),
            "family": protein.get("family", "other"),
        },
        "ligand": {
            "name": ligand.get("name", "ligand"),
            "smiles": ligand.get("smiles", ""),
            "protonation_ph": float(ligand.get("protonation_ph", 7.4)),
        },
        "pocket": {
            "method": pocket.get("method", "residues"),
            "residues": list(pocket.get("residues") or []),
            "center": pocket.get("center"),
            "box": pocket.get("box"),
        },
        "reference": {
            "pdb_id": reference.get("pdb_id"),
            "chain": reference.get("chain"),
            "ligand_ccd": reference.get("ligand_ccd"),
            "apo_pdb_id": reference.get("apo_pdb_id"),
        },
        "docking": {
            "mode": docking.get("mode", "hybrid"),
            "exhaustiveness": int(docking.get("exhaustiveness", 16)),
            "num_poses": int(docking.get("num_poses", 10)),
            "flexible_residues": docking.get("flexible_residues", "auto"),
        },
        "md": {
            "forcefield": md.get("forcefield", "amber14"),
            "ligand_forcefield": md.get("ligand_forcefield", "openff-2.1.0"),
            "minimise_steps": int(md.get("minimise_steps", 5000)),
            "equilibration_ps": int(md.get("equilibration_ps", 100)),
            "production_ps": int(md.get("production_ps", 1000)),
            "frame_interval_ps": int(md.get("frame_interval_ps", 10)),
            "platform": md.get("platform", "auto"),
        },
    }


def dump_campaign(campaign: dict) -> str:
    """campaign.yaml text, with the key order preserved rather than sorted."""
    return yaml.safe_dump(campaign, sort_keys=False, default_flow_style=False, width=100)


def campaign_sha256(campaign_yaml: str) -> str:
    return hashlib.sha256(campaign_yaml.encode("utf-8")).hexdigest()


# ---------------------------------------------------------------------------
# The archive
# ---------------------------------------------------------------------------

def write_bundle(job_id: str, campaign: dict, structure_path: Optional[str] = None,
                 dest_dir: Optional[Path] = None) -> Path:
    """Write `run_bundle_<job_id>.tar.gz` and return its path.

    The bundle template is copied verbatim: every run of every job gets exactly
    the code in this repository at deploy time, so a results archive can be
    traced back to a commit.
    """
    dest_dir = Path(dest_dir or (config.RUNS_DIR / job_id))
    dest_dir.mkdir(parents=True, exist_ok=True)
    root = f"run_bundle_{job_id}"
    archive = dest_dir / f"{root}.tar.gz"
    campaign_yaml = dump_campaign(campaign)
    (dest_dir / "campaign.yaml").write_text(campaign_yaml)

    with tarfile.open(archive, "w:gz") as tar:
        for src in sorted(config.BUNDLE_TEMPLATE_DIR.rglob("*")):
            if any(part in {"__pycache__", ".pixi"} for part in src.parts):
                continue
            if src.is_dir():
                continue
            tar.add(src, arcname=f"{root}/{src.relative_to(config.BUNDLE_TEMPLATE_DIR)}")
        _add_text(tar, f"{root}/campaign.yaml", campaign_yaml)
        if structure_path and Path(structure_path).exists():
            # The bundle's fold stage is skipped when this file is present, so
            # it always arrives as PDB regardless of what was fetched.
            pdb_text = _as_pdb(structure_path, campaign["protein"].get("chain"))
            _add_text(tar, f"{root}/model_apo.pdb", pdb_text)
    return archive


def _add_text(tar: tarfile.TarFile, arcname: str, text: str) -> None:
    data = text.encode("utf-8")
    info = tarfile.TarInfo(arcname)
    info.size = len(data)
    info.mode = 0o644
    tar.addfile(info, io.BytesIO(data))


def _as_pdb(path: str | Path, chain: Optional[str] = None) -> str:
    """Whatever was fetched (PDB or mmCIF), as a single-chain PDB with no waters.

    Docking and MD both want one protein chain and nothing else: a second copy
    in the asymmetric unit doubles the system size for no benefit, and the
    crystallisation additives are not part of the model being tested.
    """
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_waters()
    st.remove_ligands_and_waters()
    if chain:
        for model in st:
            for ch in [c.name for c in model if c.name != chain]:
                model.remove_chain(ch)
    st.setup_entities()
    return st.make_pdb_string()
