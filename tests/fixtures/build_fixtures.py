#!/usr/bin/env python3
"""Build the two test archives from real crystal structures.

A fixture made of noise would exercise the plumbing and none of the science, so
each archive is built from two genuine entries: one standing in for the
prediction, one for the crystal it is judged against.

    EGFR   4HJO (erlotinib) judged against 1M17 (erlotinib)
    beta2  5D5A (carazolol) judged against 2RH1 (carazolol)

The ligand is displaced to make a plausible docked pose, a twenty-frame
trajectory is generated around it, and `summarise.run` (the bundle's own code)
computes summary.json. So the fixtures test the bundle's last stage as well as
the server's first.

    python tests/fixtures/build_fixtures.py
"""

from __future__ import annotations

import csv
import json
import shutil
import sys
import tempfile
from pathlib import Path

import gemmi
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "bundle_template"))

FIXTURES = Path(__file__).resolve().parent
CACHE = FIXTURES / "_structures"
FRAMES = 20
FRAME_INTERVAL_PS = 10.0
LIGAND_RESNAME = "LIG"

CASES = [
    {
        "name": "egfr",
        "job_id": "gs_20260905_fixtureegfr",
        "prediction": {"pdb": "4HJO", "chain": "A", "ccd": "AQ4"},
        "reference": {"pdb": "1M17", "chain": "A", "ccd": "AQ4"},
        "uniprot": "P00533",
        "gene": "EGFR",
        "family": "kinase",
        "ligand_name": "erlotinib",
        "smiles": "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC",
        # A displacement large enough to be visible on the scorecard and small
        # enough to still be the same binding mode: a realistic good docking.
        "pose_shift": 0.9,
    },
    {
        "name": "beta2",
        "job_id": "gs_20260905_fixtureb2ar",
        "prediction": {"pdb": "5D5A", "chain": "A", "ccd": "CAU"},
        "reference": {"pdb": "2RH1", "chain": "A", "ccd": "CAU"},
        "uniprot": "P07550",
        "gene": "ADRB2",
        "family": "gpcr",
        "ligand_name": "carazolol",
        "smiles": "CC(C)NCC(O)COc1cccc2[nH]c3ccccc3c12",
        "pose_shift": 1.2,
    },
]


def fetch(pdb_id: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / f"{pdb_id}.cif"
    if not dest.exists():
        import requests

        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        dest.write_bytes(requests.get(url, timeout=60).content)
    return dest


def single_chain(path: Path, chain: str, ccd: str | None, keep_ligand: bool) -> gemmi.Structure:
    """One chain, no waters, no crystallisation additives, ligand optional."""
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_waters()
    for model in st:
        for name in [c.name for c in model if c.name != chain]:
            model.remove_chain(name)
    model = st[0]
    for ch in model:
        keep = []
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info and info.is_amino_acid():
                keep.append(res)
            elif keep_ligand and ccd and res.name.upper() == ccd.upper():
                res.name = LIGAND_RESNAME
                res.het_flag = "H"
                keep.append(res)
        # gemmi's Chain has no remove_residue; deleting by index while walking
        # the chain would skip residues, so the survivors are cloned first and
        # the chain is emptied from the end.
        keep = [res.clone() for res in keep]
        for index in range(len(ch) - 1, -1, -1):
            del ch[index]
        for res in keep:
            ch.add_residue(res)
    st.setup_entities()
    return st


def shift_ligand(st: gemmi.Structure, distance: float, seed: int = 0) -> gemmi.Structure:
    """Translate the ligand a fixed distance in a fixed random direction."""
    rng = np.random.default_rng(seed)
    direction = rng.normal(size=3)
    direction /= np.linalg.norm(direction)
    offset = direction * distance
    out = st.clone()
    for ch in out[0]:
        for res in ch:
            if res.name != LIGAND_RESNAME:
                continue
            for atom in res:
                atom.pos = gemmi.Position(atom.pos.x + offset[0], atom.pos.y + offset[1],
                                          atom.pos.z + offset[2])
    return out


def write_pdb(st: gemmi.Structure, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(st.make_pdb_string())


def ligand_sdf(complex_pdb: Path, smiles: str, dest: Path, n_poses: int = 5) -> list[dict]:
    """poses.sdf and the score rows, from the ligand in the complex."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    sys.path.insert(0, str(ROOT))
    from app.services.superpose import _ligand_mol

    mol = _ligand_mol(complex_pdb, LIGAND_RESNAME, smiles)
    if mol is None:
        raise RuntimeError(f"no ligand in {complex_pdb}")
    rows = []
    writer = Chem.SDWriter(str(dest))
    rng = np.random.default_rng(7)
    for rank in range(1, n_poses + 1):
        pose = Chem.Mol(mol)
        if rank > 1:
            # Decoys, so the pose file is not one molecule copied five times.
            conf = pose.GetConformer()
            offset = rng.normal(scale=0.6 * rank, size=3)
            for i in range(pose.GetNumAtoms()):
                p = conf.GetAtomPosition(i)
                conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(
                    p.x + offset[0], p.y + offset[1], p.z + offset[2]))
        score = round(-9.4 + 0.7 * (rank - 1), 3)
        pose.SetProp("_Name", f"pose{rank}")
        pose.SetProp("score", str(score))
        pose.SetProp("rank", str(rank))
        writer.write(pose)
        rows.append({"pose_id": f"pose{rank}", "score": score,
                     "gnn_affinity": round(7.8 - 0.3 * (rank - 1), 3), "rank": rank})
    writer.close()
    return rows


def write_trajectory(topology: Path, dest: Path, frames: int, drift: float) -> None:
    """A short trajectory: thermal noise on the protein, a slow drift on the ligand.

    Not a physical trajectory, and not pretending to be. It exists so the
    dynamics panels, the drift metric and the contact matrix have something with
    the right shape to run on.
    """
    import mdtraj as md

    frame = md.load(str(topology))
    ligand = frame.topology.select(f"resname {LIGAND_RESNAME}")
    rng = np.random.default_rng(11)
    coords = np.repeat(frame.xyz, frames, axis=0)
    for i in range(frames):
        coords[i] += rng.normal(scale=0.012, size=coords[i].shape)      # nm
        coords[i, ligand] += np.array([1, 0.4, -0.6]) * (drift / 10.0) * (i / max(1, frames - 1))
    traj = md.Trajectory(coords, frame.topology)
    traj.time = np.arange(frames) * FRAME_INTERVAL_PS
    traj.save_dcd(str(dest))


def build(case: dict) -> Path:
    from gobsmacked_run import schema, summarise
    from app.services.bundle import pocket_box, residues_near_ligand

    print(f"=== {case['name']}: {case['prediction']['pdb']} against {case['reference']['pdb']}")
    prediction_cif = fetch(case["prediction"]["pdb"])
    fetch(case["reference"]["pdb"])

    work = Path(tempfile.mkdtemp(prefix=f"gsfix_{case['name']}_"))
    results = work / "results"
    for sub in ("poses", "traj", "logs"):
        (results / sub).mkdir(parents=True, exist_ok=True)

    holo = single_chain(prediction_cif, case["prediction"]["chain"], case["prediction"]["ccd"], True)
    apo = single_chain(prediction_cif, case["prediction"]["chain"], case["prediction"]["ccd"], False)

    write_pdb(apo, results / "model_apo.pdb")
    write_pdb(holo, results / "complex_md_final.pdb")
    write_pdb(holo, results / "traj" / "topology.pdb")
    write_pdb(shift_ligand(holo, case["pose_shift"], seed=1), results / "complex_pose1.pdb")
    write_pdb(shift_ligand(holo, case["pose_shift"] * 0.55, seed=1), results / "complex_min.pdb")

    rows = ligand_sdf(results / "complex_md_final.pdb", case["smiles"],
                      results / "poses" / "poses.sdf")
    with open(results / "poses" / "scores.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["pose_id", "score", "gnn_affinity", "rank"])
        writer.writeheader()
        writer.writerows(rows)

    write_trajectory(results / "traj" / "topology.pdb", results / "traj" / "traj.dcd",
                     FRAMES, drift=0.35)

    pocket = residues_near_ligand(results / "complex_md_final.pdb", LIGAND_RESNAME)
    # A real campaign always carries a centre and box; without them the pose's
    # inside-the-box validity check has nothing to test against.
    box = pocket_box(results / "complex_md_final.pdb", pocket,
                     case["prediction"]["chain"])
    campaign = {
        "gobsmacked_version": "1.0",
        "job_id": case["job_id"],
        "title": f"{case['gene']} + {case['ligand_name']} fixture",
        "owner_token": "",
        "protein": {
            "uniprot": case["uniprot"], "sequence": "", "source_structure": "pdb",
            "source_id": case["prediction"]["pdb"], "chain": case["prediction"]["chain"],
            "residue_range": None, "family": case["family"], "gene": case["gene"],
        },
        "ligand": {"name": case["ligand_name"], "smiles": case["smiles"], "protonation_ph": 7.4},
        "pocket": {"method": "residues", "residues": pocket,
                   "center": box["center"], "box": box["box"]},
        "reference": {"pdb_id": case["reference"]["pdb"], "chain": case["reference"]["chain"],
                      "ligand_ccd": case["reference"]["ccd"], "apo_pdb_id": None},
        "docking": {"mode": "hybrid", "exhaustiveness": 16, "num_poses": 5,
                    "flexible_residues": "auto"},
        "md": {"forcefield": "amber14", "ligand_forcefield": "openff-2.1.0",
               "minimise_steps": 5000, "equilibration_ps": 100,
               "production_ps": FRAMES * FRAME_INTERVAL_PS,
               "frame_interval_ps": FRAME_INTERVAL_PS, "platform": "auto"},
    }
    # The sequence is fetched so the mode classifiers have the canonical
    # numbering to map KLIFS and GPCRdb positions through.
    from app.services.fetch import fetch_uniprot
    entry = fetch_uniprot(case["uniprot"])
    campaign["protein"]["sequence"] = entry["sequence"] if entry else ""
    (results / "campaign.yaml").write_text(yaml.safe_dump(campaign, sort_keys=False))

    def log(message: str) -> None:
        print("   ", message)

    summarise.run(campaign, work, results, log)
    (results / "logs" / "run.log").write_text("fixture archive: built by build_fixtures.py\n")
    schema.write_manifest(results, case["job_id"], results / "campaign.yaml",
                          {"fold": 0.0, "prep": 2.0, "dock": 41.0, "md": 63.0, "summarise": 6.0},
                          ["This archive is a test fixture built from crystal structures, "
                           "not the output of a real run."])
    missing = schema.check_complete(results)
    if missing:
        raise RuntimeError(f"fixture is incomplete: {missing}")

    dest = FIXTURES / f"{case['name']}_results.tar.gz"
    summarise.pack(results, dest, log)
    shutil.rmtree(work, ignore_errors=True)
    return dest


def main() -> int:
    for case in CASES:
        path = build(case)
        print(f"    {path.name}: {path.stat().st_size / 1024:.0f} kB\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
