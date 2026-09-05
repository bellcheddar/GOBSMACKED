"""Stage 2: prepare the receptor and build the ligand in 3D.

PDBFixer does the receptor: missing heavy atoms, missing residues inside the
chain, hydrogens at the campaign's pH. RDKit does the ligand: ETKDG conformer,
MMFF minimisation, one molecule in `ligand.sdf`.

The author numbering is preserved and written to `numbering.json`, because
every measurement the server makes later is expressed in it. PDBFixer renumbers
silently when it adds residues, and a shifted numbering would put the gatekeeper
on the wrong residue in every panel downstream.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def run(campaign: dict, work: Path, results: Path, log) -> dict[str, Any]:
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer
    from rdkit import Chem
    from rdkit.Chem import AllChem

    warnings: list[str] = []
    ligand_cfg = campaign.get("ligand") or {}
    ph = float(ligand_cfg.get("protonation_ph", 7.4))

    model = results / "model_apo.pdb"
    receptor = work / "receptor.pdb"
    log(f"prep: PDBFixer on {model.name} at pH {ph}")

    original = residue_numbers(model.read_text())
    fixer = PDBFixer(filename=str(model))
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    # Missing residues at the chain ends are absent from the construct, not from
    # the model: building them adds flopping tails that cost MD time and add
    # nothing at the binding site.
    chains = list(fixer.topology.chains())
    for key in list(fixer.missingResidues):
        chain_index, position = key
        chain = chains[chain_index]
        if position == 0 or position == len(list(chain.residues())):
            del fixer.missingResidues[key]
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with open(receptor, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

    prepared = residue_numbers(receptor.read_text())
    (results / "numbering.json").write_text(json.dumps({
        "original_first": original[0] if original else None,
        "original_last": original[-1] if original else None,
        "prepared_first": prepared[0] if prepared else None,
        "prepared_last": prepared[-1] if prepared else None,
        "preserved": original[:1] == prepared[:1],
    }, indent=2))
    if original[:1] != prepared[:1]:
        warnings.append(
            f"PDBFixer renumbered the receptor: it started at {original[0]} and now starts at "
            f"{prepared[0]}. Residue numbers in this run refer to the prepared file."
        )

    smiles = (ligand_cfg.get("smiles") or "").strip()
    if not smiles:
        raise RuntimeError("The campaign has no ligand SMILES.")
    log(f"prep: RDKit conformer for {ligand_cfg.get('name') or 'the ligand'}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"RDKit cannot parse the campaign SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D           # a fixed seed, so a rerun redocks the same ligand
    if AllChem.EmbedMolecule(mol, params) != 0:
        # Random-coordinate embedding is slower but succeeds on macrocycles and
        # heavily bridged systems where the distance-geometry start fails.
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise RuntimeError("RDKit could not generate a 3D conformer for this ligand.")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
        warnings.append("MMFF has no parameters for this ligand; UFF was used for the conformer. "
                        "The docked pose is unaffected, the starting geometry is coarser.")
    mol.SetProp("_Name", ligand_cfg.get("name") or "ligand")
    writer = Chem.SDWriter(str(work / "ligand.sdf"))
    writer.write(mol)
    writer.close()

    log(f"prep: receptor and ligand written ({mol.GetNumHeavyAtoms()} heavy atoms)")
    return {"warnings": warnings, "heavy_atoms": mol.GetNumHeavyAtoms()}


def residue_numbers(pdb_text: str) -> list[int]:
    numbers = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                numbers.append(int(line[22:26]))
            except ValueError:
                continue
    return numbers
