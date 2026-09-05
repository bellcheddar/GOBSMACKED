"""Stage 4: minimise, equilibrate and run.

Amber14 for the protein, OpenFF Sage for the ligand, TIP3P water with 0.15 M
NaCl in a 10 A padded box. Restraints hold the solute while fresh water relaxes
around it and are released over the equilibration, so the pocket is not shaken
apart before production starts.

The ligand's parameters come from the SDF PandaDock wrote, not from the PDB
complex: a PDB carries no bond orders, and every ligand force field starts from
the chemistry.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any, Optional

# Hydrogen mass repartitioning at 4 amu with a 2 fs step. HMR permits 4 fs, and
# 4 fs with a ligand that has a fast internal mode is where an MD run quietly
# stops conserving energy. The extra factor of two is not worth the risk here.
HYDROGEN_MASS_AMU = 4.0
TIMESTEP_FS = 2.0
FRICTION_PER_PS = 1.0
TEMPERATURE_K = 300.0
PRESSURE_BAR = 1.0
PADDING_NM = 1.0
IONIC_STRENGTH_M = 0.15
RESTRAINT_K = 1000.0            # kJ/mol/nm^2 on solute heavy atoms, released in steps
RESTRAINT_STEPS = 5
# What the ligand residue is called in every file this stage writes.
LIGAND_RESIDUE_NAME = "LIG"


def run(campaign: dict, work: Path, results: Path, log) -> dict[str, Any]:
    # openff-toolkit pulls in torch and the whole SMIRNOFF machinery, which on a
    # cold cache and a busy machine took sixteen minutes to import on the first
    # real run: long enough to look like a hang with nothing on screen. Say what
    # is happening before paying for it.
    log("md: importing OpenMM and the OpenFF toolkit (slow on a cold cache)")
    import numpy as np
    import openmm
    from openmm import app, unit
    from openmmforcefields.generators import SystemGenerator
    from openff.toolkit import Molecule

    md_cfg = campaign.get("md") or {}
    warnings: list[str] = []

    receptor = app.PDBFile(str(work / "receptor.pdb"))
    ligand_sdf = top_pose_sdf(results / "poses" / "poses.sdf", work / "pose1.sdf", log)
    ligand = Molecule.from_file(str(ligand_sdf), file_format="sdf", allow_undefined_stereo=True)

    log("md: building the system (Amber14 + "
        f"{md_cfg.get('ligand_forcefield', 'openff-2.1.0')})")
    generator = SystemGenerator(
        forcefields=["amber14-all.xml", "amber14/tip3p.xml"],
        small_molecule_forcefield=md_cfg.get("ligand_forcefield", "openff-2.1.0"),
        molecules=[ligand],
        forcefield_kwargs={
            "constraints": app.HBonds,
            "rigidWater": True,
            "removeCMMotion": True,
            "hydrogenMass": HYDROGEN_MASS_AMU * unit.amu,
        },
    )

    modeller = app.Modeller(receptor.topology, receptor.positions)
    ligand_topology = ligand.to_topology().to_openmm()
    # OpenFF's topology names the ligand residue UNK, and MDTraj classifies UNK
    # as a protein residue: `not protein` then selects nothing, and the ligand
    # RMSD series comes back empty with no error anywhere. Name it once, here.
    for residue in ligand_topology.residues():
        residue.name = LIGAND_RESIDUE_NAME
    ligand_positions = ligand.conformers[0].to_openmm()
    modeller.add(ligand_topology, ligand_positions)
    solute_atoms = modeller.topology.getNumAtoms()

    log(f"md: solvating {solute_atoms} solute atoms, {PADDING_NM * 10:.0f} A padding")
    modeller.addSolvent(generator.forcefield, model="tip3p",
                        padding=PADDING_NM * unit.nanometer,
                        ionicStrength=IONIC_STRENGTH_M * unit.molar, neutralize=True)
    total_atoms = modeller.topology.getNumAtoms()
    log(f"md: {total_atoms} atoms in the box")

    system = generator.create_system(modeller.topology, molecules=[ligand])
    restraint = add_restraints(system, modeller.topology, modeller.positions, solute_atoms)
    system.addForce(openmm.MonteCarloBarostat(PRESSURE_BAR * unit.bar,
                                              TEMPERATURE_K * unit.kelvin))

    integrator = openmm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin, FRICTION_PER_PS / unit.picosecond,
        TIMESTEP_FS * unit.femtoseconds)
    platform, properties = choose_platform(md_cfg.get("platform", "auto"), log, warnings)
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    steps_minimise = int(md_cfg.get("minimise_steps", 5000))
    log(f"md: minimising, {steps_minimise} steps")
    simulation.minimizeEnergy(maxIterations=steps_minimise)
    write_complex(simulation, results / "complex_min.pdb", solute_atoms)

    equil_ps = float(md_cfg.get("equilibration_ps", 100))
    equil_steps = int(equil_ps * 1000 / TIMESTEP_FS)
    if equil_steps:
        log(f"md: equilibrating {equil_ps:.0f} ps, releasing restraints in {RESTRAINT_STEPS} steps")
        simulation.context.setVelocitiesToTemperature(TEMPERATURE_K * unit.kelvin)
        per_step = max(1, equil_steps // RESTRAINT_STEPS)
        for index in range(RESTRAINT_STEPS):
            k = RESTRAINT_K * (1.0 - index / RESTRAINT_STEPS)
            simulation.context.setParameter("k_restraint", k)
            simulation.step(per_step)
        simulation.context.setParameter("k_restraint", 0.0)

    production_ps = float(md_cfg.get("production_ps", 1000))
    interval_ps = max(1.0, float(md_cfg.get("frame_interval_ps", 10)))
    production_steps = int(production_ps * 1000 / TIMESTEP_FS)
    interval_steps = int(interval_ps * 1000 / TIMESTEP_FS)
    traj_dir = results / "traj"
    traj_dir.mkdir(parents=True, exist_ok=True)

    # The trajectory holds the solute only. Water is 90 % of the atoms, none of
    # the analysis reads it, and a 1 ns solvated DCD is hundreds of megabytes to
    # upload for nothing.
    solute_indices = list(range(solute_atoms))
    if production_steps:
        log(f"md: production {production_ps:.0f} ps, frame every {interval_ps:.0f} ps")
        simulation.reporters.append(
            # enforcePeriodicBox=False: the subset written here IS the molecule
            # being measured, and wrapping it into the box splits it whenever it
            # drifts across a face. That is invisible in any single frame and
            # shows up later as a residue that appears to fluctuate by a box
            # length.
            app.DCDReporter(str(traj_dir / "traj.dcd"), interval_steps,
                            enforcePeriodicBox=False, atomSubset=solute_indices))
        simulation.reporters.append(app.StateDataReporter(
            str(results / "logs" / "md.csv"), interval_steps * 5, step=True, time=True,
            potentialEnergy=True, temperature=True, density=True, speed=True))
        (results / "logs").mkdir(parents=True, exist_ok=True)
        simulation.step(production_steps)

    write_complex(simulation, results / "complex_md_final.pdb", solute_atoms)
    write_complex(simulation, traj_dir / "topology.pdb", solute_atoms)
    log("md: done")
    return {"warnings": warnings, "atoms": total_atoms, "solute_atoms": solute_atoms,
            "frames": production_steps // interval_steps if interval_steps else 0}


def add_restraints(system, topology, positions, solute_atoms: int):
    """Positional restraints on solute heavy atoms, releasable by a parameter.

    `periodicdistance` rather than a plain distance, so an atom that crosses a
    periodic boundary is not yanked back across the whole box.
    """
    import openmm
    from openmm import unit

    restraint = openmm.CustomExternalForce(
        "k_restraint*periodicdistance(x, y, z, x0, y0, z0)^2")
    restraint.addGlobalParameter("k_restraint",
                                 RESTRAINT_K * unit.kilojoules_per_mole / unit.nanometer ** 2)
    for name in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(name)
    atoms = list(topology.atoms())
    for index in range(min(solute_atoms, len(atoms))):
        if atoms[index].element is None or atoms[index].element.symbol == "H":
            continue
        restraint.addParticle(index, positions[index].value_in_unit(unit.nanometer))
    system.addForce(restraint)
    return restraint


def choose_platform(wanted: str, log, warnings: list[str]):
    """The fastest platform available, at a precision it will actually accept.

    Apple's OpenCL is single precision only, and asking it for `mixed` fails
    with "No compatible OpenCL platform is available": an error that names the
    platform rather than the property, reads as "you have no GPU", and drops the
    run onto the CPU at roughly a tenth of the speed with nothing in the output
    to say so.
    """
    import openmm

    order = ["CUDA", "OpenCL", "CPU"] if wanted in ("auto", "", None) else [wanted]
    for name in order:
        try:
            platform = openmm.Platform.getPlatformByName(name)
        except Exception:
            continue
        properties = {}
        if name in ("CUDA", "OpenCL"):
            supported = set(platform.getPropertyNames())
            key = "Precision"
            if key in supported:
                properties[key] = "mixed" if name == "CUDA" else "single"
        log(f"md: platform {name}" + (f" ({properties})" if properties else ""))
        if name == "CPU" and wanted in ("auto", "", None):
            warnings.append("No GPU platform was available, so MD ran on the CPU. Expect "
                            "roughly a tenth of the throughput.")
        return platform, properties
    return None, {}


def top_pose_sdf(poses: Path, dest: Path, log) -> Path:
    """The rank-1 pose alone, as an SDF with its bond orders intact."""
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(poses), removeHs=False)
    first = next((m for m in supplier if m is not None), None)
    if first is None:
        raise RuntimeError("No readable pose in poses.sdf.")
    writer = Chem.SDWriter(str(dest))
    writer.write(first)
    writer.close()
    return dest


def write_complex(simulation, dest: Path, solute_atoms: int) -> None:
    """The solute's current coordinates, without the water.

    Waters are stripped by rebuilding a Modeller over the solute atoms rather
    than by filtering the PDB afterwards: an atom-index filter over a topology
    that has been through addSolvent is exactly the kind of off-by-one that
    produces a file which opens fine and is wrong.
    """
    from openmm import app

    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=False)
    positions = state.getPositions()
    modeller = app.Modeller(simulation.topology, positions)
    solvent = [atom for index, atom in enumerate(simulation.topology.atoms())
               if index >= solute_atoms]
    if solvent:
        modeller.delete(solvent)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "w") as fh:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)
