"""The bundle's summarise stage, on the traps a real trajectory sprang.

These are unit tests on the selection and framing logic rather than on a
trajectory, because the failures were both in deciding WHAT to measure, not in
the measuring.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bundle_template"))

mdtraj = pytest.importorskip("mdtraj", reason="mdtraj is a bundle dependency")

from gobsmacked_run import summarise  # noqa: E402


def topology_with(ligand_name: str):
    """A two-residue topology: one alanine and one ligand of six heavy atoms."""
    import mdtraj as md

    top = md.Topology()
    chain = top.add_chain()
    ala = top.add_residue("ALA", chain)
    for name in ("N", "CA", "C", "O", "CB"):
        top.add_atom(name, md.element.carbon if name != "N" else md.element.nitrogen, ala)
    lig = top.add_residue(ligand_name, chain)
    for i in range(6):
        top.add_atom(f"C{i}", md.element.carbon, lig)
    top.add_atom("H1", md.element.hydrogen, lig)
    return top


@pytest.mark.parametrize("name", ["LIG", "UNK", "UNL", "AQ4"])
def test_the_ligand_is_found_whatever_it_is_called(name):
    """UNK is the one that matters: MDTraj calls it a protein residue, so
    `not protein` selects nothing and the ligand series comes back empty with
    no error anywhere."""
    top = topology_with(name)
    selected = summarise.select_ligand(top)
    assert len(selected) == 6, f"{name}: expected the six heavy ligand atoms"
    assert all(top.atom(int(i)).residue.name == name for i in selected)
    assert all(top.atom(int(i)).element.symbol != "H" for i in selected)


def test_unk_really_does_look_like_protein_to_mdtraj():
    """The premise of the fix, asserted rather than assumed: if MDTraj ever
    stops calling UNK a protein residue, this test says so and the workaround
    can be revisited."""
    top = topology_with("UNK")
    unk = [r for r in top.residues if r.name == "UNK"][0]
    assert unk.is_protein, "MDTraj no longer classifies UNK as protein"
    assert len(top.select("not protein and not water")) == 0


def test_a_water_or_ion_is_never_taken_for_the_ligand():
    top = topology_with("HOH")
    assert len(summarise.select_ligand(top)) == 0


def test_rmsf_is_reported_with_the_frame_it_was_measured_in():
    """md.rmsf does not align: the alignment basis changes the answer by a
    factor of twenty, so the archive records which one was used."""
    import mdtraj as md

    top = topology_with("LIG")
    xyz = np.zeros((4, top.n_atoms, 3), dtype=np.float32)
    for frame in range(4):
        xyz[frame] = np.arange(top.n_atoms * 3).reshape(top.n_atoms, 3) * 0.01
        xyz[frame, 1] += frame * 0.02
    traj = md.Trajectory(xyz, top)
    ca = np.array([a.index for a in top.atoms if a.name == "CA"])
    result = summarise.rmsf(traj, ca, top)
    assert result["aligned_on"] == "all protein Ca"
    assert len(result["residues"]) == len(result["values"])
