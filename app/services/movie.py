"""The trajectory, played back: an MP4 of the pose moving inside the pocket.

The dynamics panel already reports what the trajectory did as five traces, and
the traces are the measurement. This is the same run shown rather than plotted,
because a ligand walking out of its site is obvious in two seconds of video and
takes a paragraph of RMSD to describe.

Three decisions worth stating, because each one changes what the clip means:

* **Every frame is superposed on the protein's own Ca atoms in frame 0**, so the
  camera sits on the protein and the ligand is the thing seen to move. Aligning
  on the pocket instead would hold the site still and make the whole protein
  swim, which is a different and much less readable claim.
* **The ligand is drawn in phosphor cyan**, the site's colour for the predicted
  state, and the residues lining its site in near-white. Everything else is the
  muted grey the site uses for context. Nothing here is coloured amber: no
  experimental structure appears in this clip, and amber means experiment
  everywhere else on the page.
* **It plays forwards then backwards.** A trajectory does not loop, so cutting
  frame 99 straight back to frame 0 reads as a glitch rather than as data. The
  ping-pong is stated in the caption so nobody counts it as twice the sampling.

Rendering happens here, on the droplet, rather than in the bundle, so that
archives already uploaded gain the clip on their next analysis and the bundle's
locked environment does not have to grow an encoder. It is bounded the same way
`dynamics.ligand_rmsd_to_reference` is bounded, by an atom x frame budget, and
it declines rather than risking the box when ffmpeg is missing or the budget is
past.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any, Optional

import numpy as np

from .dynamics import BUDGET_ATOM_FRAMES, _trajectory_ligand

WIDTH, HEIGHT = 720, 540       # even on both axes: h264 chroma subsampling needs it
DPI = 100
FPS = 15
CRF = "30"                     # flat background compresses hard; this lands near 400 kB
POCKET_CUTOFF_NM = 0.5         # 5 A from any ligand heavy atom, measured on frame 0
BOND_CUTOFF_NM = 0.18          # heavy-atom bond, generous enough for C-S at 1.8 A

BACKGROUND = "#243044"         # --panel
PROTEIN = "#7f8fa6"            # --grey, darkened: this is context, not the subject
POCKET = "#e8edf5"             # --text
LIGAND = "#5de1e6"             # --phos, the predicted state
CAPTION = "#9fb0c7"            # --muted


def render(topology: str | Path, trajectory: str | Path, mp4: str | Path,
           poster: str | Path, ligand_resname: Optional[str] = None,
           times_ps: Optional[list[float]] = None) -> dict[str, Any]:
    """Write the clip and its poster. Returns what the panel needs, or an error."""
    topology, trajectory = Path(topology), Path(trajectory)
    mp4, poster = Path(mp4), Path(poster)

    if not trajectory.exists() or not topology.exists():
        return {"error": "No trajectory in the archive, so there is nothing to play."}
    if shutil.which("ffmpeg") is None:
        return {"error": "ffmpeg is not installed on this server, so the clip was not encoded."}

    import mdtraj as md

    traj = md.load(str(trajectory), top=md.load_topology(str(topology)))
    if traj.n_atoms * traj.n_frames > BUDGET_ATOM_FRAMES:
        return {"error": f"That trajectory is {traj.n_atoms * traj.n_frames:,} atom-frames, past "
                         f"this server's {BUDGET_ATOM_FRAMES:,} budget, so it was not rendered."}
    if traj.n_frames < 2:
        return {"error": "One frame is not a trajectory: there is nothing to animate."}

    top = traj.topology
    ligand = _trajectory_ligand(top, ligand_resname)
    protein_ca = np.array([a.index for a in top.atoms
                           if a.name == "CA" and a.residue.name.upper() != "LIG"], dtype=int)
    if protein_ca.size < 4:
        return {"error": "Fewer than four Ca atoms in the trajectory: nothing to draw a chain with."}

    traj.superpose(traj, frame=0, atom_indices=protein_ca)
    xyz = traj.xyz * 10.0                                    # nm to Angstrom, for the caption's sake

    chains = _ca_chains(top)
    pocket = _pocket_residues(xyz[0], protein_ca, ligand, top)
    bonds = _ligand_bonds(xyz[0], ligand)
    rotation = _camera(xyz[0][protein_ca])
    limits = _limits(xyz @ rotation.T)

    frames = _draw(xyz @ rotation.T, chains, pocket, ligand, bonds, limits, times_ps, top)
    _encode(frames, mp4)
    _poster(frames[0], poster)

    return {"mp4": mp4.name, "poster": poster.name, "frames": int(traj.n_frames),
            "seconds": round(len(frames) * 2 / FPS, 1),
            "pocket_residues": len(pocket)}


def _ca_chains(top) -> list[list[int]]:
    """Ca indices per chain, in residue order, so the ribbon never jumps a break."""
    from .superpose import STANDARD_RESIDUES

    chains: list[list[int]] = []
    for chain in top.chains:
        run: list[int] = []
        previous: Optional[int] = None
        for residue in chain.residues:
            if residue.name.upper() not in STANDARD_RESIDUES:
                continue
            ca = next((a.index for a in residue.atoms if a.name == "CA"), None)
            if ca is None:
                continue
            # A gap in the numbering is a gap in the chain. Splining across it
            # draws a straight bar through the middle of the protein, which
            # looks like a helix nobody modelled.
            if previous is not None and residue.resSeq - previous > 1 and run:
                chains.append(run)
                run = []
            run.append(ca)
            previous = residue.resSeq
        if len(run) > 1:
            chains.append(run)
    return [c for c in chains if len(c) > 1]


def _pocket_residues(frame, protein_ca, ligand, top) -> set[int]:
    """Ca indices of residues with any heavy atom within 5 A of the ligand."""
    if ligand.size == 0:
        return set()
    lig_xyz = frame[ligand]
    pocket: set[int] = set()
    for ca in protein_ca:
        residue = top.atom(int(ca)).residue
        heavy = [a.index for a in residue.atoms if a.element.symbol != "H"]
        if not heavy:
            continue
        d = np.linalg.norm(frame[heavy][:, None, :] - lig_xyz[None, :, :], axis=-1)
        if d.min() <= POCKET_CUTOFF_NM * 10.0:
            pocket.add(int(ca))
    return pocket


def _ligand_bonds(frame, ligand) -> list[tuple[int, int]]:
    """Heavy-atom bonds by distance.

    The topology comes from a PDB written by OpenMM, which carries no CONECT
    records for a ligand merged in from an OpenFF topology, so MDTraj knows of
    no bonds at all inside it. Distance is enough: nothing in a relaxed complex
    puts two non-bonded heavy atoms inside 1.8 A.
    """
    if ligand.size < 2:
        return []
    pos = frame[ligand]
    d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
    i, j = np.where((d > 0.1) & (d <= BOND_CUTOFF_NM * 10.0))
    return [(int(ligand[a]), int(ligand[b])) for a, b in zip(i, j) if a < b]


def _camera(ca_frame) -> np.ndarray:
    """A rotation that puts the protein's widest two axes on screen.

    The principal axes of the Ca cloud, largest variance first. Looking down the
    smallest one shows the most protein and, because it is derived from the
    structure rather than fixed, does not leave some folds edge-on.
    """
    centred = ca_frame - ca_frame.mean(axis=0)
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    rotation = vt
    if np.linalg.det(rotation) < 0:          # keep it a rotation, not a reflection
        rotation[2] *= -1
    return rotation


def _limits(xyz) -> tuple[float, float, float, float]:
    """A window sized on what is on screen across every frame.

    Only the two axes the camera can see. Sizing on the largest of all three
    would let the depth axis, which nobody is looking down, decide the zoom and
    leave the protein small in the middle of an empty panel.
    """
    lo = xyz.reshape(-1, 3).min(axis=0)
    hi = xyz.reshape(-1, 3).max(axis=0)
    centre = (lo + hi) / 2
    span_x, span_y = float(hi[0] - lo[0]), float(hi[1] - lo[1])
    half_y = max(span_y, span_x * HEIGHT / WIDTH) / 2 * 1.04
    return float(centre[0]), float(centre[1]), half_y * WIDTH / HEIGHT, half_y


def _spline(points: np.ndarray) -> np.ndarray:
    """Smooth a Ca trace into something that reads as a ribbon rather than a zigzag."""
    from scipy.interpolate import splev, splprep

    if len(points) < 4:
        return points
    try:
        tck, _ = splprep(points.T, s=0, k=3)
    except (TypeError, ValueError):
        return points
    u = np.linspace(0, 1, len(points) * 6)
    return np.array(splev(u, tck)).T


def _draw(xyz, chains, pocket, ligand, bonds, limits, times_ps, top) -> list[np.ndarray]:
    """Render every frame to an RGB array. Matplotlib, because it is already here."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    cx, cy, half_x, half_y = limits
    fig = plt.figure(figsize=(WIDTH / DPI, HEIGHT / DPI), dpi=DPI, facecolor=BACKGROUND)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor(BACKGROUND)
    ax.set_xlim(cx - half_x, cx + half_x)
    ax.set_ylim(cy - half_y, cy + half_y)
    ax.set_aspect("equal")
    ax.set_axis_off()

    # Everything is drawn in the camera's xy plane and sorted by z, which is a
    # painter's algorithm: near things drawn last. Real depth sorting per
    # segment would cost more than it shows at this size.
    out: list[np.ndarray] = []
    for f in range(xyz.shape[0]):
        frame = xyz[f]
        for artist in list(ax.collections) + list(ax.texts):
            artist.remove()

        segments, colours, widths = [], [], []
        for chain in chains:
            path = _spline(frame[chain][:, :2])
            depth = frame[chain][:, 2]
            seg = np.stack([path[:-1], path[1:]], axis=1)
            near_pocket = np.array([c in pocket for c in chain], dtype=float)
            # Six spline points per residue, so a residue's colour spans six segments.
            spread = np.interp(np.linspace(0, len(chain) - 1, len(path)),
                               np.arange(len(chain)), near_pocket)[:-1]
            segments.extend(seg)
            # The site is drawn brighter and thicker than the rest of the fold.
            # Colour alone did not separate them: near-white against the muted
            # grey is a smaller step on screen than it looks in the palette.
            colours.extend([POCKET if v > 0.5 else PROTEIN for v in spread])
            widths.extend([2.9 if v > 0.5 else 1.7 for v in spread])
        ax.add_collection(LineCollection(segments, colors=colours, linewidths=widths,
                                         capstyle="round", zorder=2, alpha=0.95))

        if bonds:
            lig_seg = [[frame[a][:2], frame[b][:2]] for a, b in bonds]
            ax.add_collection(LineCollection(lig_seg, colors=LIGAND, linewidths=5.0,
                                             capstyle="round", zorder=4))
        if ligand.size:
            ax.scatter(frame[ligand][:, 0], frame[ligand][:, 1], s=26, c=LIGAND,
                       edgecolors="none", zorder=5)

        if times_ps and f < len(times_ps):
            ax.text(0.02, 0.03, f"{times_ps[f]:.0f} ps", transform=ax.transAxes,
                    color=CAPTION, fontsize=11, family="monospace", zorder=6)

        fig.canvas.draw()
        buf = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
        out.append(buf)
    plt.close(fig)
    return out


def _encode(frames: list[np.ndarray], dest: Path) -> None:
    """Forwards then backwards, so the loop has no cut in it."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    order = list(range(len(frames))) + list(range(len(frames) - 2, 0, -1))
    height, width, _ = frames[0].shape
    cmd = ["ffmpeg", "-y", "-loglevel", "error",
           "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{width}x{height}",
           "-r", str(FPS), "-i", "-",
           "-c:v", "libx264", "-preset", "slow", "-crf", CRF,
           "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(dest)]
    proc = subprocess.run(cmd, input=b"".join(frames[i].tobytes() for i in order),
                          capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(f"ffmpeg failed: {proc.stderr.decode('utf-8', 'replace')[:400]}")


def _poster(frame: np.ndarray, dest: Path) -> None:
    """The still the panel shows before anything is downloaded."""
    from PIL import Image

    dest.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(frame).save(dest, format="WEBP", quality=82, method=5)
