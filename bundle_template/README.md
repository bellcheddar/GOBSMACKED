# GOBSMACKED run bundle

Everything in this directory was written for one campaign by
[gobsmacked.mdeller.com](https://gobsmacked.mdeller.com). It runs on your
machine, uses your GPU, and talks to nothing.

```bash
pixi install
pixi run gobsmacked          # writes results/results.tar.gz
```

Then upload `results/results.tar.gz` on the Analyze tab of the site that
generated this bundle.

## What runs

| Stage | What it does | Typical cost |
|---|---|---|
| `fold` | ESMFold, **skipped** when `model_apo.pdb` is present (usual case) | 0 to 3 min |
| `prep` | PDBFixer on the receptor, RDKit conformer for the ligand | under a minute |
| `dock` | PandaDock, in the mode `campaign.yaml` asks for | 3 to 10 min |
| `md` | OpenMM: minimise, equilibrate with restraint release, production | 10 to 20 min per ns |
| `summarise` | MDTraj: RMSDs, RMSF, pocket volume, contact matrix, then packs the archive | 1 to 3 min |

Each stage writes `work/<stage>.done` when it finishes, so a rerun picks up
where it stopped. To redo one stage and everything after it:

```bash
pixi run gobsmacked --stage dock
pixi run gobsmacked --list        # what would run
```

## What it writes

```
results/
  manifest.json          schema version, engine versions, timings, warnings
  campaign.yaml          echoed back, exactly as run
  model_apo.pdb          the prepared starting structure
  plddt.json             per-residue pLDDT, when this bundle folded the sequence
  poses/poses.sdf        every docked pose, with scores in SD tags
  poses/scores.csv       pose_id, score, GNN affinity, rank
  complex_pose1.pdb      top pose merged with the receptor
  complex_min.pdb        after minimisation
  complex_md_final.pdb   the last MD frame
  traj/traj.dcd          production trajectory, solute only, pocket-aligned
  traj/topology.pdb
  traj/summary.json      per-frame RMSDs, RMSF, pocket volume, contact matrix
  logs/run.log
```

The archive is validated on the way back in, and an incomplete one is refused
with the missing file named, so a failed stage never turns into a puzzling
analysis.

## If something fails

* **No GPU found.** MD falls back to the CPU at roughly a tenth of the
  throughput. Lower `production_ps` in `campaign.yaml` before running, rather
  than waiting it out.
* **PandaDock cannot fetch its GNN checkpoint.** The run continues in `dock`
  mode, ranked by the empirical function, and says so in the warnings.
* **ESMFold runs out of memory.** `fold.py` already chunks by sequence length;
  a card under 12 GB may still need a shorter construct. Trimming to the domain
  is usually the right answer anyway.
* **The ligand has no MMFF parameters.** UFF is used for the starting conformer
  and the run continues. Docking and MD are unaffected: they use their own force
  fields.

## Licences

GOBSMACKED is MIT. PandaDock (MIT), OpenMM (MIT), PDBFixer (MIT),
openmmforcefields and OpenFF (MIT), MDTraj (LGPL-2.1), RDKit (BSD-3) and ESM
(MIT) are installed by pixi from their own channels. PLIP is GPL-2.0 and is
deliberately **not** part of this bundle: it runs server-side only.
