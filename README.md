# 🔬 GOBSMACKED

> **Fold, dock, relax and annotate a protein-ligand complex, then check it against the experimental structure.**

[![live](https://img.shields.io/badge/live-gobsmacked.mdeller.com-00d084?logo=icloud&logoColor=white)](https://gobsmacked.mdeller.com) ![python](https://img.shields.io/badge/python-3.12.3-3776AB?logo=python&logoColor=white) ![flask](https://img.shields.io/badge/flask-3.1.3-000000?logo=flask&logoColor=white) ![gunicorn](https://img.shields.io/badge/gunicorn-26.2.0-499848?logo=gunicorn&logoColor=white) ![nginx](https://img.shields.io/badge/nginx-1.24-009639?logo=nginx&logoColor=white) ![sqlite](https://img.shields.io/badge/sqlite-3-003B57?logo=sqlite&logoColor=white) ![rdkit](https://img.shields.io/badge/rdkit-2026.3.6-3838AB) ![biotite](https://img.shields.io/badge/biotite-1.6.0-467FF7) ![gemmi](https://img.shields.io/badge/gemmi-0.7.5-467FF7) ![mdtraj](https://img.shields.io/badge/mdtraj-1.11.1-467FF7) ![plip](https://img.shields.io/badge/PLIP-3.0.1-9b51e0) ![pandamap](https://img.shields.io/badge/PandaMap-4.3.0-9b51e0) ![pandadock](https://img.shields.io/badge/PandaDock-4.1.1-9b51e0) ![openmm](https://img.shields.io/badge/OpenMM-8.6.0-00897B) ![boltz](https://img.shields.io/badge/Boltz--2-2.2.1-00897B) ![esmfold](https://img.shields.io/badge/ESMFold-v1-00897B) ![tmtools](https://img.shields.io/badge/TM--align-0.3.0-00897B) ![molstar](https://img.shields.io/badge/Mol*-5.11.0-467FF7) ![plotly](https://img.shields.io/badge/Plotly.js-2.35.2-3F4F75?logo=plotly&logoColor=white) ![tests](https://img.shields.io/badge/pytest-160%20passing-00d084) ![data](https://img.shields.io/badge/data-RCSB%20%C2%B7%20UniProt%20%C2%B7%20AlphaFold%20DB%20%C2%B7%20KLIFS%20%C2%B7%20GPCRdb%20%C2%B7%20InterPro-467FF7) ![licence](https://img.shields.io/badge/licence-MIT-lightgrey) ![author](https://img.shields.io/badge/author-Marc%20C.%20Deller%2C%20D.Phil.-1C244B)

<table>
<tr>
<td>🌐 <b>App</b></td><td><a href="https://gobsmacked.mdeller.com" target="_blank" rel="noopener noreferrer">gobsmacked.mdeller.com</a></td>
<td>✉️ <b>Contact</b></td><td><a href="mailto:marc@marcdeller.com">marc@marcdeller.com</a></td>
<td>🐙 <b>GitHub</b></td><td><a href="https://github.com/bellcheddar/GOBSMACKED" target="_blank" rel="noopener noreferrer">bellcheddar/GOBSMACKED</a></td>
</tr>
</table>

---

![The scorecard for EGFR plus erlotinib judged against crystal structure 1M17: a grade B dial at 83.5, six graded gauges each with a sentence explaining what to do about it, the relaxed complex in Mol*, PandaMap's 2D interaction diagram, and the kinase switch list reporting DFG-in, alphaC-out and a Type I binding mode matching the crystal](docs/screenshots/scorecard.png)

**GOBSMACKED** (Ground-truth Overlay for Binding Sites, Modes And Complex Kinetics/Dynamics) takes a protein and a ligand, folds and docks and relaxes them, and then does the thing most docking pipelines skip: it goes and finds the crystal structure, superposes on the binding pocket, and tells you how close you got and why.

**Why it matters:** a docking score is a ranking, not a measurement, and a pretty predicted complex looks exactly the same whether it is right or wrong. GOBSMACKED answers three separate questions about one prediction: how close the pose lands to the crystal (PIER REVIEW), whether molecular dynamics recovers the induced fit that an apo-like predicted pocket is missing (HOLOGRAM), and whether the binding mode the prediction implies is the binding mode the crystal shows (GATEKEEPER). It is useful for: anyone validating a docking protocol before trusting it on a target with no structure, anyone asking whether ESMFold plus docking is good enough for a particular pocket, and anyone who wants the answer as a graded scorecard rather than as a folder of PDB files.

---

## 🧭 How it works

The heavy compute does not run on the server. ESMFold, the PandaDock GNN and OpenMM need a GPU and several gigabytes; the host is a shared CPU droplet with 3.8 GB. So the work is split in two, with an archive passing between them.

![The pipeline: Prepare on the droplet fetches, annotates, picks the pocket and the reference, and emits run_bundle.tar.gz; Run on your GPU folds, preps, docks, runs MD and summarises into results.tar.gz; Analyze on the droplet superposes, runs the interaction analysis, grades, classifies the binding mode and draws the dynamics](docs/screenshots/pipeline.png)

| Stage | Where | What happens |
|---|---|---|
| **Prepare** | droplet, CPU | Resolve the input, fetch the best available structure, annotate the family, pick the pocket, choose a reference crystal, emit `run_bundle.tar.gz` |
| **Run** | your machine, GPU | `pixi run gobsmacked`: fold (if needed), prep, dock, minimise and run MD, score the affinity before and after, summarise, emit `results.tar.gz` |
| **Analyze** | droplet, CPU | Validate the archive, superpose on the pocket, run PLIP and PandaMap, grade, classify the binding mode, draw the trajectory |

Nothing in the bundle contacts the server. The campaign file goes in, the results archive comes back, and both are validated against a schema so a failed stage never turns into a puzzling analysis.

## 🚀 Quick start

Open [gobsmacked.mdeller.com](https://gobsmacked.mdeller.com), paste a UniProt accession and a SMILES, pick a pocket, and download the bundle. Then, on a machine with a GPU:

```bash
curl -fL "https://gobsmacked.mdeller.com/runs/<job-id>/bundle" | tar xz \
  && cd run_bundle_<job-id> && ./run.sh
```

Prepare prints that line with the job filled in and a copy button. One step, and the only
prerequisites are `curl` and `bash`: the bundle carries its own `pixi.lock`, and `run.sh`
installs pixi if the machine has not got it, builds the environment from the lock (no solve, so
the versions are the ones the bundle was tested with) and runs all six stages. The environment
lives in `.pixi/` inside the bundle directory, so nothing is installed system-wide.

Upload `results.tar.gz` on the Analyze tab when it finishes. It is written at the top level of the
bundle directory, next to `run.sh`, rather than inside `results/`: the archive is the one file the
run exists to produce and nobody should have to go looking for it.

To run the web application locally:

```bash
git clone https://github.com/bellcheddar/GOBSMACKED.git
cd GOBSMACKED
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python -r requirements.txt
make serve                   # http://127.0.0.1:8009
```

## 🔧 Prepare

Four panels, each unlocking the next.

| Panel | Accepts | What it does |
|---|---|---|
| **1. Protein** | UniProt accession, PDB ID, raw sequence, or a dropped `.pdb` / `.cif`, plus an optional domain range | Resolves the canonical sequence and picks a starting structure by a fixed, reported priority: your upload > a named PDB entry > AlphaFold DB > ESM Atlas > fold it in the bundle. **From Pfam** fills the range from the domain the pocket sits in |
| **2. Annotation** | (automatic) | InterPro supplies Pfam domains, which route the target: `PF00069`/`PF07714` to KLIFS, `PF00001`/`PF00002`/`PF00003`/`PF10324` to GPCRdb, anything else to UniProt features alone |
| **3. Ligand and pocket** | SMILES, plus residues clicked in Mol* or on the sequence track | RDKit validates and draws the ligand; the docking box is the selection's extent plus 8 Å, floored at 18 Å per side |
| **4. Reference** | (automatic, or a typed PDB ID) | RCSB entries mapped to the accession, ranked by Morgan Tanimoto against your ligand, with resolution and a 2D depiction |

![The Prepare tab: the four input panels on the left and right, an empty Mol* scope in the centre waiting for a structure, and the seven-cell stage strip across the top showing what the bundle will run](docs/screenshots/prepare.png)

**Trim to the domain.** An AlphaFold model is of the whole precursor. EGFR's is 1,210 residues, of which the kinase domain is 253: solvating the other 957 costs an order of magnitude in MD time, adds two disordered tails that wander through the box, and tells you nothing about the pocket.

**The reference's site is renumbered onto your model.** 1M17 numbers EGFR from the mature protein and UniProt from the precursor, 24 apart. Pasting the crystal's residue list straight into the pocket picker would select real residues that are the wrong ones, which is worse than an error, so both structures are aligned by sequence first.

**The docking box is sized from the ligand, not from the residues around it.** Erlotinib's own extent plus padding is 32 x 21 x 23 A; the shell of its 51 contacting residues gives 48 x 42 x 42, a search volume five times larger that samples the true site less densely for no benefit. The centre still comes from the model, because the crystal is in its own coordinate frame and only the extent transfers.

**Reference selection prefers the same ligand over a sharper crystal.** EGFR and erlotinib make the case: 1M17 holds erlotinib itself at 2.6 Å, while the best sub-2.5 Å entry holds gefitinib at Tanimoto 0.41. Judging an erlotinib pose against a gefitinib crystal because the crystal is 0.9 Å sharper would measure the wrong thing.

**Visibility** is chosen here and defaults to public. A private run is issued a 32-character owner key, shown once, stored only as its sha256, and carried inside the bundle so uploading results needs no typing. Private runs are **absent** from the Runs table without the key rather than greyed out: a greyed row would leak that the run exists, and how many there are.

## ⚗️ Run

Six stages, each idempotent and resumable from a `.done` marker.

| Stage | Tool | Notes |
|---|---|---|
| `fold` | ESMFold | Skipped when the bundle carries a model, which is the usual case. Chunk size scales with sequence length; pocket residues below pLDDT 70 raise a warning that reaches the scorecard |
| `prep` | PDBFixer, RDKit | Missing atoms, hydrogens at the campaign pH, waters and heteroatoms removed. Terminal missing residues are deliberately not built: they are absent from the construct, not from the model |
| `dock` | PandaDock | `hybrid` (search plus SE(3) GNN rescoring), `flex` (induced fit) or `dock` (empirical only). Falls back from `hybrid` to `dock` when the GNN checkpoint cannot be fetched, and says so |
| `md` | OpenMM, OpenFF | Amber14 plus OpenFF Sage, TIP3P with 0.15 M NaCl and 10 Å padding, restraints released over the equilibration, 2 fs with hydrogen mass repartitioning. The DCD holds the solute only |
| `affinity` | Boltz-2 | Optional, on by default. The docked pose and frames sampled from the last fifth of the trajectory, each scored by Boltz-2's affinity head with the structure module bypassed. The MSA is computed once per target and cached, so only the first pose queries the server |
| `summarise` | MDTraj | Per-frame ligand and backbone RMSD, per-residue RMSF, pocket volume by voxel counting, a residue-by-frame contact matrix, then packs the archive |

```bash
./run.sh --list          # what would run, and what is already done
./run.sh --stage dock    # rerun from dock onward
```

Target on one consumer GPU: under 30 minutes for a 300-residue domain with the default 1 ns
production, and roughly an hour with the affinity stage on, which is the slowest of the six once
the MD is short. The runner prints a per-stage plan with wall-clock estimates and a finishing time
before it starts, then a progress bar and a running ETA for each stage as it goes.

**Warnings are triaged rather than dumped.** Every third-party warning the run emits was read
once, and each is either fixed at the source, suppressed by an exact-match filter with the reason
recorded next to it, or left visible because it is telling you something. What reaches the terminal
is what a reader should act on.

## 📊 Analyze

Upload the archive and the whole pipeline runs inside the request: ingest, superpose, interactions, modes, dynamics, scorecard. Seconds, not minutes.

**Superposition is on pocket Cα atoms, never on the whole chain.** A model can be excellent at the binding site and 6 Å out at a disordered terminus; superposing whole chains spreads that error into the pocket and inflates the ligand RMSD, which is the one number this app exists to report honestly. The whole-chain TM-score is reported alongside, as context.

**Residue numbering is never assumed to match.** The reference chain and the model are aligned by sequence first and every measurement walks that mapping. 1M17 numbers EGFR from the mature protein, 24 lower than UniProt: comparing raw numbers would make every contact look lost and the interaction overlap read as zero.

### 🎯 The scorecard

| Metric | A | B | C | D | F | Weight |
|---|---|---|---|---|---|---|
| Ligand RMSD, best of pose 1 and MD-final | ≤ 1.0 Å | ≤ 2.0 | ≤ 3.0 | ≤ 4.0 | > 4.0 | 30 |
| PLIP interaction Jaccard, best of pose 1 and MD-final | ≥ 0.75 | ≥ 0.55 | ≥ 0.40 | ≥ 0.25 | < 0.25 | 20 |
| Pocket Cα RMSD, MD-final | ≤ 0.8 | ≤ 1.2 | ≤ 1.8 | ≤ 2.5 | > 2.5 | 15 |
| χ1 agreement, MD-final (within 40°) | ≥ 0.85 | ≥ 0.70 | ≥ 0.55 | ≥ 0.40 | < 0.40 | 10 |
| MD stability: ligand drift, last window minus first | ≤ 0.5 Å | ≤ 1.0 | ≤ 1.5 | ≤ 2.5 | > 2.5 | 10 |
| Rescue: pocket Cα RMSD before MD minus after | ≥ +0.5 Å | ≥ +0.2 | ≥ 0 | ≥ −0.3 | < −0.3 | 10 |
| Pose validity: clashes, bond lengths, chirality, inside the box | pass | | | | fail | 5 |

The composite **GOBSMACK score** is the weighted mean of those grades. Three rules keep it honest:

- **A metric that could not be measured drops out** and the remaining weights are renormalised, which the card states rather than hiding in the arithmetic.
- **A validity failure caps the composite at 40.** A pose with a 1.8 Å clash is not a B whatever else it scored.
- **A run with no reference gets no composite at all.** There is nothing to verify it against, and scoring it on the two metrics that survive would be a grade for something nobody checked.

Every gauge carries one plain sentence saying what the number means and what to do about it, because "χ1 agreement 0.42" is not actionable and "most pocket side chains are in the wrong rotamer, which is what an apo-like predicted pocket looks like: try flex docking plus a longer equilibration" is.

### 🧬 The overlay

Model in grey, MD-final in phosphor, the crystal in amber, all superposed on the pocket, with the ten most displaced side chains and their χ1 angles listed underneath.

![The three-way overlay for the beta-2 adrenergic receptor: the relaxed complex in cyan and crystal structure 2RH1 in amber, superposed on 53 pocket Ca atoms at TM 0.999, with carazolol drawn in pink, above a table of the most displaced pocket side chains with their chi1 angles in the model and in the crystal](docs/screenshots/overlay.png)

### 🧲 Affinity, before and after MD

Boltz-2's affinity head scores the docked pose and the relaxed complex, so the panel answers a
question the rest of the scorecard cannot: did relaxing the pose change what the model thinks of
it. Each is reported as pIC50, the raw `affinity_pred_value` (log10 of IC50 in micromolar, lower
is stronger) and the binder probability, with the change between them. The post-MD figure is a
mean over several frames from the last fifth of the trajectory, with its standard deviation, so a
prediction that swings by a log unit across the sampled window says so rather than presenting one
frame as the answer.

**It is deliberately outside the composite score.** Every graded metric on the scorecard has a
crystal structure to be right or wrong about. A predicted affinity has none, and folding an
unverifiable number into a score whose whole argument is verification would break that argument.
It sits beside the grade, not inside it.

**A missing affinity is a missing panel, not a failed run.** The stage declines rather than
raising: no ligand it can parameterise, no MSA, the affinity head writing nothing for a pose. The
reason is carried into the archive and shown on the card, and the other five stages are untouched.

### 🎲 Every pose, not just the top one

Docking returns ten poses and the scorecard grades one of them, which hides the case that matters
most: the run that found the right answer and ranked it fourth. The overlay panel draws all ten at
once and tabulates, per pose, the docking score, the in-place RMSD to the top pose, the centroid
separation, a best-fit shape RMSD, and the closest heavy-atom approach to the receptor.

The distinction between the two RMSDs is the useful part. In-place RMSD asks how far this pose
sits from that one; best-fit shape RMSD asks whether it is the same conformer put somewhere else.
A pair that is 6 Å apart in place and 0.4 Å after superposition is one pose docked into two sites,
which is a search problem. A pair that is 6 Å apart both ways is two genuinely different bound
conformations, which is not.

### 🔑 The binding mode

Two families get a real answer and everything else gets an honest one.

**Kinases** are labelled from the KLIFS 85-residue pocket, mapped onto the target sequence region by region. DFG-in / out / inter follows the Modi and Dunbrack distance criteria; αC-in / out follows the β3-Lys to αC-Glu salt bridge; the Type I / I½ / II / allosteric label follows which subpockets the ligand occupies.

**GPCRs** are labelled from GPCRdb generic numbering: orthosteric, vestibule, intracellular or lipid-facing from the contacts, and active-like or inactive-like from the TM3-TM6 distance, with the NPxxY RMSD, the toggle switch χ1, the PIF motif and the sodium site reported alongside.

Both classifiers run on the prediction and on the crystal with the same code, so a difference in the label is a difference in the structure and not a difference in method.

### 🧪 Thresholds that were measured, not quoted

Two numbers in this repository were placed by measuring structures with this code rather than by taking a value from a paper that used a different atom pair:

- **The TM3-TM6 activation cut is 13 Å**, from six structures: inactive rhodopsin 1GZM 8.7, A2A 3EML 9.7, β2AR 2RH1 11.2; active metarhodopsin 3PQR 14.7, A2A 5G53 18.5, β2AR 3SN6 19.0. (Quoted "ionic lock" distances of 3 to 4 Å are guanidinium-to-carboxylate, not Cα-to-Cα, and are not comparable.)
- **A subpocket counts as occupied at two contacts, not one.** Erlotinib in 1M17 grazes exactly one αC residue at 4 Å, and a one-contact rule labels a textbook Type I inhibitor as Type I½.

## 🧱 Stack

| Component | Role | Licence | Where it runs |
|---|---|---|---|
| Flask, gunicorn, nginx, SQLite | Web application and store | BSD / MIT / public domain | droplet |
| biotite, tmtools, gemmi, NumPy, SciPy | Alignment, superposition, TM-score | BSD / MIT / MPL | droplet |
| RDKit | SMILES, depiction, fingerprints, symmetry-aware RMSD | BSD-3 | droplet and bundle |
| PLIP | Interaction fingerprints | GPL-2.0 | droplet only, as a subprocess |
| PandaMap | 2D interaction maps, empirical ΔG | MIT | droplet |
| PandaDock | Docking (hybrid search plus SE(3) GNN rescoring) | MIT | bundle |
| Boltz-2 | Affinity head, on the docked pose and on MD frames | MIT | bundle |
| ESMFold, OpenFold | Folding when no model exists | MIT / Apache-2.0 | bundle |
| OpenMM, PDBFixer, openmmforcefields | Preparation, minimisation, MD | MIT / LGPL | bundle |
| MDTraj | Trajectory analysis | LGPL-2.1 | bundle and droplet |
| Mol*, Plotly.js | 3D views and traces | MIT | browser |

**PLIP is GPL-2.0 and is run as a subprocess, never imported.** What crosses the boundary is an XML file, so no GPL code is linked into this MIT-licensed application, and PLIP is deliberately absent from the run bundle. Full attribution, with every DOI checked against Crossref, is in [THIRD_PARTY.md](THIRD_PARTY.md) and on the app's About page: both are generated from `software.yaml`, so they cannot drift.

## 📁 Repository layout

```
app/                     the Flask application (droplet only, no torch)
  routes/                prepare, analyze, runs, about
  services/              fetch, annotate, references, bundle, ingest, superpose,
                         interactions, scorecard, modes, dynamics, affinity,
                         poses, movie
  templates/  static/    Jinja2, the instrument-panel CSS, Mol* and Plotly
bundle_template/         copied verbatim into every run bundle
  run.sh                 one-step bootstrap: installs pixi, builds from the lock, runs
  pixi.lock              1,969 pinned package entries across three environments
                         (default, fold, affinity) and two platforms
  run.py                 the six-stage runner with resume markers
  gobsmacked_run/        fold, prep, dock, md, affinity, summarise, schema,
                         console, noise
design/                  the visual contract this app is built from
deploy/                  systemd units, nginx site, provision and deploy scripts
scripts/                 prune, DOI checker, THIRD_PARTY generator
tests/                   160 tests, plus two fixture archives built from crystals
software.yaml            the single source of truth for attribution
```

## 🔬 A worked example

One real run, end to end, on an M1 Max: EGFR plus erlotinib, prepared on the live site,
judged against 1M17.

![The results page for the real EGFR run: a grade D dial at 53.5, ligand RMSD 7.83 A graded F with the note that this usually means the wrong site rather than a scoring failure, MD drift 0.25 A graded A, the relaxed complex in Mol* with its PLIP interactions, and the kinase switch list reporting no hinge contact against the crystal's Type I](docs/screenshots/real-run.png)

| Stage | Wall clock | What happened |
|---|---|---|
| fetch | seconds | AlphaFold model `AF-P00533-F1`, mean pLDDT 75.9 |
| annotate | seconds | Pfam `PF07714` 714-966, KLIFS gatekeeper Thr790 |
| prepare | seconds | 1M17's site renumbered onto the model (+24), box sized on erlotinib at 32 x 21 x 23 A |
| prep | 3 s | trimmed 1,210 residues to the 253-residue kinase domain, all 51 pocket residues kept |
| dock | 937 s | PandaDock empirical search, 10 poses, best score −15.5 kcal/mol |
| md | 3,224 s | 58,266 atoms, OpenCL single precision, ~20 ns/day, 100 ps equilibration plus 500 ps production |
| summarise | 106 s | 100 frames |
| analyse | 28 s | on the droplet |

**GOBSMACK 53.5, grade D.** The interesting part is what that decomposes into:

| Metric | Value | Grade |
|---|---|---|
| Ligand RMSD | 7.83 Å | F |
| PLIP overlap | 0.33 | D |
| Pocket Cα RMSD | 1.66 Å | C |
| χ1 agreement | 0.69 | C |
| Drift, last window | 0.25 Å | A |
| MD rescue | −0.75 Å | F |

The pose is stable and physically valid: it does not move over 500 ps, it has no clash, its
stereocentres are intact and it sits inside the box. It is also 8 Å from where erlotinib
actually binds, it never touches the hinge, and the mode classifier calls it allosteric
against the crystal's Type I. A docking score of −15.5 kcal/mol says nothing about any of
that, which is the entire argument for scoring a prediction against the structure rather
than against itself.

Two caveats stated rather than buried: this run used the empirical scorer alone, because the
SE(3) GNN rescoring needs PyTorch and this machine had 12 GB of disk left, and 500 ps is a
short production run. Both are reasons to expect a worse pose, and neither changes what the
scorecard is reporting.

## 🧪 An experiment: five starting structures, one campaign

The worked example above scores a D and says the pose is 8 Å from where erlotinib
actually binds. The obvious reading is that the starting model was not good enough. That
is a testable claim, so it was tested.

Five runs of the same campaign, differing in one variable only: the structure docking
starts from. Same ligand, same pocket, same box (31.7 x 21 x 23 Å), same seed, same
exhaustiveness, same MD protocol, all judged against the same crystal. The first is a
control that should be unbeatable, because it docks erlotinib back into the very crystal
it came from.

| # | Starting structure | Pocket Cα RMSD to 1M17 | Top-ranked pose | Best pose (its rank) | GOBSMACK |
|---|---|---|---|---|---|
| 1 | **1M17 crystal, self-dock** | **0.00 Å** | **8.77 Å** | 1.51 Å (rank 8) | D 52.8 |
| 2 | 4HJO crystal, cross-dock | 1.53 Å | 4.29 Å | 2.29 Å (rank 5) | D 50.5 |
| 3 | AlphaFold DB `AF-P00533-F1` | 2.05 Å | 7.73 Å | 4.81 Å (rank 8) | D 55.0 |
| 4 | ESMFold, ESM Atlas | 1.83 Å | 8.09 Å | 5.44 Å (rank 9) | D 49.0 |
| 5 | Boltz-2 co-folded, ligand stripped | 1.28 Å | 1.40 Å | 1.40 Å (rank 1) | B 79.0 |

**The control settles it. A perfect receptor produced one of the worst top-ranked
poses.** Docking erlotinib into its own crystal, with a pocket that is correct by
construction, put an 8.77 Å pose first. The crystal pose was found: it is in the same
list of ten, at 1.51 Å. It was placed eighth.

So the sampling is not the bottleneck and the receptor is not the bottleneck. Within
this pose set the limitation is **which pose gets ranked first**, which is a
well-documented weakness of fast empirical scoring functions in general and not a
property of any one implementation. They are built to be quick enough to search millions
of placements, and that speed is bought against exactly this discrimination.

Two secondary results, both of which change how the other numbers should be read:

- **Starting-model quality does not predict pose quality.** With the 0.00 Å point
  included there is no monotonic relationship down the first two columns: the best
  possible receptor sits third from the worst on top-ranked pose. Reporting "the model
  was too rough" would have been a comfortable and wrong conclusion, and four of the five
  runs on their own would have supported it.
- **The predicted affinity does not separate right poses from wrong ones.** pIC50 ranged
  over 6.44 to 6.66 across poses spanning 1.4 to 8.8 Å, and the run with the *worst*
  top-ranked pose returned the *highest* affinity of the five. This is the measured
  argument for keeping affinity beside the composite score rather than inside it.

The one run that ranked correctly is worth stating carefully rather than celebrating.
Its receptor came from co-folding the protein *with* erlotinib and then removing the
ligand, so its pocket is already ligand-adapted in a way an apo or predicted-apo pocket
is not. It is one run, and a single B among four Ds is a hypothesis, not a result.

**This is the entire argument for the app.** Every one of these five runs produced a
stable, physically valid complex with a confident docking score, and four of them were
wrong by 4 to 9 Å. Nothing internal to a docking run distinguishes them. Only the
comparison against the experimental structure does.

## 🧫 Testing

```bash
make test                # 160 tests, about 90 seconds
make check-refs          # every DOI in software.yaml, checked against Crossref
make third-party         # regenerate THIRD_PARTY.md from software.yaml
python tests/fixtures/build_fixtures.py    # rebuild the two fixture archives
```

The fixtures are built from real crystals rather than from noise: 4HJO judged against 1M17, and 5D5A against 2RH1, with the ligand displaced to make a plausible docked pose and the bundle's own `summarise` stage computing the trajectory summary. So they exercise the last stage of the bundle as well as the first stage of the server.

End to end, EGFR plus erlotinib scores 83.5 (B) against 1M17 and labels Type I, DFG-in on both sides; β2AR plus carazolol scores 94.0 (A) against 2RH1 and labels orthosteric, inactive-like on both.

## 🌐 Deployment

```bash
cp .env.example .env         # DROPLET_SSH, SERVER_NAME, BIND_ADDR
bash deploy/deploy.sh        # rsync, reinstall dependencies, restart, verify the live page
```

First time only, on the droplet as root:

```bash
sudo SERVER_NAME=gobsmacked.mdeller.com bash /opt/gobsmacked/deploy/provision.sh
```

That installs the service user, the virtual environment, the systemd unit, the nightly prune timer, the nginx site and a Let's Encrypt certificate. The prune drops the results archive of public runs older than 90 days while keeping the scorecard, the report and the final structures, so the Runs row stays useful. Private runs are never pruned.

## ✅ To Do

Roadmap for GOBSMACKED, in dependency order. Suggestions welcome.

- [x] **Prepare, without references.** Input resolution across UniProt, RCSB, AlphaFold DB, ESM Atlas and uploads, with the priority reported rather than silently applied. Pfam routing, KLIFS and GPCRdb clients behind a 30-day cache, the Mol* pocket picker and the SVG sequence track
- [x] **Map the KLIFS pocket onto an arbitrary sequence.** The 85-character pocket string carries no residue numbers and only its regions are contiguous, so it is placed region by region, longest run first, with a bounded difflib pass for regions the alignment shifted. 83 of 85 positions map for EGFR, putting the gatekeeper on Thr790, the hinge on Met793 and DFG on Asp855
- [x] **The run bundle.** Five resumable stages in a self-contained pixi environment, with the GNN fallback, the platform choice and the solute-only trajectory
- [x] **Analyze, with verification.** Pocket-Cα superposition, symmetry-corrected ligand RMSD, χ1 agreement, PLIP fingerprints compared in reference numbering, PandaMap, the scorecard and the dynamics panels
- [x] **Both binding-mode classifiers, on day one.** Kinase and GPCR, each validated against structures whose labels are known: 1M17 Type I DFG-in, 1IEP Type II DFG-out, 2RH1 orthosteric inactive-like, 3SN6 active-like
- [x] **Runs and ownership.** Public and private runs, an owner key stored only as a hash, unguessable job IDs, and private runs absent from listings rather than redacted in them
- [x] **The About page.** A hand-drawn pipeline schematic, the grade thresholds, and a software table generated from the same file as THIRD_PARTY.md with every DOI checked against Crossref
- [x] **Deploy to gobsmacked.mdeller.com.** Provisioned, certificated and serving over HTTP/2, with the nightly prune timer armed. The static location deliberately sets no `access_log`: nginx.conf gives the droplet the `vhost` format that appends the requested host, and naming `combined` in the location would silently zero this app's visit count on the mdeller.com launcher
- [x] **Run the real loop once, end to end.** EGFR plus erlotinib, prepared on the live site, run on an M1 Max, uploaded and scored: see the worked example below. It found nine bugs the crystal-built fixtures could not, six of which failed silently
- [x] **PLIP interactions drawn in Mol\*, and the pocket as sticks.** Mol*'s viewer build exports no shape builder, so each interaction is loaded as a tiny structure of two-atom fragments joined by CONECT records: a run of them along PLIP's own endpoints reads as a dashed line, one file and one colour per interaction type. The lines are therefore the interactions the table lists rather than a second opinion computed by the viewer, which would quietly disagree with it
- [x] **An apo reference in the overlay.** Optional on Panel 4, fetched and stripped of its ligands at analysis time, and drawn in purple as a fourth toggle: the shape the pocket has with nothing bound, which is the shape a predicted model tends to resemble
- [x] **The trajectory played back beside the score.** Rendered on the droplet from the same trajectory the dynamics traces come from, so archives uploaded before it existed gain a clip on re-analysis: a Cα trace superposed on the protein, the pose in phosphor cyan, the residues lining its site in white, and forwards-then-backwards playback because a trajectory does not loop and a cut from the last frame to the first reads as a glitch rather than as data
- [x] **Affinity, before and after MD.** Built as stage 5 of six, on by default. Boltz-2's affinity head scores the docked pose and frames sampled from the last fifth of the trajectory, following the Boltzina pattern of feeding an existing pose straight to the affinity module with the structure module bypassed. Reported as pIC50, the raw log10(IC50/µM) and the binder probability for both, with the change between them and the spread across the sampled frames, and deliberately outside the composite score: every graded metric has a crystal to be right or wrong about, and a predicted affinity has none. It lives in its own pixi environment with `no-default-feature = true`, because boltz pins a torch that the docking environment must not inherit, and it never fails the run: a pose it cannot score becomes a missing panel with the reason attached
- [x] **Every pose in the overlay, with the numbers that separate them.** All ten drawn at once, and per pose the docking score, in-place RMSD to the top pose, centroid separation, best-fit shape RMSD and closest approach to the receptor. In-place against best-fit is the pair that matters: the same conformer in two sites is a search problem, two different conformers is not
- [x] **Make flex and hybrid docking actually run.** Four bugs, found only by using them rather than by reading them. `flex` treats `-o` as a filename prefix and writes to `<prefix>_results`, so every flex run appeared to produce nothing; `hybrid` accepts neither `--seed` nor `-e` and exits on either; the search radius was half the box's *longest* side, giving a sphere that reached well outside the box it was supposed to describe; and the GNN checkpoint had moved to a `v4` name and release URL, so `hybrid` silently fell back to the empirical scorer on every run
- [x] **Five starting structures, one campaign.** The experiment above: same ligand, pocket, box, seed and MD protocol, varying only the structure docking starts from, with a self-dock control that is correct by construction. It found that the top-ranked pose, not the sampling and not the receptor, is what limits the result, and that neither starting-model quality nor predicted affinity separates a right pose from a wrong one
- [ ] **Rank scoring functions on a fixed pose set.** The five runs left 50 poses whose distance to the crystal is already known, so the ranking question can be asked directly and cheaply: rescore the same poses with several independent scoring functions and ask which one puts a near-native pose first, per receptor type. No docking and no MD, because the search already found the answer every time
- [ ] **STEVEDORE: multi-ligand SAR series.** Score a congeneric series against one reference and correlate with ChEMBL affinity, which turns a single verification into a protocol assessment
- [ ] **DOCKYARD: ingest poses from other engines.** Boltz-2, Vina and DiffDock all produce poses this scorecard could grade, and the comparison is more interesting than any single engine's self-report
- [ ] **Cryptic pocket detection.** The pocket volume trace already shows a pocket opening and closing during MD; naming that as a finding rather than a plot is the next step

---

## 👤 Author

**Marc C. Deller, D.Phil.**  
Structural biologist & drug discovery scientist  

<table>
<tr>
<td>🌐</td><td><a href="https://marcdeller.com" target="_blank" rel="noopener noreferrer">marcdeller.com</a></td>
<td>✉️</td><td><a href="mailto:marc@marcdeller.com">marc@marcdeller.com</a></td>
<td>🐙</td><td><a href="https://github.com/bellcheddar/GOBSMACKED" target="_blank" rel="noopener noreferrer">github.com/bellcheddar/GOBSMACKED</a></td>
</tr>
</table>

Released under the MIT licence. The app name is a joke; the numbers are not.
