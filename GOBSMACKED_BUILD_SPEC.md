# GOBSMACKED: Build Specification for Claude Code

> **Ground-truth Overlay for Binding Sites, Modes And Complex Kinetics/Dynamics**
> Fold, dock, relax and annotate a protein-ligand complex, then check it against the experimental structure.

Author: Marc C. Deller, D.Phil. (marc@marcdeller.com, GitHub bellcheddar)
Target host: gobsmacked.mdeller.com (Flask/gunicorn behind nginx, single DigitalOcean droplet, CPU only, ~3.8 GB RAM)
Repo: github.com/bellcheddar/GOBSMACKED
Build budget: one day, four phases. Each phase must end in a runnable, deployable state.

---

## 0. Read this first

You are building a two-stage web app. The droplet cannot run ESMFold, the PandaDock GNN or OpenMM MD, so all heavy compute is packaged into a downloadable **run bundle** that the user executes on a local GPU machine. The user then uploads the **results archive** for analysis.

```
PREPARE (droplet, CPU)                     RUN (user's GPU box)                 ANALYZE (droplet, CPU)
--------------------------------           ------------------------------       --------------------------------
UniProt / sequence / PDB in                pixi run gobsmacked                  upload results.tar.gz
fetch pre-folded model if any              ESMFold (if needed)                  PLIP + PandaMap (server)
annotate (UniProt, Pfam, KLIFS, GPCRdb)    PDBFixer prep                        superpose vs reference
pocket picker (Mol*)                       RDKit ligand 3D                      scorecard + mode classifier
reference structure picker (RCSB)          PandaDock hybrid/flex                Mol* three-way overlay (client)
emit run_bundle.tar.gz                     OpenMM minimise + short MD           Plotly trajectory plots (client)
                                           MDTraj summary -> results.tar.gz
```

Three prior apps supply reusable code. Copy, do not reinvent:

| Source | Reuse |
|---|---|
| BoltzMaker (github.com/bellcheddar/BoltzMaker) | YAML campaign manifest pattern, pixi packaging, Mol* pocket picker, PLIP wrapper, SSE comparison, MPS crash recovery |
| FlexAppeal | OpenMM minimise/MD driver, MDTraj RMSD/RMSF analytics, Plotly trajectory panels |
| AlphaFraud (github.com/bellcheddar/AlphaFraud) | Flask/gunicorn/nginx/SQLite deployment layout, RCSB and AFDB fetch helpers, superposition and scoring helpers (biotite, tmtools) |

House rules: British English (licence, colour, minimise, analyse). No em dashes: use colons or parentheses. Forbidden words: groundbreaking, revolutionary, paradigm-shifting, cutting-edge, seamless, robust (as filler). marcdeller.com brand theme on every page (see Section 9).

---

## 1. Concept

GOBSMACKED unifies three questions about a predicted protein-ligand complex:

1. **Verify (PIER REVIEW):** how close does fold + dock + relax land to the crystal? A reproducibility scorecard, graded A to F.
2. **Induced fit (HOLOGRAM):** ESMFold gives an apo-like pocket. Does flexible docking plus MD relaxation move the pocket toward the holo structure? A three-way overlay (predicted model, MD-relaxed complex, experimental holo) and a rescue score.
3. **Binding mode (GATEKEEPER):** Pfam routes the annotation layer. Kinases get KLIFS numbering, DFG/alphaC states and Type I / I½ / II labels. GPCRs get GPCRdb generic numbering, orthosteric vs allosteric site and microswitch states. The predicted mode label is compared to the experimental mode label.

Success is a user being able to paste a UniProt ID and a SMILES, download a bundle, run it, upload the results and see, within one page, whether the prediction is trustworthy and why.

---

## 2. Stack and licences

| Component | Role | Licence | Where it runs |
|---|---|---|---|
| Flask 3, gunicorn, nginx, certbot | Web app | BSD/MIT | droplet |
| SQLite | Job and annotation cache | public domain | droplet |
| biotite, tmtools, numpy, scipy | Superposition, scoring | BSD/MIT | droplet |
| RDKit | SMILES to 3D, MCS, Tanimoto, PoseBusters-style checks | BSD-3 | droplet + bundle |
| PLIP | Interaction fingerprints | GPL-2.0 | droplet only, subprocess. Never vendor into the bundle or a binary. |
| PandaMap (pip `pandamap`) | 2D interaction map, empirical dG | MIT | droplet |
| PandaDock (pip `pandadock[gnn]`) | Docking (hybrid = pose search + SE(3) GNN rescoring) | MIT | bundle |
| ESMFold (facebookresearch/esm, OpenFold dep) | Folding when no model exists | MIT / Apache-2.0 | bundle |
| OpenMM, PDBFixer, openmmforcefields (OpenFF) | Prep, minimise, MD | MIT + LGPL / MIT / MIT | bundle |
| MDTraj | Trajectory analysis | LGPL-2.1 | bundle + droplet. Do not use MDAnalysis (GPL-2.0). |
| Mol* | 3D viewer | MIT | browser |
| Plotly.js 2.35.2 | Plots | MIT | browser |
| KLIFS, GPCRdb, InterPro, UniProt, RCSB, AFDB, ESM Atlas REST APIs | Annotation and reference structures | open, attribution | droplet, cached |

**Environments (decided).** Run bundle: **pixi** with a committed `pixi.lock` (conda-forge: openmm, openmmforcefields, ambertools, pdbfixer, rdkit, mdtraj, pytorch + pyg with CUDA feature, esm via pip). Droplet app: **uv** with `pyproject.toml` + `uv.lock` (Flask stack, biotite, rdkit, mdtraj, plip, pandamap). Never mix the two: the droplet must not need torch.

Repository licence: MIT. Add a `THIRD_PARTY.md` listing the above with links. Cite PandaDock (Panda 2026, bioRxiv 10.64898/2026.08.19.745667), PandaMap (Panda 2025), PLIP, KLIFS, GPCRdb and ESMFold in the About tab.

---

## 3. Repository layout

```
GOBSMACKED/
  README.md                      house-standard README (marcs-vibe-coding template)
  LICENSE                        MIT
  THIRD_PARTY.md
  pyproject.toml                 droplet app deps only
  app/
    __init__.py                  create_app(), blueprints, SQLite init
    config.py
    routes/
      prepare.py                 /prepare, /api/fetch, /api/annotate, /api/references, /api/bundle
      analyze.py                 /analyze, /api/upload, /api/results/<job_id>
      runs.py                    /runs, /runs/<job_id>, /api/runs, /api/runs/<job_id>/visibility
      about.py                   /about, static pipeline schematic + software table
    services/
      fetch.py                   UniProt, AFDB, ESM Atlas, RCSB helpers (from AlphaFraud)
      annotate.py                Pfam/InterPro, UniProt features, KLIFS, GPCRdb; family router
      references.py              RCSB holo search by UniProt + ligand Tanimoto
      bundle.py                  writes run_bundle.tar.gz
      ingest.py                  validates results.tar.gz against schema
      superpose.py               pocket superposition, RMSDs, chi1 (biotite)
      interactions.py            PLIP subprocess, PandaMap API, fingerprint Jaccard
      scorecard.py               grades and the composite GOBSMACK score
      modes.py                   kinase and GPCR binding-mode classifiers
      dynamics.py                MDTraj summaries, rescue score
    db.py                        SQLite schema and helpers
    templates/                   Jinja2, one base.html with brand header
    static/
      css/gobsmacked.css
      js/pocket_picker.js        from BoltzMaker
      js/viewer.js               Mol* wrapper, three-way overlay, sync
      js/plots.js                Plotly panels
      icon/gobsmacked.svg, .png  from marcs-vibe-icon
  bundle_template/               copied verbatim into every run bundle
    pixi.toml
    run.py
    gobsmacked_run/
      fold.py                    ESMFold or skip
      prep.py                    PDBFixer, protonation, ligand 3D (RDKit)
      dock.py                    pandadock hybrid / flex
      md.py                      OpenMM minimise + MD (from FlexAppeal)
      summarise.py               MDTraj -> summary.json
      schema.py                  results schema version
  design/
    gobsmacked_design_4_instrument_panel.html   visual contract for Section 9
  deploy/
    gobsmacked.service           gunicorn systemd unit
    nginx.conf
  tests/
    test_ingest.py, test_scorecard.py, test_modes.py, fixtures/
```

---

## 4. Data contracts

### 4.1 campaign.yaml (written by Prepare, read by run.py)

```yaml
gobsmacked_version: "1.0"
job_id: "gs_20260905_ab12cd"
protein:
  uniprot: "P00533"                      # or null
  sequence: "MRPSGTAGAALLALLAALCPASRALEEKK..."
  source_structure: "afdb"                # afdb | esm_atlas | pdb | user_pdb | fold
  source_id: "AF-P00533-F1"               # PDB id, AFDB id, or null when folding
  chain: "A"
  residue_range: [696, 1022]              # optional trim to the domain
  family: "kinase"                        # kinase | gpcr | other (from Pfam router)
ligand:
  name: "erlotinib"
  smiles: "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"
  protonation_ph: 7.4
pocket:
  method: "residues"                      # residues | center_box | reference_ligand
  residues: ["A:718", "A:745", "A:790", "A:855"]
  center: [x, y, z]                       # filled in by Prepare from residues
  box: [22, 22, 22]
reference:
  pdb_id: "1M17"
  chain: "A"
  ligand_ccd: "AQ4"
  apo_pdb_id: null                        # optional apo reference
docking:
  mode: "hybrid"                          # hybrid | flex | dock
  exhaustiveness: 16
  num_poses: 10
  flexible_residues: "auto"               # auto = pocket residues within 4 A of box centre
md:
  forcefield: "amber14"
  ligand_forcefield: "openff-2.1.0"
  minimise_steps: 5000
  equilibration_ps: 100
  production_ps: 1000
  frame_interval_ps: 10
  platform: "auto"                        # CUDA | OpenCL | CPU | auto
```

### 4.2 results.tar.gz (written by run.py, consumed by Analyze)

```
results/
  manifest.json          {schema: "1.0", job_id, campaign_sha256, engine_versions{}, timings{}, warnings[]}
  campaign.yaml          echoed back
  model_apo.pdb          the input or folded structure, prepared (hydrogens, no ligand)
  plddt.json             per-residue pLDDT when folded (else null)
  poses/
    poses.sdf            all docked poses, PandaDock scores in SD tags
    scores.csv           pose_id, pandadock_score, gnn_affinity, rank
  complex_pose1.pdb      top pose merged with receptor, before MD
  complex_min.pdb        after minimisation
  complex_md_final.pdb   last MD frame
  traj/
    traj.dcd             production trajectory, protein + ligand only, aligned on pocket
    topology.pdb
    summary.json         per-frame ligand RMSD to pose1, protein RMSD, per-residue RMSF, pocket volume, contact matrix (residue x frame, 4 A heavy-atom)
  logs/
    run.log
```

`ingest.py` must reject archives that fail schema validation with a message naming the missing file. Extra files are ignored.

### 4.3 SQLite

```
jobs(job_id PK, created, updated, title, uniprot, ligand_name, ligand_smiles, family, reference_pdb, status, visibility, owner_token, campaign_yaml TEXT, results_path, scorecard_json TEXT, gobsmack_score REAL, grade TEXT, mode_predicted TEXT, mode_reference TEXT, mode_match INTEGER)
  -- status: prepared | results_uploaded | analysed | failed
  -- visibility: public | private (set at Prepare time, editable by the owner)
  -- owner_token: random 32-char secret issued at Prepare time; stored hashed (sha256). Proves ownership for private runs and visibility changes.
annotation_cache(key PK, source, fetched, payload TEXT)      -- key = source:identifier, TTL 30 days
reference_cache(uniprot PK, fetched, payload TEXT)
```

---

## 5. Prepare tab

Route `/prepare`. Single page, four panels, top to bottom, each unlocking the next.

**Panel 1: Protein.** Text box accepting a UniProt accession, raw sequence or a PDB ID; file drop for PDB/mmCIF. On submit `/api/fetch` returns: canonical sequence, UniProt features, and the best available structure in this priority order: user upload > PDB entry (if given) > AFDB model > ESM Atlas model > "fold in bundle". Show which was chosen, mean pLDDT if applicable, and a Mol* preview.

**Panel 2: Annotation.** `/api/annotate` calls InterPro for Pfam domains and routes on them:
- `PF00069`, `PF07714` (Pkinase, PK_Tyr_Ser-Thr) → `family = kinase`, fetch KLIFS: kinase ID from UniProt, 85-residue pocket numbering mapped to sequence positions, DFG-in/out and alphaC-in/out state of every KLIFS structure for the kinase.
- `PF00001`, `PF00002`, `PF00003`, `PF10324` (7tm_1/2/3, 7TM_GPCR_Srsx) → `family = gpcr`, fetch GPCRdb: generic residue numbering (Ballesteros-Weinstein), ligand-site residues from the family's structures, microswitch positions.
- anything else → `family = other`, UniProt features only (binding sites, active sites, transmembrane).
Render a sequence track (SVG, no library) with UniProt features, Pfam domains and family-specific positions. Clicking a residue highlights it in the Mol* preview and adds it to the pocket selection.

**Panel 3: Ligand and pocket.** SMILES box with RDKit validation and a 2D depiction (RDKit.js from CDN, or server-side SVG). Pocket picker lifted from BoltzMaker: click residues in Mol*, or type a residue list, or pick "use reference ligand site" once a reference is chosen. Compute box centre and size from selected residues (centroid, extent + 8 Å, min 18 Å per side).

**Panel 4: Reference structure.** `/api/references` queries RCSB for entries mapped to the UniProt (or 90 % identity to the sequence), pulls the CCD ligands of each, computes Tanimoto (Morgan r=2, 2048 bits) against the input SMILES, and lists entries sorted by similarity with resolution, ligand name and a 2D depiction. Default selection: highest Tanimoto with resolution < 2.5 Å. Allow "none" (verification disabled) and an optional apo reference. Also allow the user to type a PDB ID directly.

**Visibility.** A two-state control on the Generate bundle panel, defaulting to **Public**: "Public: this run appears in the Runs tab for anyone" / "Private: only people with the run link and owner key can see it". Same pattern as BoltzMaker. On generate, the server issues an `owner_token`, embeds it in `campaign.yaml` under `owner_token`, shows it once in the UI with a copy button, and stores only its sha256. The token rides inside the results archive (campaign.yaml is echoed back), so uploading results to a private run needs no extra typing; changing visibility later or deleting a run requires pasting the token.

**Generate bundle.** `/api/bundle` writes `campaign.yaml`, copies `bundle_template/`, includes `model_apo.pdb` when a structure was fetched, and returns `run_bundle_<job_id>.tar.gz`. Show the three commands the user needs:

```
tar xzf run_bundle_gs_xxx.tar.gz && cd run_bundle_gs_xxx
pixi install
pixi run gobsmacked          # writes results/results.tar.gz
```

---

## 6. Run bundle (run.py)

Sequential stages, each idempotent and resumable (a `.done` marker per stage, same pattern as BoltzMaker MPS recovery). `--stage` flag reruns from a given stage.

1. **fold**: skip if `model_apo.pdb` present. Otherwise ESMFold via `esm.pretrained.esmfold_v1()`; chunk size auto from sequence length; write `plddt.json`. Warn if any pocket residue has pLDDT < 70.
2. **prep**: PDBFixer (missing atoms, hydrogens at `protonation_ph`), remove waters and heteroatoms, renumber consistently with campaign residue numbering (keep author numbering in a mapping file). Ligand: RDKit from SMILES, ETKDG conformer, MMFF minimise, write `ligand.sdf`.
3. **dock**: `pandadock hybrid -r receptor.pdb -l ligand.sdf --center ... --box ... -o poses/` (or `flex` / `dock` per campaign). Parse output into `poses.sdf` and `scores.csv`. Merge top pose into `complex_pose1.pdb`.
4. **md** (from FlexAppeal): OpenMM system with Amber14 + OpenFF ligand params via openmmforcefields `SystemGenerator`; TIP3P solvent, 0.15 M NaCl, 10 Å padding; minimise `minimise_steps`; NVT then NPT equilibration with heavy-atom restraints released over `equilibration_ps`; production `production_ps` at 300 K, 2 fs, HMR. Save DCD every `frame_interval_ps`. Write `complex_min.pdb` and `complex_md_final.pdb`.
5. **summarise**: MDTraj. Align on pocket Cα. Per-frame: ligand heavy-atom RMSD to pose1, protein Cα RMSD, pocket Cα RMSD. Per-residue RMSF. Pocket volume per frame (simple grid: count 1 Å voxels within box that are > 1.4 Å from any protein heavy atom and enclosed, no fpocket dependency). Contact matrix (pocket residues × frames, any heavy atom within 4 Å of ligand). Write `summary.json`, then `results.tar.gz`.

Time budget target on a single consumer GPU: under 30 minutes for a 300-residue domain with the default 1 ns production. Print a wall-clock estimate before starting.

---

## 7. Analyze tab

Route `/analyze`. Upload `results.tar.gz`. Server runs ingest → superpose → interactions → scorecard → modes → dynamics, stores `scorecard_json`, and renders a single results page with a sticky sub-nav: **Scorecard · Complex · Overlay · Dynamics · Mode · Report**.

### 7.1 Superposition (superpose.py)

Reference chain and predicted model are aligned on **pocket Cα atoms** (residues within 8 Å of the reference ligand), not on the whole chain, so fold errors far from the site do not inflate ligand RMSD. Also compute whole-chain TM-score (tmtools) for context. Report:

- `lig_rmsd_pose1`, `lig_rmsd_min`, `lig_rmsd_md_final`: ligand heavy-atom RMSD to reference ligand after pocket superposition, symmetry-corrected (RDKit `GetBestRMS` on the reference ligand graph).
- `pocket_ca_rmsd_model`, `pocket_ca_rmsd_md_final`
- `pocket_sc_rmsd_*`: pocket side-chain heavy-atom RMSD
- `chi1_agreement`: fraction of pocket residues whose χ1 is within 40° of the reference (model and MD-final)
- `tm_score_chain`

If no reference was chosen, these fields are null and the Scorecard shows an "unverified" banner.

### 7.2 Interactions (interactions.py)

Run PLIP (subprocess, XML out) on `complex_pose1.pdb`, `complex_md_final.pdb` and the reference holo. Convert each to a fingerprint set of `(interaction_type, residue_number)` tuples using reference numbering via the mapping file. Compute Jaccard(predicted, reference) for pose1 and MD-final. Run PandaMap on the same three files for the 2D maps (PNG) and its empirical dG. Show a three-column interaction table with ticks and crosses per residue.

### 7.3 Scorecard (scorecard.py)

| Metric | A | B | C | D | F |
|---|---|---|---|---|---|
| Ligand RMSD, best of (pose1, MD-final) | ≤ 1.0 Å | ≤ 2.0 | ≤ 3.0 | ≤ 4.0 | > 4.0 |
| Pocket Cα RMSD, MD-final | ≤ 0.8 | ≤ 1.2 | ≤ 1.8 | ≤ 2.5 | > 2.5 |
| χ1 agreement, MD-final | ≥ 0.85 | ≥ 0.70 | ≥ 0.55 | ≥ 0.40 | < 0.40 |
| PLIP Jaccard, best of (pose1, MD-final) | ≥ 0.75 | ≥ 0.55 | ≥ 0.40 | ≥ 0.25 | < 0.25 |
| Pose validity (PoseBusters-style: clashes < 2.2 Å, bond lengths, chirality preserved, ligand inside box) | pass | | | | fail |
| MD stability: ligand RMSD drift (last 200 ps mean minus first 200 ps mean) | ≤ 0.5 Å | ≤ 1.0 | ≤ 1.5 | ≤ 2.5 | > 2.5 |
| Rescue score: pocket Cα RMSD-to-ref (model) minus (MD-final) | ≥ +0.5 Å | ≥ +0.2 | ≥ 0 | ≥ −0.3 | < −0.3 |

Composite **GOBSMACK score** (0 to 100): weighted mean of graded metrics (ligand RMSD 30, Jaccard 20, pocket Cα 15, χ1 10, stability 10, rescue 10, validity gate 5). A validity fail caps the composite at 40. Render as a big grade tile plus a radar plot. Add one plain-English sentence per metric explaining what it means and what to do about it (for example "Ligand is in the right pocket but flipped: check symmetry and hinge H-bonds").

### 7.4 Complex view

Mol* with `complex_md_final.pdb` loaded, ligand as ball-and-stick, pocket residues as sticks, PLIP interactions drawn as dashed lines coloured by type (reuse BoltzMaker's PLIP-to-Mol* layer). Toggle buttons: pose1 / minimised / MD-final. PandaMap 2D PNG alongside.

### 7.5 Overlay view (HOLOGRAM)

Three states superposed on pocket Cα: **Model** (grey), **MD-final** (marcdeller blue), **Reference** (amber), plus optional **Apo reference** (purple). Ligands from MD-final and reference both shown. Side chains of pocket residues coloured by displacement between MD-final and reference (white to red, 0 to 3 Å). Checkbox per state. Camera locked to the pocket. Below it: a table of the ten most-displaced pocket residues with model, MD-final and reference χ1.

### 7.6 Dynamics view

Plotly panels from `summary.json` (FlexAppeal layouts): ligand RMSD to pose1 and **to reference ligand** over time (two traces, this second one is the induced-fit story); protein and pocket Cα RMSD; per-residue RMSF with pocket residues shaded; pocket volume over time with the reference pocket volume as a horizontal line; contact persistence heatmap (residue × time) with a persistence bar (fraction of frames in contact) and the reference contacts marked.

### 7.7 Mode view (GATEKEEPER)

`modes.py` classifies both the predicted complex and the reference:

**Kinase.** Using KLIFS pocket numbering mapped onto the sequence:
- DFG state from the D of DFG (KLIFS 81) side-chain position relative to the αC glutamate (KLIFS 24) and gatekeeper (KLIFS 45); αC state from the K(17)–E(24) salt-bridge distance. Use the Dunbrack/Modi–Dunbrack style distance criteria (D1 = Cα(DFG-Phe) to Cα(αC-Glu + 4), D2 = Cζ(DFG-Phe) to Cα(K17)); label DFG-in / DFG-out / DFG-inter.
- Subpocket occupancy: which KLIFS regions (hinge, gatekeeper, back pocket I/II, front pocket, solvent-exposed, DFG, αC) have ligand heavy atoms within 4 Å.
- Type label: I (hinge, DFG-in, no back pocket II), I½ (hinge plus back pocket, DFG-in), II (back pocket II, DFG-out), allosteric (no hinge contact).
- Report hinge H-bond count and gatekeeper contact.

**GPCR.** Using GPCRdb generic numbering:
- Site label: orthosteric (contacts to the family's consensus orthosteric residues), extracellular vestibule, intracellular (contacts with TM6/TM7 cytoplasmic ends, ICL), lipid-facing.
- Microswitch states measured on the model and reference: TM3–TM6 distance (3.50 to 6.30 Cα), NPxxY (7.49 to 7.53) RMSD to inactive template, PIF/toggle (6.48) χ angles. Label active-like / inactive-like.

**Other.** Report contact residues with UniProt feature names (binding site, active site, metal).

Output a side-by-side card: predicted mode vs reference mode, with a match/mismatch verdict that feeds the Report tab as a categorical check separate from the numeric scorecard.

### 7.8 Report

Single printable page: scorecard, mode verdict, three-way overlay screenshot (Mol* `getImage`), interaction table, key dynamics plots, engine versions and timings, citations. "Download report (HTML)" and "Download all (zip)" buttons.

---

### 7.9 Runs tab

Route `/runs`. Every Prepare creates a row; every Analyze fills it in. The page is a table in the deck's panel style, one row per run, newest first, with a search box (UniProt, ligand name, PDB ID, job ID) and filters for family (kinase / GPCR / other), status and grade.

Columns: job ID (link to `/runs/<job_id>`), created, protein (UniProt + name), ligand (name + 2D thumbnail from RDKit SVG), family, reference PDB, status as a stage-strip miniature (seven 6 px cells coloured by state), GOBSMACK score and grade tile, mode verdict LED (green match / amber differ / grey unverified). Rows for private runs are listed only when the viewer holds the owner token (kept in `localStorage` after Prepare or after pasting it once); otherwise they are absent, not greyed, so private runs leak nothing, not even their count.

`/runs/<job_id>` is the Analyze results page for that run (Scorecard, Complex, Overlay, Dynamics, Mode, Report), plus a run header with: title (editable by owner), visibility badge and toggle (owner only, requires token), "Download bundle", "Download results", "Download report", "Delete run" (owner only, confirm by retyping the job ID). A private run's URL is unguessable (job IDs are 12 random base32 characters) but the page still checks the token before rendering anything other than "This run is private. Paste the owner key to view it."

`/api/runs` returns JSON for the table (public rows, plus private rows whose token hash matches the supplied `X-Owner-Token` header). `/api/runs/<job_id>/visibility` PATCH flips public/private with the token. Results archives and bundles live under `data/runs/<job_id>/`; a nightly systemd timer prunes `results.tar.gz` for public runs older than 90 days but keeps `scorecard_json`, the report HTML and the final PDBs so the row stays useful.

The stage strip on `/runs` shows aggregate counts (runs prepared, results uploaded, analysed, failed) instead of one job's stages.

### 7.10 About tab

Route `/about`. Static page, same panel style, three parts.

**Pipeline schematic.** An inline SVG (hand-drawn in the repo, not generated at runtime) showing the two hosts and seven stages as the stage strip does, with data flowing left to right: inputs (UniProt / sequence / PDB / SMILES) → Prepare on the droplet (fetch, annotate, pocket, reference) → `run_bundle.tar.gz` → Run on the user's GPU (fold, prep, dock, MD, summarise) → `results.tar.gz` → Analyze on the droplet (superpose, interactions, scorecard, modes, dynamics) → Runs and Report. Each stage box lists the software it uses in `--muted` text. Colour the droplet lane slate and the GPU lane `--panel-2`; files as amber chips. Under the schematic, three short paragraphs: what GOBSMACKED measures, why superposition is on pocket Cα atoms, and what the grades mean (the thresholds table from 7.3, reproduced).

**Software and references.** One table, grouped by stage, with columns: tool, what it does here, version pinned in the bundle or app, licence, links. Links column carries the GitHub repository and the paper DOI where one exists. Populate from `THIRD_PARTY.md` so the two never drift (generate the table at build time from a single `software.yaml`). Required entries:

| Tool | Repository | Reference |
|---|---|---|
| ESMFold (esm) | github.com/facebookresearch/esm | Lin et al. 2023, Science, doi:10.1126/science.ade2574 |
| OpenFold | github.com/aqlaboratory/openfold | Ahdritz et al. 2024, Nat Methods, doi:10.1038/s41592-024-02272-z |
| PandaDock | github.com/pritampanda15/PandaDock | Panda 2026, bioRxiv, doi:10.64898/2026.08.19.745667 |
| PandaMap | github.com/pritampanda15/PandaMap | Panda 2025, PandaMap (PyPI) |
| PLIP | github.com/pharmai/plip | Adasme et al. 2021, NAR, doi:10.1093/nar/gkab294 |
| OpenMM | github.com/openmm/openmm | Eastman et al. 2024, J Phys Chem B, doi:10.1021/acs.jpcb.3c04565 |
| PDBFixer | github.com/openmm/pdbfixer | (OpenMM project) |
| openmmforcefields / OpenFF | github.com/openmm/openmmforcefields, github.com/openforcefield/openff-toolkit | Boothroyd et al. 2023, JCTC, doi:10.1021/acs.jctc.3c00039 |
| MDTraj | github.com/mdtraj/mdtraj | McGibbon et al. 2015, Biophys J, doi:10.1016/j.bpj.2015.08.015 |
| RDKit | github.com/rdkit/rdkit | rdkit.org |
| biotite | github.com/biotite-dev/biotite | Kunzmann & Hamacher 2018, BMC Bioinf, doi:10.1186/s12859-018-2367-z |
| tmtools / TM-align | github.com/jvkersch/tmtools | Zhang & Skolnick 2005, NAR, doi:10.1093/nar/gki524 |
| Mol* | github.com/molstar/molstar | Sehnal et al. 2021, NAR, doi:10.1093/nar/gkab314 |
| Plotly.js | github.com/plotly/plotly.js | plotly.com |
| KLIFS | klifs.net | Kanev et al. 2021, NAR, doi:10.1093/nar/gkaa895 |
| GPCRdb | gpcrdb.org, github.com/protwis | Pándy-Szekeres et al. 2023, NAR, doi:10.1093/nar/gkac1013 |
| InterPro / Pfam | ebi.ac.uk/interpro | Paysan-Lafosse et al. 2023, NAR, doi:10.1093/nar/gkac993 |
| UniProt | uniprot.org | The UniProt Consortium 2023, NAR, doi:10.1093/nar/gkac1052 |
| RCSB PDB | rcsb.org | Burley et al. 2023, NAR, doi:10.1093/nar/gkac1077 |
| AlphaFold DB | alphafold.ebi.ac.uk | Varadi et al. 2024, NAR, doi:10.1093/nar/gkad1011 |
| ESM Metagenomic Atlas | esmatlas.com | Lin et al. 2023 (as ESMFold) |
| PoseBusters (criteria) | github.com/maabuu/posebusters | Buttenschoen et al. 2024, Chem Sci, doi:10.1039/D3SC04185A |
| Flask, gunicorn, nginx, SQLite, pixi, uv | respective repositories | version only |

Verify every DOI resolves at build time (a `make check-refs` target that HEAD-requests each) and fix any that do not; do not ship a reference that has not been checked.

**Credits and licence.** Author, contact, repo link, MIT licence, how to cite GOBSMACKED itself (repo URL and version), and one line acknowledging that PLIP is GPL-2.0 and run server-side only.

Both tabs sit in the top navigation next to Prepare and Analyze: **Prepare · Analyze · Runs · About**, rendered as four text links in the brand line's centre, current tab in `--phos`.

---

## 8. Four phases, one day

**Phase 1 (2 h): Prepare without references.** Flask skeleton from AlphaFraud, brand base template, `/prepare` panels 1 to 3, fetch and annotate services with SQLite cache, pocket picker, `/api/bundle` emitting a valid `campaign.yaml` and the untouched bundle template. Includes the visibility control and owner token, the `jobs` row insert, and a first-pass `/runs` table listing prepared runs. Definition of done: paste P00533 + erlotinib SMILES, pick the ATP site, choose Private, download a bundle that `pixi install`s, and see the run listed in Runs only when the token is present.

**Phase 2 (2 h): The bundle runs.** `run.py` stages 1 to 5 with the resume markers. Test on EGFR kinase domain (1M17 chain A, erlotinib) with `production_ps: 100` for speed. Definition of done: `results.tar.gz` validates against `ingest.py`.

**Phase 3 (2 h): Analyze without verification.** Upload, ingest, PLIP, PandaMap, Complex view, Dynamics view. Definition of done: the EGFR results render every panel with no reference selected.

**Phase 4 (2 h): Verify, overlay, modes, ship.** Panel 4 reference search, `superpose.py`, `scorecard.py`, Overlay view, `modes.py` with **both kinase and GPCR paths complete** (GPCR annotation is a day-one requirement, not a stretch), Report, `/runs/<job_id>` result pages with owner controls, the About tab with schematic and checked reference table, README, icon, systemd + nginx, nightly prune timer, deploy to gobsmacked.mdeller.com. Definition of done: EGFR/erlotinib scores against 1M17 and is labelled Type I, DFG-in, in both predicted and reference; β2-adrenergic receptor (P07550) with carazolol scores against 2RH1 and is labelled orthosteric, inactive-like, in both.

To protect the GPCR path within the day: build `annotate.py`'s GPCRdb client in Phase 1 alongside KLIFS (same cache pattern), and write `tests/test_modes.py` fixtures for 2RH1 before Phase 4 starts. Second fixture archive: β2AR/carazolol, 20 frames.

---

## 9. UI design brief: "Instrument panel"

The reference mockup is `design/gobsmacked_design_4_instrument_panel.html` (copy it into the repo under `design/` and treat it as the visual contract). GOBSMACKED is deliberately not another flat white page with blue cards: it is a dark control room for reading a pipeline's instruments. Build every page from this brief, not from the marcdeller.com HTML app template used by the other apps.

### 9.1 Tokens

```css
:root {
  --slate:   #1b2433;   /* page background: slate blue, never near-black */
  --panel:   #243044;   /* cards */
  --panel-2: #2c3a50;   /* nested surfaces, gauges, toggles */
  --line:    #3a4a63;   /* borders, gauge tracks */
  --phos:    #5de1e6;   /* phosphor cyan: primary data colour, predicted / MD-final */
  --amber:   #ffb454;   /* experimental reference, warnings */
  --red:     #ff5c5c;   /* failure, thresholds */
  --green:   #7ee2a8;   /* pass, stage complete, mode match */
  --grey:    #9fb0c7;   /* AFDB / ESMFold input model, muted text */
  --purple:  #c39cff;   /* optional apo reference */
  --text:    #e8edf5;
  --muted:   #9fb0c7;
}
```

Typeface: **Sora** (Google Fonts) for everything, weights 300/400/600/800. No second face; code and residue IDs use Sora at 400 with `font-variant-numeric: tabular-nums`. No all-caps labels, no eyebrow labels, no middle-dot meta strings.

Grade colours: A green, B phos, C amber, D #ff8a5c, F red. Gauge fill uses the grade colour.

### 9.2 Layout

Every page opens with the **stage strip**: seven equal cells (Fetch, Annotate, Fold, Dock, MD, Verify, Mode), each with a 3 px top border coloured by state (green done, amber attention, red failed, grey pending or skipped) and one line of status text. On Prepare the strip shows which stages the bundle will run; on Analyze it shows what the archive contained and timings. It is the app's navigation and progress indicator in one.

Below the strip, a three-column **deck** at ≥ 1100 px: left 320 px (Score), centre fluid (Scope: Mol* and Plotly), right 340 px (Mode switches). Collapses to one column on narrow screens in the order Score, Scope, Mode. Cards: `--panel` background, 1 px `--line` border, 6 px radius, 18 px padding. Card headings are 13 px, weight 600, `--muted`, sentence case.

Prepare uses the same deck: left column holds the input and ligand forms, centre the sequence track and Mol* pocket picker, right the annotation and reference lists. The "Generate bundle" button is phos-filled, slate text, and is the only filled button on the page.

### 9.3 Components

- **Score dial**: semicircular SVG arc, `--line` track, grade-coloured fill proportional to the GOBSMACK score, grade letter centred in the arc at 44 px weight 800, score number below at 56 px weight 800, one-line label in `--muted`. Under it, a 2 × 3 grid of gauges: label 11 px muted, value 20 px weight 600, 4 px bar.
- **Scope**: Mol* container on `#111826` with a faint 32 px cyan grid (`rgba(93,225,230,.06)`), a HUD line top-left in phos 11 px stating the superposition basis and TM-score, and toggle buttons bottom-left (Model / MD final / Reference / Apo ref) styled as `--panel-2` chips with a phos border when pressed (`aria-pressed`). State colours in Mol* follow the tokens: model grey, MD-final phos, reference amber, apo purple. Displacement colouring for side chains: `--panel-2` to `--red`.
- **Traces**: Plotly on transparent backgrounds, Sora 12 px, line colour phos, thresholds as red dashed lines labelled at the left edge, reference values as amber dashed lines. Grid lines `--line`. No mode bar logo, no lasso or select.
- **Switch list**: rows with an LED (10 px circle with a matching 8 px glow) and a right-aligned value. Green = present and matching the reference, amber = differs from the reference, grey = absent on both. Kinase rows: DFG, αC helix, hinge H-bonds, gatekeeper contact, back pocket I, back pocket II, front pocket, solvent-exposed. GPCR rows: site (orthosteric / vestibule / intracellular / lipid-facing), TM3–TM6 distance, NPxxY RMSD, toggle switch χ, PIF motif, sodium site contact. Below the list, a two-cell "Predicted / Crystal" pair showing the mode label in phos, then a full-width verdict bar: green tint and border for "Binding mode matches", amber for "Binding mode differs" with one sentence saying how.
- **Sequence track** (Prepare): horizontal SVG, residues as 6 px cells, Pfam domains as bands under the track, KLIFS or GPCRdb positions as phos ticks with labels on hover, selected pocket residues filled phos. Clicking a cell toggles it in the pocket set and highlights it in Mol*.
- **Report**: same components, light-on-dark preserved when printing (set `color-adjust: exact`); a "Download report" chip and a "Download all" chip in the stage strip's right end.

### 9.4 Motion and states

One motion only: when Analyze finishes, the stage strip cells switch to green left to right at 80 ms intervals and the score dial fills once. Everything else changes instantly. Respect `prefers-reduced-motion`. Empty Scope shows "Upload a results archive to see the complex" with the upload control inside the scope area. Errors are one sentence naming what failed and the fix, in the same panel style with a red top border.

### 9.5 Brand

A single 12 px brand line above the stage strip: "Marc C. Deller, D.Phil." linking to marcdeller.com on the left, "GOBSMACKED · job id" centred, marc@marcdeller.com on the right, in `--muted` with a `--line` bottom border. No gradient header, no logo dot, no footer beyond a one-line citation strip on the Report page. The app icon (marcs-vibe-icon skill) sits at 24 px beside the app name.

Tone of all copy: plain, slightly dry, no marketing adjectives. The app name is a joke; the numbers are not.

---

## 10. Testing and acceptance

- `tests/fixtures/` holds a pre-computed EGFR results archive (small: 20 frames) so Analyze tests run on the droplet in seconds.
- `test_ingest.py`: valid archive passes; missing `summary.json` fails with a named error.
- `test_scorecard.py`: reference complex scored against itself gives grade A across the board and GOBSMACK = 100; a pose translated 5 Å gives ligand RMSD F.
- `test_modes.py`: 1M17 labelled Type I, DFG-in, αC-in; 1IEP (imatinib, Abl) labelled Type II, DFG-out.
- `test_runs.py`: private run absent from `/api/runs` without the token, present with it; visibility PATCH rejected with a wrong token; `/runs/<job_id>` for a private run renders only the key prompt.
- `make check-refs` passes with every DOI in `software.yaml` resolving.
- Manual: the full Prepare → bundle → run → Analyze loop on EGFR/erlotinib completes and the composite is ≥ 70 when the AFDB model is used as input.

---

## 11. Out of scope for v1 (note in README as roadmap)

Multi-ligand SAR series and ChEMBL affinity correlation (STEVEDORE), multi-engine ingestion of Boltz-2/Vina/DiffDock poses (DOCKYARD), FEP hand-off, batch leaderboard mode, cryptic-pocket detection beyond the simple volume trace.
