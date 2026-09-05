# Third-party software

GOBSMACKED itself is MIT licensed. This file lists everything it uses, what it
uses it for, and under what licence.

**PLIP is GPL-2.0.** It is invoked as a subprocess and never imported, so no GPL
code is linked into this application, and it is not vendored into the run
bundle. Everything else here is permissively licensed or public data.

This file is generated from `software.yaml` by `make third-party`. Edit that file,
not this one.


## Prepare

| Tool | Role | Licence | Repository | Reference |
|---|---|---|---|---|
| UniProt | Canonical sequence, features, protein naming | CC BY 4.0 | [www.uniprot.org](https://www.uniprot.org) | [The UniProt Consortium 2023, Nucleic Acids Research](https://doi.org/10.1093/nar/gkac1052) |
| InterPro / Pfam | Domain assignment, and the family router that follows from it | CC0 | [www.ebi.ac.uk/interpro](https://www.ebi.ac.uk/interpro) | [Paysan-Lafosse et al. 2023, Nucleic Acids Research](https://doi.org/10.1093/nar/gkac993) |
| KLIFS | 85-residue kinase pocket numbering, DFG and alphaC states | open, attribution | [klifs.net](https://klifs.net) | [Kanev et al. 2021, Nucleic Acids Research](https://doi.org/10.1093/nar/gkaa895) |
| GPCRdb | Ballesteros-Weinstein generic numbering, segments, microswitches | Apache-2.0 | [github.com/protwis](https://github.com/protwis) | [Pandy-Szekeres et al. 2023, Nucleic Acids Research](https://doi.org/10.1093/nar/gkac1013) |
| RCSB PDB | Reference structure search, entry metadata, coordinates | CC0 | [www.rcsb.org](https://www.rcsb.org) | [Burley et al. 2023, Nucleic Acids Research](https://doi.org/10.1093/nar/gkac1077) |
| AlphaFold DB | Predicted starting structures and their pLDDT | CC BY 4.0 | [alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk) | [Varadi et al. 2024, Nucleic Acids Research](https://doi.org/10.1093/nar/gkad1011) |
| ESM Metagenomic Atlas | Folding short sequences without a GPU | MIT | [esmatlas.com](https://esmatlas.com) | [Lin et al. 2023, Science](https://doi.org/10.1126/science.ade2574) |
| RDKit | SMILES validation, depiction, fingerprints, symmetry-aware RMSD | BSD-3-Clause | [github.com/rdkit/rdkit](https://github.com/rdkit/rdkit) | rdkit.org |

## Run

| Tool | Role | Licence | Repository | Reference |
|---|---|---|---|---|
| ESMFold (esm) | Folding the sequence when no model exists | MIT | [github.com/facebookresearch/esm](https://github.com/facebookresearch/esm) | [Lin et al. 2023, Science](https://doi.org/10.1126/science.ade2574) |
| OpenFold | ESMFold's structure module | Apache-2.0 | [github.com/aqlaboratory/openfold](https://github.com/aqlaboratory/openfold) | [Ahdritz et al. 2024, Nature Methods](https://doi.org/10.1038/s41592-024-02272-z) |
| PDBFixer | Missing atoms, protonation, cleaning before MD | MIT | [github.com/openmm/pdbfixer](https://github.com/openmm/pdbfixer) | The OpenMM project |
| PandaDock | Pose search and SE(3) GNN rescoring | MIT | [github.com/pritampanda15/PandaDock](https://github.com/pritampanda15/PandaDock) | [Panda 2026, bioRxiv](https://doi.org/10.64898/2026.08.19.745667) |
| OpenMM | Minimisation, equilibration and production MD | MIT | [github.com/openmm/openmm](https://github.com/openmm/openmm) | [Eastman et al. 2024, Journal of Physical Chemistry B (OpenMM 8)](https://doi.org/10.1021/acs.jpcb.3c06662) |
| openmmforcefields / OpenFF | Ligand parameters (OpenFF 2.1.0) alongside Amber14 | MIT | [github.com/openforcefield/openff-toolkit](https://github.com/openforcefield/openff-toolkit) | [Boothroyd et al. 2023, Journal of Chemical Theory and Computation](https://doi.org/10.1021/acs.jctc.3c00039) |
| MDTraj | Trajectory analysis, in the bundle and on the server | LGPL-2.1 | [github.com/mdtraj/mdtraj](https://github.com/mdtraj/mdtraj) | [McGibbon et al. 2015, Biophysical Journal](https://doi.org/10.1016/j.bpj.2015.08.015) |

## Analyze

| Tool | Role | Licence | Repository | Reference |
|---|---|---|---|---|
| PLIP | Interaction fingerprints, run as a subprocess only | GPL-2.0 | [github.com/pharmai/plip](https://github.com/pharmai/plip) | [Adasme et al. 2021, Nucleic Acids Research](https://doi.org/10.1093/nar/gkab294) |
| PandaMap | 2D interaction maps and an empirical binding-energy estimate | MIT | [github.com/pritampanda15/PandaMap](https://github.com/pritampanda15/PandaMap) | Panda 2025, PandaMap |
| biotite | Sequence alignment behind every residue-number mapping | BSD-3-Clause | [github.com/biotite-dev/biotite](https://github.com/biotite-dev/biotite) | [Kunzmann and Hamacher 2018, BMC Bioinformatics](https://doi.org/10.1186/s12859-018-2367-z) |
| tmtools / TM-align | Whole-chain TM-score, reported as context for the pocket fit | MIT | [github.com/jvkersch/tmtools](https://github.com/jvkersch/tmtools) | [Zhang and Skolnick 2005, Nucleic Acids Research](https://doi.org/10.1093/nar/gki524) |
| PoseBusters (criteria) | The pose-validity checks, reimplemented on this app's own geometry | BSD-3-Clause | [github.com/maabuu/posebusters](https://github.com/maabuu/posebusters) | [Buttenschoen et al. 2024, Chemical Science](https://doi.org/10.1039/D3SC04185A) |
| gemmi | Structure input and output, mmCIF and PDB | MPL-2.0 | [github.com/project-gemmi/gemmi](https://github.com/project-gemmi/gemmi) |  |

## Browser

| Tool | Role | Licence | Repository | Reference |
|---|---|---|---|---|
| Mol* | Every 3D view, including the three-way overlay | MIT | [github.com/molstar/molstar](https://github.com/molstar/molstar) | [Sehnal et al. 2021, Nucleic Acids Research](https://doi.org/10.1093/nar/gkab314) |
| Plotly.js | Trajectory panels | MIT | [github.com/plotly/plotly.js](https://github.com/plotly/plotly.js) | plotly.com |

## Hosting

| Tool | Role | Licence | Repository | Reference |
|---|---|---|---|---|
| Flask, gunicorn, nginx, SQLite | The web application and its store | BSD / MIT / public domain | [flask.palletsprojects.com](https://flask.palletsprojects.com) |  |
| pixi and uv | The bundle's environment, and the server's | BSD-3-Clause / Apache-2.0 | [github.com/prefix-dev/pixi](https://github.com/prefix-dev/pixi) |  |

