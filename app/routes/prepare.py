"""Prepare: four panels, each unlocking the next, ending in a run bundle.

The page holds no server-side session. Each panel posts what it has and gets
back what it needs, and only `/api/bundle` writes anything: one `jobs` row, one
campaign file and one archive. That keeps a half-finished Prepare from leaving
rows behind, and makes every endpoint independently testable.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from flask import (Blueprint, jsonify, render_template, request,
                   send_file, url_for)
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem.Draw import rdMolDraw2D

from .. import config, db
from ..services import annotate as annotate_svc
from ..services import bundle as bundle_svc
from ..services import fetch as fetch_svc
from ..services import references as ref_svc

bp = Blueprint("prepare", __name__)

SAFE_NAME = re.compile(r"^[A-Za-z0-9._-]+$")

PREPARE_STAGES = [
    ("Fetch", "waiting for an input"),
    ("Annotate", "family unknown"),
    ("Fold", "decided once a structure is chosen"),
    ("Dock", "PandaDock hybrid"),
    ("MD", "OpenMM, 1 ns"),
    ("Verify", "pick a reference"),
    ("Mode", "set by the family"),
]


@bp.route("/prepare")
def prepare_page():
    stages = [{"name": n, "text": t, "state": "pending"} for n, t in PREPARE_STAGES]
    return render_template("prepare.html", tab="prepare", stages=stages)


# ---------------------------------------------------------------------------
# Panel 1: protein
# ---------------------------------------------------------------------------

@bp.post("/api/fetch")
def api_fetch():
    """Resolve whatever the user typed or dropped into a sequence and a model."""
    uploaded_path = None
    if "file" in request.files and request.files["file"].filename:
        f = request.files["file"]
        suffix = Path(f.filename).suffix.lower()
        if suffix not in (".pdb", ".cif", ".ent", ".mmcif"):
            return jsonify({"error": "Upload a .pdb or .cif file."}), 400
        uploaded_path = config.STRUCT_CACHE / f"upload_{db.new_job_id()}{suffix}"
        uploaded_path.parent.mkdir(parents=True, exist_ok=True)
        f.save(uploaded_path)

    # The panel posts multipart when a file rides along and JSON otherwise, so
    # both shapes are read rather than assuming one.
    payload = request.get_json(silent=True) or {}
    text = request.form.get("query") or payload.get("query") or ""
    pdb_id = request.form.get("pdb_id") or payload.get("pdb_id") or ""

    try:
        result = fetch_svc.resolve_protein(text, pdb_id=pdb_id, uploaded=uploaded_path,
                                           job_hint=db.new_job_id())
    except Exception as exc:      # a public API being down is not a 500 for the user
        return jsonify({"error": f"Could not resolve that input: {exc}"}), 502

    out = result.to_dict()
    if result.structure_path:
        out["structure_url"] = url_for("prepare.api_structure",
                                       name=Path(result.structure_path).name)
        out["chains"] = fetch_svc.structure_chains(result.structure_path)
        out["ligands"] = fetch_svc.structure_ligands(result.structure_path)
        out["structure_name"] = Path(result.structure_path).name
        # The sequence track counts from 1; a crystal numbers however its
        # depositor chose (1M17 runs 24 low against UniProt, from the signal
        # peptide). Clicking residue 790 on the track has to select residue 766
        # in that structure, so the mapping is computed here rather than assumed
        # in the browser.
        out["numbering"] = _numbering(result.structure_path, result.chain, result.sequence)
    return jsonify(out)


def _translate_residues(reference_path: str, residues: list[str], model_path,
                        model_chain: str | None) -> list[str]:
    """Reference residue labels ("A:766") in the model's own numbering.

    Both structures are read as chains and aligned by sequence, which is the
    same route every measurement on the Analyze side takes. A residue with no
    counterpart is dropped rather than carried across unchanged.
    """
    from ..services import superpose as sup_svc

    reference = sup_svc.load_chain(reference_path)
    model = sup_svc.load_chain(model_path, model_chain)
    if reference is None or model is None:
        return []
    mapping = sup_svc.align_numbering(model, reference)          # model -> reference
    inverse = {ref: mod for mod, ref in mapping.items()}
    chain_name = model.name
    out = []
    for label in residues:
        try:
            number = int(str(label).split(":")[-1])
        except ValueError:
            continue
        if number in inverse:
            out.append(f"{chain_name}:{inverse[number]}")
    return out


def _numbering(structure_path: str, chain: str, sequence: str) -> dict[str, int]:
    """Sequence position (1-based, as a string key) -> residue number in the model."""
    from ..services import modes as modes_svc
    from ..services import superpose as sup_svc

    loaded = sup_svc.load_chain(structure_path, chain)
    if loaded is None or not sequence:
        return {}
    return {str(k): v for k, v in modes_svc.map_sequence_to_structure(loaded, sequence).items()}


@bp.get("/api/structure/<name>")
def api_structure(name: str):
    """Serve a cached structure to Mol*. Names only: no paths, no traversal."""
    if not SAFE_NAME.match(name):
        return jsonify({"error": "bad name"}), 400
    path = config.STRUCT_CACHE / name
    if not path.exists():
        return jsonify({"error": "not cached"}), 404
    return send_file(path, mimetype="chemical/x-pdb" if path.suffix == ".pdb" else "chemical/x-cif")


# ---------------------------------------------------------------------------
# Panel 2: annotation
# ---------------------------------------------------------------------------

@bp.post("/api/annotate")
def api_annotate():
    payload = request.get_json(silent=True) or {}
    accession = payload.get("uniprot")
    sequence = payload.get("sequence", "")
    gene = payload.get("gene", "")
    features = payload.get("features") or []
    if not sequence:
        return jsonify({"error": "No sequence to annotate."}), 400
    try:
        return jsonify(annotate_svc.annotate(accession, sequence, gene, features))
    except Exception as exc:
        return jsonify({"error": f"Annotation failed: {exc}"}), 502


# ---------------------------------------------------------------------------
# Panel 3: ligand and pocket
# ---------------------------------------------------------------------------

@bp.post("/api/ligand")
def api_ligand():
    """Validate a SMILES and draw it. The depiction is server-side RDKit, not a
    CDN library: the same code draws the reference ligands, so two structures
    that look alike here really are alike."""
    payload = request.get_json(silent=True) or {}
    smiles = (payload.get("smiles") or "").strip()
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return jsonify({"error": "RDKit cannot parse that SMILES."}), 400
    return jsonify({
        "smiles": Chem.MolToSmiles(mol),
        "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(Descriptors.MolWt(mol), 1),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "rotatable": Descriptors.NumRotatableBonds(mol),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "svg": depict(mol),
    })


def depict(mol, width: int = 260, height: int = 150) -> str:
    """A 2D depiction in the panel palette, as inline SVG."""
    mol = Chem.Mol(mol)
    Chem.rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.bondLineWidth = 1.4
    # Atom colours that read on a slate panel: the default palette is tuned for
    # white paper and disappears here.
    opts.updateAtomPalette({
        6: (0.91, 0.93, 0.96), 7: (0.36, 0.88, 0.90), 8: (1.0, 0.36, 0.36),
        16: (1.0, 0.71, 0.33), 17: (0.49, 0.89, 0.66), 9: (0.49, 0.89, 0.66),
        35: (0.76, 0.61, 1.0), 53: (0.76, 0.61, 1.0), 15: (1.0, 0.71, 0.33),
    })
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


@bp.post("/api/pocket")
def api_pocket():
    """Centre and box from a residue list, or from a reference ligand's site."""
    payload = request.get_json(silent=True) or {}
    name = payload.get("structure_name", "")
    if not SAFE_NAME.match(name or ""):
        return jsonify({"error": "No structure to measure the pocket on."}), 400
    path = config.STRUCT_CACHE / name
    if not path.exists():
        return jsonify({"error": "That structure is no longer cached; fetch it again."}), 404

    residues = payload.get("residues") or []
    if payload.get("from_ligand"):
        residues = bundle_svc.residues_near_ligand(path, payload["from_ligand"])
        if not residues:
            return jsonify({"error": f"No residues within "
                                     f"{config.POCKET_RADIUS:.0f} A of {payload['from_ligand']}."}), 400
    box = bundle_svc.pocket_box(path, residues, payload.get("chain", "A"))
    if not box:
        return jsonify({"error": "None of those residues are in the structure."}), 400

    # When the site came from a reference ligand, take the box SIZE from that
    # ligand and keep the centre from the residues. The centre has to come from
    # the model, because the crystal is in its own coordinate frame; the extent
    # is frame-independent and is the only honest size for the search.
    sized_from = payload.get("size_from_reference") or {}
    if sized_from.get("pdb_id") and sized_from.get("ligand_ccd"):
        fetched = fetch_svc.fetch_pdb(sized_from["pdb_id"])
        if fetched:
            extent = bundle_svc.ligand_extent(fetched["path"], sized_from["ligand_ccd"])
            if extent:
                box["box"] = extent
                box["sized_from"] = f"{sized_from['pdb_id']} {sized_from['ligand_ccd']}"
    return jsonify({"residues": residues, **box})


# ---------------------------------------------------------------------------
# Panel 4: references
# ---------------------------------------------------------------------------

@bp.post("/api/references")
def api_references():
    payload = request.get_json(silent=True) or {}
    accession = payload.get("uniprot")
    smiles = payload.get("smiles", "")
    if not accession:
        return jsonify({"entries": [], "default": None,
                        "note": "Reference search needs a UniProt accession. "
                                "Type a PDB ID instead, or continue unverified."})
    try:
        ranked = ref_svc.rank_references(accession, smiles)
    except Exception as exc:
        return jsonify({"error": f"RCSB search failed: {exc}"}), 502
    default = ref_svc.default_reference(ranked)
    for entry in ranked[:60]:
        lig = entry.get("best_ligand")
        if lig and lig.get("smiles"):
            mol = Chem.MolFromSmiles(lig["smiles"])
            lig["svg"] = depict(mol, 150, 90) if mol else ""
    return jsonify({
        "entries": ranked[:60],
        "default": (default or {}).get("pdb_id"),
        "total": len(ranked),
    })


@bp.post("/api/reference_site")
def api_reference_site():
    """Fetch a reference entry and return its ligand's contact residues.

    The residues come back in the reference's own numbering, which is what the
    pocket picker needs when the model and the crystal share a numbering (they
    do when the model is AFDB or the sequence is the canonical isoform).
    """
    payload = request.get_json(silent=True) or {}
    pdb_id = (payload.get("pdb_id") or "").strip().upper()
    ccd = (payload.get("ligand_ccd") or "").strip().upper()
    got = fetch_svc.fetch_pdb(pdb_id)
    if not got:
        return jsonify({"error": f"RCSB has no entry {pdb_id}."}), 404
    ligands = fetch_svc.structure_ligands(got["path"])
    if not ccd and ligands:
        ccd = ligands[0]["ccd"]
    reference_residues = bundle_svc.residues_near_ligand(got["path"], ccd) if ccd else []

    # Translated into the MODEL's numbering when the caller says which model it
    # is holding. 1M17 numbers EGFR from the mature protein and UniProt from the
    # precursor, 24 apart, so pasting the crystal's residue list straight into
    # the pocket picker selects real residues that are the wrong ones: a
    # silently misplaced docking box rather than an error.
    residues = reference_residues
    mapped = None
    model_name = payload.get("structure_name")
    if reference_residues and model_name and SAFE_NAME.match(model_name):
        model_path = config.STRUCT_CACHE / model_name
        if model_path.exists():
            mapped = _translate_residues(got["path"], reference_residues,
                                         model_path, payload.get("chain"))
            if mapped:
                residues = mapped

    return jsonify({
        "pdb_id": pdb_id,
        "ligand_ccd": ccd,
        "residues": residues,
        "reference_residues": reference_residues,
        "renumbered": bool(mapped) and mapped != reference_residues,
        "ligands": ligands,
        "chains": fetch_svc.structure_chains(got["path"]),
        "resolution": got.get("resolution"),
        "structure_url": url_for("prepare.api_structure", name=Path(got["path"]).name),
        "structure_name": Path(got["path"]).name,
    })


# ---------------------------------------------------------------------------
# Generate the bundle
# ---------------------------------------------------------------------------

@bp.post("/api/bundle")
def api_bundle():
    payload = request.get_json(silent=True) or {}
    protein = payload.get("protein") or {}
    ligand = payload.get("ligand") or {}
    pocket = payload.get("pocket") or {}
    reference = payload.get("reference") or {}

    if not protein.get("sequence"):
        return jsonify({"error": "No protein sequence: complete Panel 1 first."}), 400
    if not ligand.get("smiles"):
        return jsonify({"error": "No ligand SMILES: complete Panel 3 first."}), 400
    if not pocket.get("center"):
        return jsonify({"error": "No pocket: select residues or use a reference ligand's site."}), 400

    job_id = db.new_job_id()
    owner_token = db.new_owner_token()
    visibility = "private" if payload.get("visibility") == "private" else "public"
    title = (payload.get("title") or "").strip()[:120]

    campaign = bundle_svc.build_campaign(
        job_id=job_id, protein=protein, ligand=ligand, pocket=pocket,
        reference=reference, docking=payload.get("docking") or {},
        md=payload.get("md") or {}, owner_token=owner_token, title=title,
    )

    structure_name = protein.get("structure_name")
    structure_path = None
    if structure_name and SAFE_NAME.match(structure_name):
        candidate = config.STRUCT_CACHE / structure_name
        if candidate.exists():
            structure_path = str(candidate)

    try:
        archive = bundle_svc.write_bundle(job_id, campaign, structure_path)
    except Exception as exc:
        return jsonify({"error": f"Could not write the bundle: {exc}"}), 500

    db.insert_job(
        job_id=job_id,
        title=title or f"{protein.get('uniprot') or 'sequence'} + {ligand.get('name') or 'ligand'}",
        uniprot=protein.get("uniprot"),
        protein_name=protein.get("protein_name", ""),
        ligand_name=ligand.get("name", ""),
        ligand_smiles=ligand.get("smiles", ""),
        family=protein.get("family", "other"),
        reference_pdb=reference.get("pdb_id"),
        status="prepared",
        visibility=visibility,
        owner_hash=db.hash_token(owner_token),
        campaign_yaml=bundle_svc.dump_campaign(campaign),
    )

    return jsonify({
        "job_id": job_id,
        "owner_token": owner_token,
        "visibility": visibility,
        "bundle_url": url_for("runs.download_bundle", job_id=job_id),
        "bundle_name": archive.name,
        "run_url": url_for("runs.run_page", job_id=job_id),
        # `pixi run` solves and installs the environment it needs, so a separate
        # install step is not the first thing to tell anyone to do. It is still
        # in the bundle's README for anyone who wants to fetch the environment
        # before going offline.
        "commands": [
            f"tar xzf {archive.name} && cd run_bundle_{job_id}",
            "pixi run gobsmacked",
        ],
    })
