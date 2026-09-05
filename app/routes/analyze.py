"""Analyze: take a results archive apart and say whether the prediction held.

The whole pipeline runs inside the upload request, in this order:

    ingest -> superpose -> interactions -> modes -> dynamics -> scorecard

It is deliberately synchronous. Every step is seconds on a CPU (PLIP on three
complexes, one Kabsch fit, one small trajectory read), the alternative is a job
queue and a polling UI for a wait shorter than most page loads, and a failure
can be reported in the response that caused it rather than found in a log.
"""

from __future__ import annotations

import json
import shutil
import time
from pathlib import Path
from typing import Any, Optional

from flask import (Blueprint, abort, jsonify, render_template, request,
                   url_for)

from .. import config, db
from ..services import annotate as annotate_svc
from ..services import dynamics as dyn_svc
from ..services import fetch as fetch_svc
from ..services import ingest as ingest_svc
from ..services import interactions as inter_svc
from ..services import modes as modes_svc
from ..services import scorecard as score_svc
from ..services import superpose as sup_svc

bp = Blueprint("analyze", __name__)

ANALYZE_STAGES = ["Fetch", "Annotate", "Fold", "Dock", "MD", "Verify", "Mode"]


# What each rail cell will be filled in from, before an archive arrives.
ANALYZE_WAITING = {
    "Fetch": "the structure the bundle started from",
    "Annotate": "the family, from the campaign",
    "Fold": "whether ESMFold ran or a model was supplied",
    "Dock": "poses and scores from the archive",
    "MD": "the trajectory summary",
    "Verify": "the reference named in the campaign",
    "Mode": "family labels, once the archive is read",
}


@bp.route("/analyze")
def analyze_page():
    stages = [{"name": n, "state": "pending", "text": ANALYZE_WAITING[n]}
              for n in ANALYZE_STAGES]
    return render_template(
        "analyze.html", tab="analyze", stages=stages,
        strip_actions='<span class="rail-note">Needed: a results.tar.gz from a run</span>')


@bp.post("/api/upload")
def api_upload():
    """Accept results.tar.gz, analyse it, and answer with the run's URL."""
    if "file" not in request.files or not request.files["file"].filename:
        return jsonify({"error": "No file in the upload."}), 400
    upload = request.files["file"]
    if not upload.filename.endswith((".tar.gz", ".tgz")):
        return jsonify({"error": "Upload the results.tar.gz the bundle wrote."}), 400

    staging = config.DATA_DIR / "incoming"
    staging.mkdir(parents=True, exist_ok=True)
    temp_archive = staging / f"upload_{int(time.time() * 1000)}.tar.gz"
    upload.save(temp_archive)

    try:
        # Peek at the campaign before choosing a directory: the job id lives
        # inside the archive, not in the filename.
        peek = ingest_svc.extract(temp_archive, staging / temp_archive.stem)
    except ingest_svc.IngestError as exc:
        temp_archive.unlink(missing_ok=True)
        return jsonify({"error": str(exc)}), 400

    job_id = peek.job_id
    row = db.get_job(job_id) if job_id else None
    if row is None:
        shutil.rmtree(peek.root, ignore_errors=True)
        temp_archive.unlink(missing_ok=True)
        return jsonify({"error": f"This archive belongs to job {job_id or '(unnamed)'}, which "
                                 f"this server has no record of. Bundles and results must come "
                                 f"from the same instance."}), 404
    if row["visibility"] == "private" and not db.token_matches(row, peek.owner_token):
        shutil.rmtree(peek.root, ignore_errors=True)
        temp_archive.unlink(missing_ok=True)
        return jsonify({"error": "This run is private and the archive carries no matching "
                                 "owner key."}), 403

    run_dir = config.RUNS_DIR / job_id
    run_dir.mkdir(parents=True, exist_ok=True)
    results_dir = run_dir / "results"
    if results_dir.exists():
        shutil.rmtree(results_dir, ignore_errors=True)
    shutil.move(str(peek.root), str(results_dir))
    shutil.move(str(temp_archive), str(run_dir / "results.tar.gz"))
    results = ingest_svc.Results(root=results_dir, manifest=peek.manifest,
                                 campaign=peek.campaign, summary=peek.summary,
                                 present=peek.present, missing_optional=peek.missing_optional,
                                 warnings=list(peek.warnings))
    if not ingest_svc.campaign_matches(results, row["campaign_yaml"] or ""):
        results.warnings.append(
            "The campaign in this archive differs from the one this server issued. "
            "The settings reported below are the ones actually run."
        )

    db.update_job(job_id, status="results_uploaded",
                  results_path=str(run_dir / "results.tar.gz"))
    try:
        card = analyse(row, results)
    except Exception as exc:                       # a bad archive must not 500
        db.update_job(job_id, status="failed", error=str(exc)[:500])
        return jsonify({"error": f"Analysis failed: {exc}"}), 500

    db.update_job(
        job_id, status="analysed", scorecard_json=json.dumps(card),
        gobsmack_score=card["scorecard"]["score"], grade=card["scorecard"]["grade"],
        mode_predicted=(card["modes"]["predicted"] or {}).get("label"),
        mode_reference=(card["modes"]["reference"] or {}).get("label"),
        mode_match=_as_int(card["modes"]["verdict"].get("match")),
        error=None,
    )
    return jsonify({"job_id": job_id, "url": url_for("runs.run_page", job_id=job_id),
                    "score": card["scorecard"]["score"], "grade": card["scorecard"]["grade"]})


def _as_int(value) -> Optional[int]:
    return None if value is None else int(bool(value))


# ---------------------------------------------------------------------------
# The pipeline
# ---------------------------------------------------------------------------

def analyse(row, results: ingest_svc.Results) -> dict[str, Any]:
    """Everything the results page shows, as one JSON-serialisable dict."""
    campaign = results.campaign
    protein = campaign.get("protein") or {}
    ligand = campaign.get("ligand") or {}
    pocket = campaign.get("pocket") or {}
    reference = campaign.get("reference") or {}
    sequence = protein.get("sequence", "")
    family = protein.get("family", "other")
    started = time.time()

    card: dict[str, Any] = {
        "job_id": results.job_id,
        "campaign": campaign,
        "manifest": results.manifest,
        "warnings": list(results.warnings),
        "reference": {"pdb_id": reference.get("pdb_id"),
                      "ligand_ccd": reference.get("ligand_ccd"),
                      "apo_pdb_id": reference.get("apo_pdb_id")},
    }

    # --- structures ------------------------------------------------------
    md_final = results.path("complex_md_final.pdb")
    pose1 = results.path("complex_pose1.pdb")
    minimised = results.path("complex_min.pdb")
    model_apo = results.path("model_apo.pdb")
    # Per file, not once. PandaDock names the ligand residue in the pose complex
    # and OpenMM names it in the relaxed one, and they need not agree: filtering
    # PLIP on the wrong name returns an empty site rather than an error, which
    # would read as "the docked pose makes no interactions at all".
    ligand_names = {state: _ligand_resname(path) for state, path in
                    (("pose1", pose1), ("minimised", minimised), ("md_final", md_final))
                    if path is not None}
    ligand_resname = ligand_names.get("md_final") or ligand_names.get("pose1")
    card["ligand_resname"] = ligand_resname
    card["ligand_resnames"] = ligand_names

    reference_pdb = None
    if reference.get("pdb_id"):
        fetched = fetch_svc.fetch_pdb(reference["pdb_id"])
        if fetched:
            reference_pdb = inter_svc.as_pdb(
                fetched["path"], results.root / "reference.pdb",
                chain=reference.get("chain") or None)
            card["reference"]["resolution"] = fetched.get("resolution")
        else:
            card["warnings"].append(
                f"RCSB has no entry {reference['pdb_id']}, so this run is unverified.")

    # The apo reference is optional and purely for the overlay: it is what the
    # pocket looked like with nothing bound, which is the shape a predicted
    # model tends to resemble.
    if reference.get("apo_pdb_id"):
        apo = fetch_svc.fetch_pdb(reference["apo_pdb_id"])
        if apo:
            inter_svc.as_pdb(apo["path"], results.root / "apo_reference.pdb",
                             chain=reference.get("apo_chain") or None, keep_ligand=False)
        else:
            card["warnings"].append(
                f"RCSB has no entry {reference['apo_pdb_id']}, so the overlay has no apo state.")

    # --- superposition ---------------------------------------------------
    geometry: dict[str, Any] = {}
    if reference_pdb:
        for label, path in (("pose1", pose1), ("md_final", md_final),
                            ("minimised", minimised), ("model", model_apo)):
            if path is None:
                continue
            got = sup_svc.compare(
                path, reference_pdb,
                reference_ccd=reference.get("ligand_ccd"),
                reference_chain=reference.get("chain") or None,
                model_ligand_ccd=ligand_names.get(label),
                ligand_smiles=ligand.get("smiles", ""),
            )
            geometry[label] = got
    card["geometry"] = geometry

    final = geometry.get("md_final") or {}
    first = geometry.get("pose1") or {}
    model_geo = geometry.get("model") or {}
    card["displacements"] = final.get("displacements") or []

    # --- put every state in one coordinate frame ---------------------------
    #
    # Until this point each file is correct in its own frame and none of them
    # agree: the crystal is in the crystal's frame, the AlphaFold model is in
    # AlphaFold's, and OpenMM has translated the complex into its periodic box.
    # Handing those to the viewer draws three structures scattered across the
    # scene. Everything below, including PLIP and the 2D map, runs on the
    # superposed copies, so the coordinates the page shows are the coordinates
    # the numbers were measured on.
    states = {"model": model_apo, "pose1": pose1, "minimised": minimised, "md_final": md_final}
    # Kept because the docking box is defined in the model's own frame: the
    # inside-the-box check has to be made there, not on a superposed copy.
    original_paths = dict(states)
    target = reference_pdb or model_apo
    target_chain = reference.get("chain") if reference_pdb else None
    pocket_target = final.get("pocket_residues_reference") if reference_pdb else None
    superposed: dict[str, Any] = {}
    aligned_paths: dict[str, Any] = {}
    for label, path in states.items():
        if path is None:
            continue
        if reference_pdb is None and label == "model":
            aligned_paths[label] = path      # everything else is aligned onto it
            superposed[label] = {"file": Path(path).name, "basis": "reference frame"}
            continue
        fit = sup_svc.align_onto(path, target, target_chain=target_chain,
                                 target_residues=pocket_target)
        if fit is None:
            aligned_paths[label] = path
            continue
        dest = results.root / f"superposed_{label}.pdb"
        sup_svc.write_transformed(path, dest, fit["rotation"], fit["translation"])
        aligned_paths[label] = dest
        superposed[label] = {"file": dest.name, "atoms": fit["atoms"],
                             "basis": fit["basis"], "rmsd": fit["rmsd"]}
    card["superposed"] = superposed
    card["superposed_onto"] = (card["reference"].get("pdb_id") if reference_pdb
                               else "the prepared model")

    # From here on, the aligned copies are the files everything reads.
    model_apo = aligned_paths.get("model", model_apo)
    pose1 = aligned_paths.get("pose1", pose1)
    minimised = aligned_paths.get("minimised", minimised)
    md_final = aligned_paths.get("md_final", md_final)
    ligand_names = {state: _ligand_resname(path) for state, path in
                    (("pose1", pose1), ("minimised", minimised), ("md_final", md_final))
                    if path is not None}
    ligand_resname = ligand_names.get("md_final") or ligand_names.get("pose1")
    card["ligand_resname"] = ligand_resname
    card["ligand_resnames"] = ligand_names

    # --- interactions ----------------------------------------------------
    fingerprints: dict[str, Any] = {}
    jaccards: dict[str, Optional[float]] = {}
    plip: dict[str, Any] = {}
    reference_fp: set = set()
    if reference_pdb:
        ref_result = inter_svc.run_plip(reference_pdb, reference.get("ligand_ccd"))
        plip["reference"] = ref_result
        reference_fp = inter_svc.fingerprint(ref_result)

    for label, path in (("pose1", pose1), ("md_final", md_final)):
        if path is None:
            continue
        result = inter_svc.run_plip(path, ligand_names.get(label))
        plip[label] = result
        mapping = {int(k): v for k, v in (geometry.get(label, {}).get("numbering_map") or {}).items()}
        fp = inter_svc.fingerprint(result, mapping or None)
        fingerprints[label] = sorted(f"{t}:{n}" for t, n in fp)
        jaccards[label] = inter_svc.jaccard(fp, reference_fp) if reference_fp else None

    card["plip"] = {k: {"counts": v.get("counts", {}), "n": v.get("n", 0),
                        "error": v.get("error")} for k, v in plip.items()}
    # Geometry for the 3D view, drawn from PLIP's own endpoints so the dashed
    # lines in the scope are the interactions the table lists rather than a
    # second opinion computed by the viewer.
    card["interaction_lines"] = {
        state: inter_svc.write_interaction_lines(plip[state], results.root, f"lines_{state}")
        for state in ("pose1", "md_final") if state in plip and not plip[state].get("error")
    }
    # The pocket residues as their own file, so the Complex view can draw them
    # as sticks over the cartoon.
    pocket_residues = pocket.get("residues") or []
    card["pocket_sticks"] = {}
    for state, path in (("pose1", pose1), ("minimised", minimised), ("md_final", md_final)):
        if path is None:
            continue
        written = inter_svc.write_pocket_sticks(path, pocket_residues,
                                                results.root / f"pocket_{state}.pdb")
        if written:
            card["pocket_sticks"][state] = written.name
    card["jaccard"] = jaccards
    if reference_fp:
        md_fp = {tuple(x.split(":", 1)) for x in fingerprints.get("md_final", [])}
        md_fp = {(t, int(n)) for t, n in md_fp}
        card["interaction_table"] = inter_svc.comparison_table(md_fp, reference_fp)
    else:
        card["interaction_table"] = []

    # --- 2D map ----------------------------------------------------------
    if md_final:
        card["pandamap"] = inter_svc.pandamap(md_final, results.root / "pandamap_md_final.png",
                                              ligand_names.get("md_final"))
    if reference_pdb:
        card["pandamap_reference"] = inter_svc.pandamap(
            reference_pdb, results.root / "pandamap_reference.png", reference.get("ligand_ccd"))

    # --- binding mode ----------------------------------------------------
    card["modes"] = classify_modes(family, protein, sequence, md_final, reference_pdb,
                                   reference, ligand_names.get("md_final"), plip)

    # --- dynamics --------------------------------------------------------
    dynamics = dyn_svc.summarise(results.summary)
    if reference_pdb and final.get("numbering_map") and results.path("traj/traj.dcd"):
        trace = dyn_svc.ligand_rmsd_to_reference(
            results.root / "traj" / "topology.pdb", results.root / "traj" / "traj.dcd",
            reference_pdb, final["numbering_map"], final.get("pocket_residues_reference") or [],
            reference_chain=reference.get("chain") or None,
            reference_ccd=reference.get("ligand_ccd"),
            ligand_resname=ligand_names.get("md_final"), smiles=ligand.get("smiles", ""))
        dynamics["ligand_rmsd_reference"] = trace
    card["dynamics"] = dynamics

    # --- scorecard -------------------------------------------------------
    ligand_rmsd = _best(final.get("ligand_rmsd"), first.get("ligand_rmsd"))
    validity_on = "md_final" if md_final else "pose1"
    validity_path = original_paths.get(validity_on)
    validity = score_svc.check_validity(
        validity_path, ligand.get("smiles", ""), pocket.get("center"), pocket.get("box"),
        ligand_ccd=ligand_names.get(validity_on)) if validity_path else {}
    card["scorecard"] = score_svc.composite({
        "ligand_rmsd": ligand_rmsd,
        "plip_jaccard": _best(jaccards.get("md_final"), jaccards.get("pose1"), lowest=False),
        "pocket_ca_rmsd": final.get("pocket_ca_rmsd"),
        "chi1_agreement": final.get("chi1_agreement"),
        "md_drift": dynamics.get("drift"),
        "rescue": dyn_svc.rescue(model_geo.get("pocket_ca_rmsd") or first.get("pocket_ca_rmsd"),
                                 final.get("pocket_ca_rmsd")),
    }, validity)
    card["scorecard"]["tm_score"] = final.get("tm_score")
    card["scorecard"]["pocket_ca_atoms"] = final.get("pocket_ca_atoms")
    card["timings"] = {"analysis_s": round(time.time() - started, 1),
                       **(results.manifest.get("timings") or {})}
    card["stages"] = build_stages(card, results)
    return card


def classify_modes(family: str, protein: dict, sequence: str, md_final, reference_pdb,
                   reference: dict, ligand_resname: Optional[str], plip: dict) -> dict[str, Any]:
    """Run the family's classifier on the prediction and on the crystal."""
    out: dict[str, Any] = {"family": family, "predicted": None, "reference": None}
    accession = protein.get("uniprot")
    md_plip = plip.get("md_final") or {}
    plip_rows = md_plip.get("interactions") if md_plip.get("interactions") else None

    if family == "kinase" and accession:
        klifs = annotate_svc.klifs_annotation(accession, sequence,
                                              protein.get("gene", "") or "")
        pocket_map = (klifs or {}).get("pocket_map") or {}
        if pocket_map:
            if md_final:
                out["predicted"] = modes_svc.classify_kinase(
                    md_final, sequence, pocket_map, ligand_ccd=ligand_resname,
                    plip_rows=plip_rows)
            if reference_pdb:
                out["reference"] = modes_svc.classify_kinase(
                    reference_pdb, sequence, pocket_map,
                    chain_name=reference.get("chain") or None,
                    ligand_ccd=reference.get("ligand_ccd"))
            out["klifs"] = {"kinase_id": (klifs or {}).get("kinase_id"),
                            "mapped": len(pocket_map),
                            "named": (klifs or {}).get("named_positions")}
        else:
            out["note"] = "KLIFS has no pocket for this accession, so no kinase labels."

    elif family == "gpcr" and accession:
        gp = annotate_svc.gpcrdb_annotation(accession)
        generic = (gp or {}).get("generic") or {}
        if generic:
            if md_final:
                out["predicted"] = modes_svc.classify_gpcr(
                    md_final, sequence, generic, ligand_ccd=ligand_resname)
            if reference_pdb:
                out["reference"] = modes_svc.classify_gpcr(
                    reference_pdb, sequence, generic,
                    chain_name=reference.get("chain") or None,
                    ligand_ccd=reference.get("ligand_ccd"))
                if md_final and out["predicted"] and not out["predicted"].get("error"):
                    out["npxxy_rmsd"] = modes_svc.npxxy_rmsd(
                        md_final, reference_pdb,
                        out["predicted"].get("positions") or {},
                        out["reference"].get("positions") or {},
                        reference_chain=reference.get("chain") or None)
            out["gpcrdb"] = {"entry_name": (gp or {}).get("entry_name"),
                             "numbered": len(generic)}
        else:
            out["note"] = "GPCRdb has no numbering for this accession, so no receptor labels."

    else:
        out["note"] = ("This target is neither a kinase nor a GPCR, so the mode panel lists "
                       "contacts against UniProt features rather than a family label.")

    out["verdict"] = modes_svc.compare_modes(out.get("predicted") or {}, out.get("reference"))
    return out


def build_stages(card: dict, results: ingest_svc.Results) -> list[dict]:
    """The stage strip for a finished run: what the archive contained."""
    manifest = results.manifest
    timings = manifest.get("timings") or {}
    campaign = results.campaign
    protein = campaign.get("protein") or {}
    folded = (results.root / "plddt.json").exists()

    def timing(key: str) -> str:
        v = timings.get(key)
        return f"{v:.0f} s" if isinstance(v, (int, float)) else ""

    stages = [
        {"name": "Fetch", "state": "ready",
         "text": f"{protein.get('source_structure', '?')} {protein.get('source_id') or ''}".strip()},
        {"name": "Annotate",
         "state": "ready" if card["modes"].get("family") != "other" else "optional",
         "text": (f"{card['modes'].get('family')}: "
                  + ("KLIFS pocket numbering" if card["modes"].get("family") == "kinase"
                     else "GPCRdb numbering" if card["modes"].get("family") == "gpcr"
                     else "UniProt features only"))},
        {"name": "Fold", "state": "ready",
         "text": (f"ESMFold, {timing('fold')}".strip() if folded
                  else "skipped, a model was supplied")},
        {"name": "Dock", "state": "ready",
         "text": f"{(campaign.get('docking') or {}).get('mode', 'hybrid')}"
                 f"{', ' + timing('dock') if timing('dock') else ''}"},
        {"name": "MD",
         "state": "ready" if card["dynamics"]["frames"] else "optional",
         "text": (f"{(campaign.get('md') or {}).get('production_ps', '?')} ps, "
                  f"{card['dynamics']['frames']} frames"
                  if card["dynamics"]["frames"] else "no trajectory in this archive")},
    ]
    verified = card["scorecard"].get("verified")
    fit = (card.get("superposed") or {}).get("md_final") or {}
    stages.append({
        "name": "Verify",
        "state": "ready" if verified else "optional",
        "text": (f"{card['reference'].get('pdb_id')}, {fit.get('atoms', '?')} pocket Ca"
                 + (f", fit {fit['rmsd']} A" if fit.get("rmsd") is not None else "")) if verified
                else "unverified: no reference was chosen",
    })
    verdict = card["modes"]["verdict"]
    stages.append({
        "name": "Mode",
        "state": ("ready" if verdict.get("match")
                  else "optional" if verdict.get("match") is False else "pending"),
        "text": (verdict.get("verdict", "") if verdict.get("match") is not None
                 else "unverified: nothing to compare the label to"),
    })
    return stages


def _ligand_resname(complex_path) -> Optional[str]:
    """The docked ligand's residue name, read from the complex the bundle wrote."""
    if complex_path is None:
        return None
    import gemmi

    from ..services.superpose import is_ligand_residue

    st = gemmi.read_structure(str(complex_path))
    st.setup_entities()
    st.remove_waters()
    best = None
    for ch in st[0]:
        for res in ch:
            if not is_ligand_residue(res):
                continue
            if best is None or len(res) > best[1]:
                best = (res.name, len(res))
    return best[0] if best else None


def _best(*values, lowest: bool = True) -> Optional[float]:
    present = [v for v in values if v is not None]
    if not present:
        return None
    return min(present) if lowest else max(present)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def render_results(row, card: Optional[dict], owned: bool = False):
    """The results page, shared by /runs/<job_id> and the post-upload redirect."""
    if card is None:
        stages = [{"name": n, "state": "pending", "text": ANALYZE_WAITING[n]}
                  for n in ANALYZE_STAGES]
        return render_template("run_waiting.html", tab="runs", row=row, stages=stages,
                               job_id=row["job_id"], owned=owned)
    return render_template("results.html", tab="runs", row=row, card=card,
                           card_json=json.dumps(card), stages=card.get("stages"),
                           job_id=row["job_id"], owned=owned, animate=True)
