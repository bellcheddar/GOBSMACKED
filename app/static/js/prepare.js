/* Prepare: four panels, each unlocking the next, ending in a bundle.

   The page holds all the state; the server holds none until Generate bundle is
   pressed. Every panel posts what it has and renders what comes back, so a
   half-finished Prepare leaves nothing behind. */
(function () {
  "use strict";

  var state = {
    protein: null,        // the /api/fetch payload
    annotation: null,
    ligand: null,
    pocket: null,
    reference: null,      // {pdb_id, chain, ligand_ccd}
    references: [],
    visibility: "public",
    numbering: {},        // sequence position -> structure residue number
    inverse: {},          // structure residue number -> sequence position
  };

  var scope = null;
  var track = null;
  var selected = {};      // "chain:number" -> true, in STRUCTURE numbering

  function $(id) { return document.getElementById(id); }

  function unlock(id, open) {
    var panel = $(id);
    if (panel) panel.classList.toggle("locked", !open);
  }

  function say(id, html, kind) {
    var node = $(id);
    if (!node) return;
    node.innerHTML = html;
    node.style.color = kind === "error" ? "var(--red)" : "";
  }

  function postJSON(url, body) {
    return fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    }).then(function (r) {
      return r.json().then(function (data) {
        if (!r.ok) throw new Error(data.error || ("request failed (" + r.status + ")"));
        return data;
      });
    });
  }

  function setStage(index, stateName, text) {
    var cells = document.querySelectorAll(".stages .stage");
    if (!cells[index]) return;
    cells[index].className = "stage " + stateName;
    cells[index].querySelector("span").textContent = text;
  }

  // --- Panel 1: protein ----------------------------------------------------

  function fetchProtein() {
    var query = $("query").value.trim();
    var pdbHint = $("pdb-hint").value.trim();
    var file = $("upload").files[0];
    if (!query && !file) {
      say("protein-result", "Type a UniProt accession, PDB ID or sequence first.", "error");
      return;
    }
    say("protein-result", "Fetching...");
    setStage(0, "warn", "fetching");

    var form = new FormData();
    form.append("query", query);
    form.append("pdb_id", pdbHint);
    if (file) form.append("file", file);

    fetch("/api/fetch", { method: "POST", body: form })
      .then(function (r) { return r.json().then(function (d) {
        if (!r.ok) throw new Error(d.error || "fetch failed"); return d; }); })
      .then(function (data) {
        state.protein = data;
        state.numbering = data.numbering || {};
        state.inverse = {};
        Object.keys(state.numbering).forEach(function (pos) {
          state.inverse[state.numbering[pos]] = parseInt(pos, 10);
        });
        renderProtein(data);
        return annotate(data);
      })
      .catch(function (err) {
        say("protein-result", err.message, "error");
        setStage(0, "fail", "failed");
      });
  }

  function renderProtein(data) {
    var bits = [];
    if (data.protein_name) bits.push("<b>" + data.protein_name + "</b>");
    if (data.accession) bits.push(data.accession + (data.gene ? " (" + data.gene + ")" : ""));
    if (data.organism) bits.push(data.organism);
    bits.push(data.sequence.length + " residues");
    var source = {
      afdb: "AlphaFold DB model", esm_atlas: "ESM Atlas model", pdb: "PDB entry",
      user_pdb: "your upload", fold: "no model: the bundle will fold it",
    }[data.source_structure] || data.source_structure;
    bits.push(source + (data.source_id ? " " + data.source_id : ""));
    if (data.mean_plddt) bits.push("mean pLDDT " + data.mean_plddt);
    if (data.resolution) bits.push(data.resolution.toFixed(2) + " A");
    say("protein-result", bits.join(" &middot; "));

    var warnings = $("protein-warnings");
    warnings.innerHTML = "";
    (data.warnings || []).forEach(function (w) {
      var li = document.createElement("li");
      li.textContent = w;
      warnings.appendChild(li);
    });

    setStage(0, "done", source + (data.source_id ? " " + data.source_id : ""));
    setStage(2, data.source_structure === "fold" ? "warn" : "done",
             data.source_structure === "fold" ? "ESMFold in the bundle" : "skipped, model supplied");

    unlock("panel-ligand", true);
    unlock("panel-track", true);
    unlock("panel-pocket", true);

    if (data.structure_url) {
      $("scope-empty").style.display = "none";
      loadScope(data.structure_url);
    } else {
      $("scope-empty").innerHTML = "<span>No structure yet: the bundle will fold this sequence, " +
        "so pick the pocket by residue number below.</span>";
    }
  }

  function loadScope(url) {
    var hud = $("scope-hud");
    if (!scope) {
      GobViewer.create("viewer").then(function (created) {
        scope = created;
        scope.onResidueClick(function (hit) { toggleResidue(hit.chain, hit.seqId); });
        return scope.load("model", url, { primary: true });
      }).then(function () {
        hud.textContent = "input model";
      }).catch(function (err) {
        $("scope-empty").style.display = "";
        $("scope-empty").innerHTML = "<span>" + err.message + "</span>";
      });
    } else {
      scope.load("model", url, { primary: true }).then(function () {
        hud.textContent = "input model";
      });
    }
  }

  // --- Panel 2: annotation -------------------------------------------------

  function annotate(protein) {
    setStage(1, "warn", "annotating");
    return postJSON("/api/annotate", {
      uniprot: protein.accession, sequence: protein.sequence,
      gene: protein.gene, features: protein.features,
    }).then(function (data) {
      state.annotation = data;
      renderAnnotation(data, protein);
    }).catch(function (err) {
      setStage(1, "fail", err.message);
    });
  }

  function renderAnnotation(data, protein) {
    var family = data.family;
    setStage(1, "done", family === "other" ? "UniProt features only" :
             (family === "kinase" ? "KLIFS pocket" : "GPCRdb numbering"));
    setStage(6, "pending", family === "kinase" ? "DFG, aC, subpockets" :
             (family === "gpcr" ? "site and microswitches" : "contacts only"));

    var rows = [];
    rows.push(switchRow("Family", family, family !== "other"));
    (data.pfam || []).slice(0, 4).forEach(function (d) {
      rows.push(switchRow(d.pfam, d.name + " " + d.start + "-" + d.end, true));
    });
    if (data.klifs) {
      var named = data.klifs.named_positions || {};
      rows.push(switchRow("KLIFS", data.klifs.kinase_id ? "kinase " + data.klifs.kinase_id : "none",
                          !!data.klifs.kinase_id));
      rows.push(switchRow("Gatekeeper", named.gatekeeper || "unmapped", !!named.gatekeeper));
      rows.push(switchRow("DFG Asp", named.dfg_asp || "unmapped", !!named.dfg_asp));
      rows.push(switchRow("Hinge", (named.hinge || []).filter(Boolean).join(", ") || "unmapped",
                          (named.hinge || []).filter(Boolean).length > 0));
    }
    if (data.gpcrdb) {
      rows.push(switchRow("GPCRdb", data.gpcrdb.entry_name, true));
      Object.keys(data.gpcrdb.switches || {}).forEach(function (name) {
        var value = data.gpcrdb.switches[name];
        rows.push(switchRow(name.replace(/_/g, " "), value || "absent", !!value));
      });
    }
    $("family-switches").innerHTML = rows.join("");
    $("annotation-notes").innerHTML = (data.notes || []).join("<br>");

    track = new SequenceTrack("track", {
      onToggle: function () { syncFromTrack(); },
      onHover: function (position) {
        var structureNumber = state.numbering[position];
        if (scope && structureNumber) scope.highlightResidue(state.protein.chain, structureNumber);
      },
    });
    track.render({
      sequence: protein.sequence, pfam: data.pfam, features: data.features,
      positions: data.positions,
    });
  }

  /* "696-1022" -> [696, 1022]. Anything else is no range at all rather than a
     half-parsed one: trimming to the wrong span would be worse than not
     trimming, and silently so. */
  function parseRange(text) {
    var match = /^\s*(\d+)\s*[-\u2013:]\s*(\d+)\s*$/.exec(text || "");
    if (!match) return null;
    var lo = parseInt(match[1], 10);
    var hi = parseInt(match[2], 10);
    return hi > lo ? [lo, hi] : null;
  }

  /* The Pfam domain the pocket sits in, or the largest one when nothing is
     selected yet. This is the button most people want: an AlphaFold model is of
     the whole precursor, and solvating 883 residues around a 327-residue kinase
     domain costs an order of magnitude in MD time for nothing. */
  function rangeFromPfam() {
    var domains = (state.annotation && state.annotation.pfam) || [];
    if (!domains.length) {
      say("protein-result", "No Pfam domains fetched yet: press Fetch first.", "error");
      return;
    }
    var positions = Object.keys(selected).map(function (key) {
      return state.inverse[parseInt(key.split(":")[1], 10)];
    }).filter(Boolean);
    var chosen = null;
    if (positions.length) {
      var mid = positions[Math.floor(positions.length / 2)];
      chosen = domains.filter(function (d) { return d.start <= mid && mid <= d.end; })[0];
    }
    if (!chosen) {
      chosen = domains.slice().sort(function (a, b) {
        return (b.end - b.start) - (a.end - a.start);
      })[0];
    }
    $("residue-range").value = chosen.start + "-" + chosen.end;
    say("protein-result", "Trimming to " + chosen.pfam + " " + chosen.name +
        ", residues " + chosen.start + " to " + chosen.end + ".");
  }

  function switchRow(label, value, on) {
    return "<li><span><i class='led " + (on ? "" : "off") + "'></i>" + label + "</span><span>" +
           value + "</span></li>";
  }

  // --- Panel 3: ligand and pocket -----------------------------------------

  function checkLigand() {
    var smiles = $("smiles").value.trim();
    if (!smiles) { say("ligand-props", "Paste a SMILES first.", "error"); return; }
    postJSON("/api/ligand", { smiles: smiles }).then(function (data) {
      state.ligand = data;
      $("ligand-depiction").innerHTML = data.svg;
      say("ligand-props", [
        data.formula, data.mw + " Da", data.heavy_atoms + " heavy atoms",
        data.rotatable + " rotatable", "logP " + data.logp,
      ].join(" &middot; "));
      unlock("panel-reference", true);
      setStage(3, "done", "PandaDock " + $("dock-mode").value);
      setStage(4, "done", $("production").value + " ps OpenMM");
      maybeUnlockRun();
    }).catch(function (err) {
      say("ligand-props", err.message, "error");
    });
  }

  function toggleResidue(chain, number) {
    var key = chain + ":" + number;
    if (selected[key]) delete selected[key]; else selected[key] = true;
    syncToTrack();
    writeResidueBox();
    computeBox();
  }

  function syncToTrack() {
    if (!track) return;
    var positions = Object.keys(selected).map(function (key) {
      return state.inverse[parseInt(key.split(":")[1], 10)];
    }).filter(Boolean);
    track.setSelected(positions);
  }

  function syncFromTrack() {
    var chain = (state.protein && state.protein.chain) || "A";
    selected = {};
    track.residues().forEach(function (position) {
      var number = state.numbering[position];
      if (number !== undefined) selected[chain + ":" + number] = true;
    });
    writeResidueBox();
    computeBox();
  }

  function writeResidueBox() {
    $("residues").value = Object.keys(selected).sort(byResidue).join(", ");
    if (scope) {
      scope.showPocket(Object.keys(selected).map(function (key) {
        var parts = key.split(":");
        return { chain: parts[0], seqId: parseInt(parts[1], 10) };
      }));
    }
  }

  function byResidue(a, b) {
    return parseInt(a.split(":")[1], 10) - parseInt(b.split(":")[1], 10);
  }

  function readResidueBox() {
    selected = {};
    $("residues").value.split(",").forEach(function (item) {
      var text = item.trim();
      if (!text) return;
      selected[text.indexOf(":") >= 0 ? text : ((state.protein && state.protein.chain) || "A") + ":" + text] = true;
    });
    syncToTrack();
  }

  function computeBox(fromLigand) {
    if (!state.protein || !state.protein.structure_name) {
      say("pocket-result", "No structure to measure a box on. The bundle will centre the box on " +
          "the residues you name once it has folded the sequence.", "error");
      state.pocket = { method: "residues", residues: Object.keys(selected) };
      maybeUnlockRun();
      return;
    }
    var body = {
      structure_name: state.protein.structure_name,
      chain: state.protein.chain,
      residues: Object.keys(selected),
    };
    if (state.reference && state.reference.pdb_id && state.reference.ligand_ccd &&
        state.pocketFromReference) {
      body.size_from_reference = { pdb_id: state.reference.pdb_id,
                                   ligand_ccd: state.reference.ligand_ccd };
    }
    if (fromLigand) body.from_ligand = fromLigand;
    if (!body.residues.length && !fromLigand) {
      say("pocket-result", "");
      return;
    }
    postJSON("/api/pocket", body).then(function (data) {
      state.pocket = {
        method: fromLigand ? "reference_ligand" : "residues",
        residues: data.residues, center: data.center, box: data.box,
      };
      if (fromLigand) {
        selected = {};
        data.residues.forEach(function (r) { selected[r] = true; });
        syncToTrack();
        $("residues").value = data.residues.join(", ");
        if (scope) scope.showPocket(data.residues.map(function (key) {
          var parts = key.split(":");
          return { chain: parts[0], seqId: parseInt(parts[1], 10) };
        }));
      }
      say("pocket-result", data.n_residues + " residues &middot; centre " +
          data.center.map(function (v) { return v.toFixed(1); }).join(", ") +
          " &middot; box " + data.box.join(" x ") + " A" +
          (data.sized_from ? " (sized on " + data.sized_from + ")" : "") +
          (data.missing && data.missing.length ? "<br>Not in the structure: " + data.missing.join(", ") : ""));
      if (scope) scope.focus(data.residues.map(function (key) {
        var parts = key.split(":");
        return { chain: parts[0], seqId: parseInt(parts[1], 10) };
      }));
      maybeUnlockRun();
    }).catch(function (err) {
      say("pocket-result", err.message, "error");
    });
  }

  // --- Panel 4: reference --------------------------------------------------

  function searchReferences() {
    if (!state.protein || !state.protein.accession) {
      $("reference-list").innerHTML = "<p class='note'>Reference search needs a UniProt " +
        "accession. Type a PDB ID above instead, or continue unverified.</p>";
      return;
    }
    $("reference-list").innerHTML = "<p class='note'>Searching RCSB...</p>";
    postJSON("/api/references", {
      uniprot: state.protein.accession,
      smiles: (state.ligand && state.ligand.smiles) || $("smiles").value.trim(),
    }).then(function (data) {
      state.references = data.entries || [];
      renderReferences(data);
    }).catch(function (err) {
      $("reference-list").innerHTML = "<p class='note' style='color:var(--red)'>" + err.message + "</p>";
    });
  }

  function renderReferences(data) {
    if (!data.entries || !data.entries.length) {
      $("reference-list").innerHTML = "<p class='note'>" +
        (data.note || "RCSB has no entries for this accession.") + "</p>";
      return;
    }
    var html = "<div class='scroll' style='max-height:420px;overflow-y:auto'><table><tbody>";
    data.entries.forEach(function (entry) {
      var lig = entry.best_ligand || {};
      html += "<tr data-pdb='" + entry.pdb_id + "' style='cursor:pointer'>" +
        "<td><b>" + entry.pdb_id + "</b><br><span class='muted'>" +
        (entry.resolution ? entry.resolution.toFixed(2) + " A" : "no resolution") + "</span></td>" +
        "<td>" + (lig.ccd ? lig.ccd : "apo") +
        (entry.tanimoto !== null && entry.tanimoto !== undefined ?
          "<br><span class='muted'>T " + entry.tanimoto + "</span>" : "") + "</td>" +
        "<td style='width:150px'>" + (lig.svg || "") + "</td></tr>";
    });
    html += "</tbody></table></div>";
    $("reference-list").innerHTML = html;
    $("reference-list").querySelectorAll("tr").forEach(function (row) {
      row.addEventListener("click", function () { chooseReference(row.getAttribute("data-pdb")); });
    });
    if (data.default) chooseReference(data.default);
  }

  function chooseReference(pdbId) {
    var entry = state.references.filter(function (e) { return e.pdb_id === pdbId; })[0];
    var ccd = entry && entry.best_ligand ? entry.best_ligand.ccd : null;
    state.reference = { pdb_id: pdbId, ligand_ccd: ccd, chain: null };
    setStage(5, "done", pdbId + (ccd ? " " + ccd : " (apo)"));
    $("reference-list").querySelectorAll("tr").forEach(function (row) {
      row.style.background = row.getAttribute("data-pdb") === pdbId ? "var(--panel-2)" : "";
    });
    postJSON("/api/reference_site", {
      pdb_id: pdbId, ligand_ccd: ccd,
      structure_name: state.protein && state.protein.structure_name,
      chain: state.protein && state.protein.chain,
    }).then(function (data) {
      state.reference.chain = (data.chains && data.chains.length) ? data.chains[0].chain : null;
      state.reference.ligand_ccd = data.ligand_ccd;
      state.reference.site_residues = data.residues;
      state.reference.structure_url = data.structure_url;
      if (scope && data.structure_url) {
        scope.load("reference", data.structure_url, { color: GobViewer.COLOURS.reference });
      }
    }).catch(function () { /* the site is optional: the reference still counts */ });
  }

  function useReferenceSite() {
    if (!state.reference || !state.reference.ligand_ccd) {
      say("pocket-result", "Pick a reference entry with a ligand first.", "error");
      return;
    }
    // The site is defined on the reference structure, so the residues come back
    // in the reference's numbering: map them onto the model through the
    // sequence positions before they become the docking box.
    if (state.protein && state.protein.structure_name) {
      computeBox(null);
    }
    postJSON("/api/reference_site", {
      pdb_id: state.reference.pdb_id, ligand_ccd: state.reference.ligand_ccd,
      // Naming the model is what lets the server hand back the site in the
      // model's numbering rather than the crystal's.
      structure_name: state.protein && state.protein.structure_name,
      chain: state.protein && state.protein.chain,
    }).then(function (data) {
      state.pocketFromReference = true;
      $("residues").value = data.residues.join(", ");
      if (data.renumbered) {
        say("pocket-result", "Site taken from " + data.pdb_id + " and renumbered onto this model: "
            + data.reference_residues.length + " crystal residues, "
            + data.residues.length + " matched.");
      }
      readResidueBox();
      writeResidueBox();
      computeBox();
    }).catch(function (err) { say("pocket-result", err.message, "error"); });
  }

  // --- Generate ------------------------------------------------------------

  function maybeUnlockRun() {
    unlock("panel-run", !!(state.ligand && state.pocket));
  }

  function generate() {
    if (!state.protein || !state.ligand) {
      say("bundle-result", "Fetch a protein and check a ligand first.", "error");
      return;
    }
    var button = $("bundle-btn");
    button.disabled = true;
    button.textContent = "Writing bundle...";
    postJSON("/api/bundle", {
      title: $("title").value.trim(),
      visibility: state.visibility,
      protein: {
        uniprot: state.protein.accession,
        sequence: state.protein.sequence,
        source_structure: state.protein.source_structure,
        source_id: state.protein.source_id,
        chain: state.protein.chain,
        family: (state.annotation && state.annotation.family) || "other",
        residue_range: parseRange($("residue-range").value),
        protein_name: state.protein.protein_name,
        gene: state.protein.gene,
        structure_name: state.protein.structure_name,
      },
      ligand: {
        name: $("ligand-name").value.trim() || "ligand",
        smiles: state.ligand.smiles,
        protonation_ph: parseFloat($("ph").value) || 7.4,
      },
      pocket: state.pocket,
      reference: Object.assign({}, state.reference || {}, {
        apo_pdb_id: ($("reference-apo").value.trim().toUpperCase() || null),
      }),
      docking: { mode: $("dock-mode").value, num_poses: parseInt($("num-poses").value, 10) },
      md: {
        production_ps: parseInt($("production").value, 10),
        frame_interval_ps: parseInt($("frame-interval").value, 10),
      },
    }).then(function (data) {
      renderBundle(data);
    }).catch(function (err) {
      say("bundle-result", "<p class='note' style='color:var(--red)'>" + err.message + "</p>");
    }).then(function () {
      button.disabled = false;
      button.textContent = "Generate bundle";
    });
  }

  function renderBundle(data) {
    try {
      var keys = JSON.parse(localStorage.getItem("gobsmacked_keys") || "{}");
      keys[data.job_id] = data.owner_token;
      localStorage.setItem("gobsmacked_keys", JSON.stringify(keys));
    } catch (err) { /* private browsing: the key is still shown below */ }

    $("bundle-result").innerHTML =
      "<h3>" + data.job_id + "</h3>" +
      "<p><a class='btn primary' href='" + data.bundle_url + "'>Download " + data.bundle_name + "</a></p>" +
      "<pre class='cmd'>" + data.commands.join("\n") + "</pre>" +
      "<p class='note'>Owner key, shown once. It is stored in this browser and rides inside the " +
      "bundle, so uploading results needs no typing. Keep a copy to change visibility or delete " +
      "the run from elsewhere.</p>" +
      "<code class='key' id='owner-key'>" + data.owner_token + "</code>" +
      "<button class='small' id='copy-key'>Copy key</button> " +
      "<a class='btn small' href='" + data.run_url + "'>Open the run page</a>";

    var copy = $("copy-key");
    if (copy) copy.addEventListener("click", function () {
      navigator.clipboard.writeText(data.owner_token).then(function () {
        copy.textContent = "Copied";
      });
    });
  }

  // --- Wiring --------------------------------------------------------------

  /* The example is a real value in a muted colour, not a placeholder: pressing
     Fetch straight away runs it. The moment a field is edited it stops looking
     like a suggestion, and the banner goes once anything has been changed. */
  function wireExample() {
    var fields = Array.prototype.slice.call(document.querySelectorAll("[data-example]"));
    fields.forEach(function (field) {
      if (!field.value) field.value = field.getAttribute("data-example");
      var promote = function () {
        field.classList.remove("example");
        var banner = $("example-banner");
        if (banner) banner.classList.add("gone");
      };
      field.addEventListener("input", promote);
      field.addEventListener("paste", promote);
      // Focus alone does not count as editing: tabbing through to look at a
      // field should not claim the example as the user's own.
      field.addEventListener("keydown", function (event) {
        if (event.key !== "Tab" && event.key !== "Shift") promote();
      });
    });
    var clear = $("clear-example");
    if (clear) {
      clear.addEventListener("click", function () {
        fields.forEach(function (field) {
          field.value = "";
          field.classList.remove("example");
        });
        $("example-banner").classList.add("gone");
        $("query").focus();
      });
    }
  }

  document.addEventListener("DOMContentLoaded", function () {
    wireExample();
    $("fetch-btn").addEventListener("click", fetchProtein);
    $("range-from-pfam").addEventListener("click", rangeFromPfam);
    $("query").addEventListener("keydown", function (e) { if (e.key === "Enter") fetchProtein(); });
    $("ligand-btn").addEventListener("click", checkLigand);
    $("pocket-btn").addEventListener("click", function () { readResidueBox(); writeResidueBox(); computeBox(); });
    $("pocket-clear").addEventListener("click", function () {
      selected = {}; state.pocketFromReference = false;
      syncToTrack(); writeResidueBox(); say("pocket-result", "");
    });
    $("pocket-from-ref").addEventListener("click", useReferenceSite);
    $("reference-search").addEventListener("click", searchReferences);
    $("reference-none").addEventListener("click", function () {
      state.reference = null;
      setStage(5, "warn", "unverified");
      $("reference-list").innerHTML = "<p class='note'>No reference: the scorecard will show an " +
        "unverified banner and the overlay will have two states instead of three.</p>";
    });
    $("reference-manual-btn").addEventListener("click", function () {
      var id = $("reference-manual").value.trim().toUpperCase();
      if (id) chooseReference(id);
    });
    $("bundle-btn").addEventListener("click", generate);

    $("visibility").querySelectorAll("button").forEach(function (button) {
      button.addEventListener("click", function () {
        state.visibility = button.getAttribute("data-value");
        $("visibility").querySelectorAll("button").forEach(function (other) {
          other.setAttribute("aria-pressed", other === button ? "true" : "false");
        });
        $("visibility-note").textContent = state.visibility === "public"
          ? "This run appears in the Runs tab for anyone."
          : "Only people with the run link and the owner key can see it.";
      });
    });

    ["dock-mode", "production"].forEach(function (id) {
      $(id).addEventListener("change", function () {
        setStage(3, "done", "PandaDock " + $("dock-mode").value);
        setStage(4, "done", $("production").value + " ps OpenMM");
      });
    });
  });
})();
