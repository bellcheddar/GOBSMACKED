/* The results page: two Mol* scopes, the Plotly panels, and the owner controls. */
(function () {
  "use strict";

  var card = null;
  var jobId = document.body.getAttribute("data-job-id");
  var complexScope = null;
  var overlayScope = null;

  function fileUrl(name) { return "/runs/" + jobId + "/file/" + name; }

  /* Run `action` after Mol* has finished its own post-load camera work. */
  function settle(action) {
    requestAnimationFrame(function () { window.setTimeout(action, 400); });
  }

  function ownerToken() {
    try {
      var keys = JSON.parse(localStorage.getItem("gobsmacked_keys") || "{}");
      return keys[jobId] || "";
    } catch (err) { return ""; }
  }

  function startComplex() {
    var host = document.getElementById("complex-viewer");
    if (!host) return;
    GobViewer.create(host).then(function (scope) {
      complexScope = scope;
      return scope.load("md_final", fileUrl("complex_md_final.pdb"),
                        { primary: true, color: GobViewer.COLOURS.md_final });
    }).then(function () {
      // Mol*'s own preset resets the camera to the whole structure after the
      // load promise resolves, so focusing immediately is silently undone.
      // One frame plus a beat lands after it.
      settle(function () { complexScope.focusLigand("md_final"); });
      if (card.geometry && card.geometry.md_final) {
        var g = card.geometry.md_final;
        document.getElementById("complex-hud").textContent =
          "MD-final complex" + (g.tm_score ? " · TM " + g.tm_score : "");
      }
    }).catch(function (err) {
      document.getElementById("complex-hud").textContent = err.message;
    });

    document.querySelectorAll("#complex-toggles button").forEach(function (button) {
      button.addEventListener("click", function () {
        var state = button.getAttribute("data-state");
        var files = { pose1: "complex_pose1.pdb", minimised: "complex_min.pdb",
                      md_final: "complex_md_final.pdb" };
        document.querySelectorAll("#complex-toggles button").forEach(function (other) {
          other.setAttribute("aria-pressed", other === button ? "true" : "false");
        });
        if (!complexScope) return;
        complexScope.clearAll().then(function () {
          return complexScope.load(state, fileUrl(files[state]),
                                   { primary: true, color: GobViewer.COLOURS[state] || GobViewer.COLOURS.md_final });
        }).then(function () {
          settle(function () { complexScope.focusLigand(state); });
          document.getElementById("complex-hud").textContent = state.replace("_", " ");
        });
      });
    });
  }

  function startOverlay() {
    var host = document.getElementById("overlay-viewer");
    if (!host) return;
    GobViewer.create(host).then(function (scope) {
      overlayScope = scope;
      var work = scope.load("md_final", fileUrl("complex_md_final.pdb"),
                            { primary: true, color: GobViewer.COLOURS.md_final });
      if (card.reference && card.reference.pdb_id) {
        work = work.then(function () {
          return scope.load("reference", fileUrl("reference.pdb"),
                            { color: GobViewer.COLOURS.reference });
        });
      }
      return work;
    }).then(function () {
      var g = (card.geometry || {}).md_final || {};
      document.getElementById("overlay-hud").textContent =
        "pocket Ca superposition · " + (g.pocket_ca_atoms || 0) + " atoms" +
        (g.tm_score ? " · TM " + g.tm_score : "");
      // The reference is drawn in the crystal's own frame, so the overlay is
      // only meaningful because the bundle's structures were written in it too.
      if (g.pocket_residues_reference && overlayScope) {
        settle(function () {
          overlayScope.focus(g.pocket_residues_reference.map(function (n) {
            return { chain: g.reference_chain, seqId: n };
          }));
        });
      }
    }).catch(function (err) {
      document.getElementById("overlay-hud").textContent = err.message;
    });

    document.querySelectorAll("#overlay-toggles button").forEach(function (button) {
      button.addEventListener("click", function () {
        var state = button.getAttribute("data-state");
        var pressed = button.getAttribute("aria-pressed") === "true";
        button.setAttribute("aria-pressed", pressed ? "false" : "true");
        if (!overlayScope) return;
        if (overlayScope.has(state)) {
          overlayScope.setVisible(state, !pressed);
        } else if (!pressed) {
          var files = { model: "model_apo.pdb", md_final: "complex_md_final.pdb",
                        reference: "reference.pdb" };
          overlayScope.load(state, fileUrl(files[state]), { color: GobViewer.COLOURS[state] });
        }
      });
    });
  }

  function ownerControls() {
    var toggle = document.getElementById("toggle-visibility");
    var remove = document.getElementById("delete-run");
    var message = document.getElementById("owner-message");
    if (toggle) {
      toggle.addEventListener("click", function () {
        var current = toggle.getAttribute("data-visibility");
        var wanted = current === "public" ? "private" : "public";
        fetch("/api/runs/" + jobId + "/visibility", {
          method: "PATCH",
          headers: { "Content-Type": "application/json", "X-Owner-Token": ownerToken() },
          body: JSON.stringify({ visibility: wanted, token: ownerToken() }),
        }).then(function (r) { return r.json(); }).then(function (data) {
          if (data.error) { message.textContent = data.error; return; }
          toggle.setAttribute("data-visibility", wanted);
          toggle.textContent = "Make " + (wanted === "public" ? "private" : "public");
          message.textContent = "This run is " + wanted + ".";
        });
      });
    }
    if (remove) {
      remove.addEventListener("click", function () {
        var typed = window.prompt("Retype the job ID to delete this run and its files:");
        if (!typed) return;
        fetch("/api/runs/" + jobId, {
          method: "DELETE",
          headers: { "Content-Type": "application/json", "X-Owner-Token": ownerToken() },
          body: JSON.stringify({ confirm: typed, token: ownerToken() }),
        }).then(function (r) { return r.json(); }).then(function (data) {
          if (data.error) { message.textContent = data.error; return; }
          window.location = "/runs";
        });
      });
    }
  }

  document.addEventListener("DOMContentLoaded", function () {
    var node = document.getElementById("run-data");
    if (!node) return;
    try { card = JSON.parse(node.textContent); } catch (err) { return; }
    GobPlots.all(card.dynamics);
    startComplex();
    startOverlay();
    ownerControls();

    // The sub-nav marks the section the reader is in rather than the last one
    // clicked, so scrolling and clicking agree.
    var sections = ["scorecard", "complex", "overlay", "dynamics", "mode", "report"];
    window.addEventListener("scroll", function () {
      var current = sections[0];
      sections.forEach(function (id) {
        var el = document.getElementById(id);
        if (el && el.getBoundingClientRect().top < 120) current = id;
      });
      document.querySelectorAll(".subnav a").forEach(function (a) {
        a.classList.toggle("current", a.getAttribute("href") === "#" + current);
      });
    }, { passive: true });
  });
})();
