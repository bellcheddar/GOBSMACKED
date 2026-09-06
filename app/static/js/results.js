/* The results page: two Mol* scopes, the Plotly panels, and the owner controls. */
(function () {
  "use strict";

  var card = null;
  var jobId = document.body.getAttribute("data-job-id");
  var complexScope = null;
  var overlayScope = null;

  function fileUrl(name) { return "/runs/" + jobId + "/file/" + name; }

  /* The superposed copy of a state, which is what every view should draw: the
     raw files are each correct in their own frame and none of them agree. */
  function stateFile(state) {
    var entry = (card.superposed || {})[state];
    if (entry && entry.file) return entry.file;
    return { model: "model_apo.pdb", pose1: "complex_pose1.pdb",
             minimised: "complex_min.pdb", md_final: "complex_md_final.pdb",
             reference: "reference.pdb", apo: "apo_reference.pdb" }[state];
  }

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
      return scope.load("md_final", fileUrl(stateFile("md_final")),
                        { primary: true, color: GobViewer.COLOURS.md_final });
    }).then(function () {
      // Mol*'s own preset resets the camera to the whole structure after the
      // load promise resolves, so focusing immediately is silently undone.
      // One frame plus a beat lands after it.
      drawInteractions("md_final").then(function () {
        settle(function () { complexScope.focusLigand("md_final"); });
      });
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
        document.querySelectorAll("#complex-toggles button").forEach(function (other) {
          other.setAttribute("aria-pressed", other === button ? "true" : "false");
        });
        if (!complexScope) return;
        complexScope.clearAll().then(function () {
          return complexScope.load(state, fileUrl(stateFile(state)),
                                   { primary: true, color: GobViewer.COLOURS[state] || GobViewer.COLOURS.md_final });
        }).then(function () {
          drawInteractions(state).then(function () {
            settle(function () { complexScope.focusLigand(state); });
          });
          document.getElementById("complex-hud").textContent = state.replace("_", " ");
        });
      });
    });
  }

  /* The PLIP layer for one state, plus the legend that explains its colours.
     Minimised has no PLIP run of its own, so the layer is simply absent there
     rather than showing another state's interactions. */
  function drawInteractions(state) {
    // Always returns a promise: the caller waits on it before framing the
    // camera, and an early bare `return` here left it framing nothing and put
    // "Cannot read properties of undefined" in the HUD.
    if (!complexScope) return Promise.resolve();
    var work = drawPocket(state);
    var lines = (card.interaction_lines || {})[state] || [];
    var legend = document.getElementById("interaction-legend");
    if (!lines.length) {
      if (legend) legend.innerHTML = state === "minimised"
        ? "<span class='muted'>PLIP is run on the docked pose and the relaxed complex, not on the minimised intermediate.</span>"
        : "";
      return work.then(function () { return complexScope.hideInteractions(); });
    }
    work = work.then(function () { return complexScope.showInteractions(lines, fileUrl); });
    if (legend) {
      legend.innerHTML = lines.map(function (entry) {
        var hex = "#" + ("000000" + entry.colour.toString(16)).slice(-6);
        return "<span style='margin-right:14px'><i style='display:inline-block;width:18px;" +
               "height:2px;background:" + hex + ";vertical-align:middle;margin-right:6px'></i>" +
               entry.type + " (" + entry.count + ")</span>";
      }).join("");
    }
    return work;
  }

  /* The pocket residues as sticks over the cartoon. */
  function drawPocket(state) {
    if (!complexScope) return Promise.resolve();
    var file = (card.pocket_sticks || {})[state];
    return complexScope.remove("pocket").then(function () {
      if (!file) return;
      // ligandColor as well as color: the pocket file holds disconnected
      // residues, so Mol* classifies them as a ligand component and would
      // otherwise paint them the ligand's pink, which is the one colour on
      // screen that has to mean the ligand.
      return complexScope.load("pocket", fileUrl(file),
                               { color: GobViewer.COLOURS.model,
                                 ligandColor: GobViewer.COLOURS.model,
                                 format: "pdb", representation: "ball-and-stick" });
    });
  }

  /* Every pose, overlaid on the receptor and the crystal ligand.

     Three loads rather than one per pose: the top pose in phosphor because it
     is the one that went to MD, the rest in grey, and the crystal ligand in
     amber. Ten separate colours would be a rainbow, and the distinction the eye
     actually needs is "the pose that went forward" against "the ones that did
     not". */
  var posesScope = null;

  /* The trajectory itself, in Mol*, with controls.

     loadTrajectory rather than loadStructureFromData: the viewer build exports
     it, and it takes a topology and a coordinates file, which is exactly what
     the archive carries. Mol*'s own UI is left on for this one panel, because
     playing a trajectory needs a play button and the site's chrome has none. */
  var motionScope = null;

  function startMotion() {
    var host = document.getElementById("motion-viewer");
    if (!host || !(card.dynamics || {}).frames) return;
    GobViewer.create(host).then(function (scope) {
      return scope.viewer.loadTrajectory({
        model: { kind: "model-url", url: fileUrl("traj/topology.pdb"), format: "pdb" },
        coordinates: { kind: "coordinates-url", url: fileUrl("traj/traj.dcd"),
                       format: "dcd", isBinary: true },
        preset: "default",
      }).then(function () {
        // loadTrajectory goes round Scope.load, so nothing themed it, nothing
        // registered it and the camera never moved: the run rendered in Mol*'s
        // default green, half out of frame. Adopt the entry by hand.
        var all = scope.plugin.managers.structure.hierarchy.current.structures;
        var entry = all[all.length - 1];
        if (!entry) return scope;
        scope.entries.motion = entry;
        scope.primary = entry;
        return scope.theme(entry, GobViewer.COLOURS.md_final).then(function () {
          // After the canvas has its final size, not before: this card stretches
          // to its row, so a reset fired on load frames a canvas that is about
          // to change height and leaves the protein sitting on the bottom edge.
          settle(function () {
            try { scope.viewer.plugin.handleResize(); } catch (err) { /* keep the view */ }
          });
          return scope;
        });
      });
    }).then(function (scope) {
      motionScope = scope;
      document.getElementById("motion-hud").textContent =
        card.dynamics.frames + " frames \u00b7 drag to turn, scroll to zoom";
      // Mol*'s own model-index animation, driven from our button rather than
      // from its control panel: the panel is most of a half-width card and
      // nearly all of it is tools this page has no use for.
      var play = document.getElementById("motion-play");
      if (!play) return;

      function setPlaying(wanted) {
        var manager = scope.plugin.managers.animation;
        try {
          if (!wanted) {
            manager.stop();
          } else {
            var animation = (manager.animations || []).filter(function (a) {
              return a.name === "built-in.animate-model-index";
            })[0];
            if (!animation) throw new Error("no model-index animation in this build");
            manager.play(animation, { mode: { name: "loop", params: { direction: "forward" } },
                                      maxFPS: 15 });
          }
          play.setAttribute("aria-pressed", wanted ? "true" : "false");
          play.textContent = wanted ? "Pause" : "Play";
          return true;
        } catch (err) {
          document.getElementById("motion-hud").textContent =
            "playback unavailable: " + err.message;
          return false;
        }
      }

      play.addEventListener("click", function () {
        setPlaying(play.getAttribute("aria-pressed") !== "true");
      });
      // Running on arrival: a trajectory panel that opens as a still picture
      // looks like a viewer that failed to load, and the whole reason this
      // replaced a video is that it moves.
      //
      // The framing comes AFTER playback starts, and is repeated once. Stepping
      // the model index rebuilds the structure, and that re-frames the camera:
      // focusing first put the ligand in view and the first animation frame
      // pulled straight back out to the whole 253-residue chain.
      settle(function () {
        setPlaying(true);
        // The camera is left where Mol* puts it. Framing the ligand here does
        // not work and I could not find out why: focusLoci is called with a
        // non-empty element loci, reports success, and the view does not move,
        // whether it runs before playback, after it, once, or five times with
        // the structure re-resolved each attempt in case the model-index
        // animation had invalidated it. Rather than ship a retry loop that
        // does nothing, this records what was tried. Drag and scroll work, so
        // the panel is usable; it just opens on the whole chain.
      });
    }).catch(function (err) {
      document.getElementById("motion-hud").textContent =
        "the trajectory could not be loaded: " + err.message;
    });
  }

  function startPoses() {
    var host = document.getElementById("poses-viewer");
    var dp = card.poses || {};
    if (!host || !dp.rows || !dp.rows.length) return;
    GobViewer.create(host).then(function (scope) {
      posesScope = scope;
      var work = scope.load("protein", fileUrl(stateFile("pose1")),
                            { primary: true, color: GobViewer.COLOURS.model });
      if (dp.rest_file) {
        work = work.then(function () {
          // ligandColor, not color: these files are nothing but ligand, and
          // `color` themes the protein component, so the poses came out in
          // default element colours -- a red and blue tangle with no way to
          // tell the chosen pose from the rest.
          return scope.load("rest", fileUrl(dp.rest_file),
                            { color: 0x7f8fa6, ligandColor: 0x7f8fa6,
                              representation: "ball-and-stick" });
        });
      }
      work = work.then(function () {
        return scope.load("top", fileUrl(dp.top_file),
                          { color: 0x5de1e6, ligandColor: 0x5de1e6,
                            representation: "ball-and-stick" });
      });
      if (dp.reference_file) {
        work = work.then(function () {
          return scope.load("reference", fileUrl(dp.reference_file),
                            { color: GobViewer.COLOURS.reference,
                              ligandColor: GobViewer.COLOURS.reference,
                              representation: "ball-and-stick" });
        });
      }
      // On the poses, not on the protein: the receptor is context and fills the
      // panel if the camera is left to frame everything.
      return work.then(function () { return scope.focusLigand("top"); });
    }).then(function () {
      var best = dp.best_by_rmsd || {};
      document.getElementById("poses-hud").textContent =
        dp.count + " poses \u00b7 top " +
        (dp.rows[0].rmsd !== null ? dp.rows[0].rmsd + " A" : "unmeasured") +
        (best.rank && best.rank !== 1 ? " \u00b7 nearest is pose " + best.rank +
         " at " + best.rmsd + " A" : "") + " from the crystal";
    }).catch(function (err) {
      document.getElementById("poses-hud").textContent = "viewer unavailable: " + err.message;
    });

    document.querySelectorAll("#poses-toggles button").forEach(function (button) {
      button.addEventListener("click", function () {
        if (!posesScope) return;
        var on = button.getAttribute("aria-pressed") !== "true";
        button.setAttribute("aria-pressed", on ? "true" : "false");
        posesScope.setVisible(button.getAttribute("data-layer"), on);
      });
    });
  }

  function startOverlay() {
    var host = document.getElementById("overlay-viewer");
    if (!host) return;
    GobViewer.create(host).then(function (scope) {
      overlayScope = scope;
      var work = scope.load("md_final", fileUrl(stateFile("md_final")),
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
      var fit = (card.superposed || {}).md_final || {};
      document.getElementById("overlay-hud").textContent =
        (fit.basis || "pocket Ca") + " superposition onto " + (card.superposed_onto || "the model") +
        " · " + (fit.atoms || g.pocket_ca_atoms || 0) + " atoms" +
        (fit.rmsd !== undefined && fit.rmsd !== null ? " · fit " + fit.rmsd + " A" : "") +
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
          overlayScope.load(state, fileUrl(stateFile(state)), { color: GobViewer.COLOURS[state] });
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
  startMotion();
  startPoses();
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
