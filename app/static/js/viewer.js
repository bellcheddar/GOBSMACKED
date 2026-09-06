/* Mol* wrapper: one scope, several structures, and a residue picker.

   Mol*'s "viewer" build exports molstar.Viewer and very little else. There is no
   selection query language on the global, so residues are addressed by walking
   the loaded model's own tables, which are reachable through the structure. That
   is the one place this file reaches past the public surface, and it is confined
   to residueLoci() and residueAt().

   Colours follow the app's tokens: model grey, MD-final phosphor, reference
   amber, apo purple. They are passed in by the caller rather than decided here,
   so the legend on the page and the colours in the scene cannot drift. */
(function () {
  "use strict";

  var OPTIONS = {
    /* No extensions. Mol* enables all of them by default, and Volumes &
       Segmentations fetches a listing from a remote server the moment a viewer
       is created: a network call from a page that otherwise talks only to its
       own origin, for a feature nothing here uses. */
    extensions: [],
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowSequence: false,
    layoutShowLog: false,
    layoutShowLeftPanel: false,
    layoutShowRemoteState: false,
    viewportShowExpand: false,
    viewportShowControls: false,
    viewportShowSettings: false,
    viewportShowSelectionMode: false,
    viewportShowAnimation: false,
    viewportShowTrajectoryControls: false,
    pdbProvider: "rcsb",
  };

  var COLOURS = {
    model: 0x9fb0c7,
    md_final: 0x5de1e6,
    reference: 0xffb454,
    apo: 0xc39cff,
    pose1: 0x7ee2a8,
  };

  // The ligand's colour in every pane: warm, and off the scale the protein
  // states are drawn from, so it never reads as a fourth state.
  var LIGAND = 0xff5c8a;

  var NON_LIGAND = { HOH: 1, DOD: 1, SO4: 1, PO4: 1, GOL: 1, EDO: 1, NA: 1, CL: 1, ZN: 1, MG: 1 };

  var STANDARD_RESIDUE = {};
  ("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL " +
   "HID HIE HIP CYX ASH GLH LYN MSE SEC PYL").split(" ").forEach(function (name) {
    STANDARD_RESIDUE[name] = 1;
  });

  function compId(hierarchy, element, residueIndex) {
    if (hierarchy.atoms && hierarchy.atoms.auth_comp_id) {
      return hierarchy.atoms.auth_comp_id.value(element);
    }
    if (hierarchy.residues && hierarchy.residues.auth_comp_id) {
      return hierarchy.residues.auth_comp_id.value(residueIndex);
    }
    return "";
  }

  function hasWebGL() {
    try {
      var canvas = document.createElement("canvas");
      return !!(window.WebGLRenderingContext &&
                (canvas.getContext("webgl") || canvas.getContext("experimental-webgl")));
    } catch (err) {
      return false;
    }
  }

  function allIndices(unit) {
    var indices = new Int32Array(unit.elements.length);
    for (var i = 0; i < indices.length; i++) indices[i] = i;
    return indices;
  }

  /* Every atom of one residue, as a Mol* Loci.

     A Loci element is {unit, indices}, where indices are offsets INTO
     unit.elements rather than element ids. Getting that wrong yields a loci that
     is accepted, highlights the wrong atoms and never errors. */
  function residueLoci(structure, chainId, seqId) {
    var elements = [];
    structure.units.forEach(function (unit) {
      var hierarchy = unit.model.atomicHierarchy;
      if (!hierarchy || !hierarchy.residueAtomSegments) return;
      var indices = [];
      for (var i = 0; i < unit.elements.length; i++) {
        var element = unit.elements[i];
        var chainIndex = hierarchy.chainAtomSegments.index[element];
        if (chainId && hierarchy.chains.auth_asym_id.value(chainIndex) !== chainId) continue;
        var residueIndex = hierarchy.residueAtomSegments.index[element];
        if (hierarchy.residues.auth_seq_id.value(residueIndex) === seqId) indices.push(i);
      }
      if (indices.length) elements.push({ unit: unit, indices: new Int32Array(indices) });
    });
    if (!elements.length) return null;
    return { kind: "element-loci", structure: structure, elements: elements };
  }

  function manyResidueLoci(structure, residues) {
    var elements = [];
    var wanted = {};
    residues.forEach(function (r) { wanted[r.chain + ":" + r.seqId] = true; });
    structure.units.forEach(function (unit) {
      var hierarchy = unit.model.atomicHierarchy;
      if (!hierarchy || !hierarchy.residueAtomSegments) return;
      var indices = [];
      for (var i = 0; i < unit.elements.length; i++) {
        var element = unit.elements[i];
        var chainIndex = hierarchy.chainAtomSegments.index[element];
        var residueIndex = hierarchy.residueAtomSegments.index[element];
        var key = hierarchy.chains.auth_asym_id.value(chainIndex) + ":" +
                  hierarchy.residues.auth_seq_id.value(residueIndex);
        if (wanted[key]) indices.push(i);
      }
      if (indices.length) elements.push({ unit: unit, indices: new Int32Array(indices) });
    });
    if (!elements.length) return null;
    return { kind: "element-loci", structure: structure, elements: elements };
  }

  /* The chain and residue number under a clicked loci, or null. */
  function residueAt(loci) {
    if (!loci || loci.kind !== "element-loci" || !loci.elements.length) return null;
    var entry = loci.elements[0];
    if (!entry.indices || !entry.indices.length) return null;
    var unit = entry.unit;
    var element = unit.elements[entry.indices[0]];
    var hierarchy = unit.model.atomicHierarchy;
    if (!hierarchy || !hierarchy.residueAtomSegments) return null;
    var chainIndex = hierarchy.chainAtomSegments.index[element];
    var residueIndex = hierarchy.residueAtomSegments.index[element];
    return {
      chain: hierarchy.chains.auth_asym_id.value(chainIndex),
      seqId: hierarchy.residues.auth_seq_id.value(residueIndex),
      name: hierarchy.atoms.auth_comp_id
        ? hierarchy.atoms.auth_comp_id.value(element)
        : hierarchy.residues.auth_comp_id.value(residueIndex),
    };
  }

  function Scope(viewer, host) {
    this.viewer = viewer;
    this.plugin = viewer.plugin;
    this.host = host;
    this.entries = {};          // name -> hierarchy structure entry
    this.primary = null;
    this.pocket = [];
  }

  Scope.prototype.data = function (name) {
    var entry = name ? this.entries[name] : this.primary;
    return entry && entry.cell.obj ? entry.cell.obj.data : null;
  };

  /* Load a structure by name. Loading the same name again replaces it, so a
     panel that reloads after a fetch does not accumulate copies. */
  Scope.prototype.load = function (name, url, options) {
    var self = this;
    options = options || {};
    return self.remove(name).then(function () {
      return fetch(url);
    }).then(function (response) {
      if (!response.ok) throw new Error("structure request failed (" + response.status + ")");
      return response.text();
    }).then(function (text) {
      var format = options.format ||
        (url.indexOf(".cif") >= 0 || text.indexOf("data_") === 0 ? "mmcif" : "pdb");
      var before = self.plugin.managers.structure.hierarchy.current.structures.length;
      return self.viewer.loadStructureFromData(text, format).then(function () {
        var all = self.plugin.managers.structure.hierarchy.current.structures;
        if (all.length <= before) throw new Error("the file held no structure");
        var entry = all[all.length - 1];
        self.entries[name] = entry;
        if (!self.primary || options.primary) self.primary = entry;
        var colour = options.color !== undefined ? options.color : COLOURS[name];
        // Representation FIRST, theme second. updateRepresentations installs a
        // fresh representation of the new type, and a fresh representation
        // carries the default colour theme: doing it the other way round
        // applied the colour and then threw it away, which is why ten poses
        // drawn as ball-and-stick came out in element colours with no way to
        // tell the chosen pose from the rest.
        var swap = Promise.resolve();
        if (options.representation) {
          var manager = self.plugin.managers.structure.component;
          swap = Promise.all((entry.components || []).map(function (component) {
            if (!component.representations.length) return null;
            return manager.updateRepresentations(
              [component], component.representations[0],
              { type: { name: options.representation, params: {} } });
          }));
        }
        return swap.then(function () {
          return self.theme(entry, colour, options.ligandColor);
        });
      });
    }).then(function () { return self; });
  };


  /* Protein in the state's colour, ligand in its own.

     Applied per component: colouring every component uniformly makes the
     ligand the same cyan as the 300 residues around it, which is the one thing
     in the scene the page is about. */
  Scope.prototype.theme = function (entry, colour, ligandColour) {
    var manager = this.plugin.managers.structure.component;
    var polymer = [], ligand = [];
    (entry.components || []).forEach(function (component) {
      var key = component.key || "";
      (key.indexOf("ligand") >= 0 ? ligand : polymer).push(component);
    });
    var work = [];
    if (polymer.length && colour !== undefined) {
      work.push(manager.updateRepresentationsTheme(
        polymer, { color: "uniform", colorParams: { value: colour } }));
    }
    if (ligand.length) {
      work.push(manager.updateRepresentationsTheme(
        ligand, { color: "uniform",
                  colorParams: { value: ligandColour === undefined ? LIGAND : ligandColour } }));
    }
    return Promise.all(work);
  };

  /* Sit on the ligand, which is what every one of these panes is about.

     Component identity is read defensively: whether comp ids live on the atom
     table or the residue table depends on how the file was parsed, and reading
     the wrong one throws inside a promise where the only visible symptom is the
     HUD showing "Cannot read properties of undefined". */
  Scope.prototype.focusLigand = function (name) {
    var data = this.data(name);
    if (!data) return false;
    var elements = [];
    data.units.forEach(function (unit) {
      var hierarchy = unit.model.atomicHierarchy;
      if (!hierarchy || !hierarchy.residueAtomSegments) return;
      var indices = [];
      for (var i = 0; i < unit.elements.length; i++) {
        var element = unit.elements[i];
        var residueIndex = hierarchy.residueAtomSegments.index[element];
        var comp = compId(hierarchy, element, residueIndex);
        if (!comp || STANDARD_RESIDUE[comp] || NON_LIGAND[comp]) continue;
        indices.push(i);
      }
      if (indices.length) elements.push({ unit: unit, indices: new Int32Array(indices) });
    });
    if (!elements.length) return false;
    // extraRadius pulls back far enough to show the pocket the ligand is in.
    // Framing the ligand alone fills the pane with a drug and clips away the
    // thing it is bound to, which is the half the reader is judging.
    this.plugin.managers.camera.focusLoci(
      { kind: "element-loci", structure: data, elements: elements },
      { extraRadius: 12, minRadius: 20, durationMs: 0 });
    return true;
  };

  Scope.prototype.remove = function (name) {
    var entry = this.entries[name];
    if (!entry) return Promise.resolve(false);
    delete this.entries[name];
    if (this.primary === entry) this.primary = null;
    try {
      return Promise.resolve(this.plugin.managers.structure.hierarchy.remove([entry]))
        .then(function () { return true; });
    } catch (err) {
      return Promise.resolve(false);
    }
  };

  Scope.prototype.setVisible = function (name, visible) {
    var entry = this.entries[name];
    if (!entry) return false;
    var components = entry.components || [];
    var hidden = [];
    components.forEach(function (component) {
      (component.representations || []).forEach(function (rep) {
        hidden.push(rep.cell.transform.ref);
      });
    });
    var manager = this.plugin.managers.structure.component;
    components.forEach(function (component) {
      manager.toggleVisibility([component], visible);
    });
    return true;
  };

  Scope.prototype.has = function (name) { return !!this.entries[name]; };

  /* PLIP's interactions, drawn as dashed lines.

     Mol*'s viewer build exports Viewer and almost nothing else: there is no
     shape builder to draw a line with. So each interaction arrives as a tiny
     structure of two-atom fragments joined by CONECT records, one file per
     interaction type, and a run of those fragments along the interaction vector
     reads as a dash. Uniform colour per file is what makes the type legible. */
  Scope.prototype.showInteractions = function (lines, urlFor) {
    var self = this;
    var previous = Object.keys(self.entries).filter(function (name) {
      return name.indexOf("interaction:") === 0;
    });
    return Promise.all(previous.map(function (name) { return self.remove(name); }))
      .then(function () {
        return (lines || []).reduce(function (chain, entry) {
          return chain.then(function () {
            return self.load("interaction:" + entry.type, urlFor(entry.file),
                             { color: entry.colour, format: "pdb",
                               representation: "line" });
          });
        }, Promise.resolve());
      });
  };

  Scope.prototype.hideInteractions = function () {
    var self = this;
    return Promise.all(Object.keys(self.entries)
      .filter(function (name) { return name.indexOf("interaction:") === 0; })
      .map(function (name) { return self.remove(name); }));
  };

  Scope.prototype.clearAll = function () {
    this.entries = {};
    this.primary = null;
    this.pocket = [];
    return Promise.resolve(this.plugin.clear());
  };

  /* Pocket residues drawn as sticks over the cartoon, in phosphor. Called on
     every selection change, so it replaces its own component rather than adding
     one each time. */
  Scope.prototype.showPocket = function (residues, colour) {
    var self = this;
    self.pocket = residues || [];
    var data = self.data();
    if (!data) return Promise.resolve(false);
    var loci = self.pocket.length ? manyResidueLoci(data, self.pocket) : null;
    self.plugin.managers.interactivity.lociSelects.deselectAll();
    if (loci) {
      self.plugin.managers.structure.selection.fromLoci("set", loci);
    }
    return Promise.resolve(true);
  };

  Scope.prototype.focus = function (residues) {
    var data = this.data();
    if (!data) return false;
    var loci = residues && residues.length ? manyResidueLoci(data, residues) : null;
    if (!loci) return this.resetCamera();
    this.plugin.managers.camera.focusLoci(loci);
    return true;
  };

  Scope.prototype.resetCamera = function () {
    var data = this.data();
    if (!data) {
      if (this.plugin.managers.camera) this.plugin.managers.camera.reset();
      return false;
    }
    this.plugin.managers.camera.focusLoci({
      kind: "element-loci", structure: data,
      elements: data.units.map(function (unit) {
        return { unit: unit, indices: allIndices(unit) };
      }),
    });
    return true;
  };

  Scope.prototype.highlightResidue = function (chain, seqId) {
    var data = this.data();
    if (!data) return;
    var loci = residueLoci(data, chain, seqId);
    if (loci) this.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci: loci });
  };

  Scope.prototype.onResidueClick = function (callback) {
    var self = this;
    this.plugin.behaviors.interaction.click.subscribe(function (event) {
      var hit = residueAt(event.current && event.current.loci);
      if (hit) callback(hit);
    });
  };

  /* An image of the current scene, for the Report page. */
  Scope.prototype.snapshot = function () {
    try {
      return this.plugin.helpers.viewportScreenshot.getImageDataUri();
    } catch (err) {
      return Promise.resolve(null);
    }
  };

  function create(hostId) {
    var host = typeof hostId === "string" ? document.getElementById(hostId) : hostId;
    if (!host) return Promise.reject(new Error("no viewer element"));
    if (!hasWebGL()) return Promise.reject(new Error("this browser has no WebGL"));
    return molstar.Viewer.create(host, OPTIONS).then(function (viewer) {
      // Mol* paints its own light background over the scope's grid. Transparent
      // lets the panel show through, which is what makes these read as
      // instruments rather than as embedded viewers.
      // Mol*'s default clear colour is a cream (#fcfbf3) that reads as a hole
      // in the panel. transparentBackground is set false deliberately: with it
      // true the canvas has no alpha to be transparent into and clears to
      // white, which looks like the prop was ignored. An opaque clear at the
      // scope's own colour is what actually works.
      try {
        viewer.plugin.canvas3d.setProps({
          transparentBackground: false,
          renderer: { backgroundColor: 0x111826 },
          camera: { helper: { axes: { name: "off", params: {} } } },
        });
      } catch (err) { /* an older build without these props still renders */ }
      var scope = new Scope(viewer, host);
      // Mol* sizes its canvas once. In a grid whose rows stretch, the container
      // changes size after layout and on every window resize, and without this
      // the canvas keeps its first size and spills out of the panel.
      var nudge = function () {
        try { viewer.plugin.handleResize(); } catch (err) { /* nothing to do */ }
      };
      window.addEventListener("resize", nudge);
      if (window.ResizeObserver) {
        var observer = new ResizeObserver(nudge);
        observer.observe(host);
        scope.observer = observer;
      }
      requestAnimationFrame(nudge);
      return scope;
    });
  }

  window.GobViewer = { create: create, COLOURS: COLOURS };
})();
