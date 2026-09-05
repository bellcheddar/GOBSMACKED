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
        if (colour !== undefined) {
          return self.plugin.managers.structure.component.updateRepresentationsTheme(
            entry.components, { color: "uniform", colorParams: { value: colour } });
        }
      });
    }).then(function () { return self; });
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
      return new Scope(viewer, host);
    });
  }

  window.GobViewer = { create: create, COLOURS: COLOURS };
})();
