/* Plotly panels, in the instrument-panel palette.

   Transparent backgrounds, Sora at 12 px, phosphor for the predicted trace and
   amber for anything measured on the crystal, thresholds as red dashed lines
   labelled at the left edge. No mode bar logo, no lasso, no select: this is a
   readout, not a data explorer. */
(function () {
  "use strict";

  var C = {
    phos: "#5de1e6", amber: "#ffb454", red: "#ff5c5c", green: "#7ee2a8",
    grey: "#9fb0c7", line: "#3a4a63", text: "#e8edf5",
  };

  var CONFIG = {
    displaylogo: false, responsive: true,
    modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d", "toggleSpikelines"],
  };

  function layout(title, xTitle, yTitle, extra) {
    var base = {
      title: { text: title, font: { size: 13, color: C.grey }, x: 0, xanchor: "left" },
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      font: { family: "Sora, sans-serif", size: 12, color: C.text },
      xaxis: { title: { text: xTitle, font: { size: 11 } }, gridcolor: C.line,
               zerolinecolor: C.line, linecolor: C.line },
      yaxis: { title: { text: yTitle, font: { size: 11 } }, gridcolor: C.line,
               zerolinecolor: C.line, linecolor: C.line },
      margin: { l: 56, r: 16, t: 34, b: 44 },
      showlegend: false,
      hovermode: "x unified",
    };
    return Object.assign(base, extra || {});
  }

  function threshold(value, label, colour) {
    return {
      type: "line", xref: "paper", x0: 0, x1: 1, y0: value, y1: value,
      line: { color: colour || C.red, width: 1, dash: "dash" },
    };
  }

  function annotation(value, label, colour) {
    return {
      xref: "paper", x: 0, y: value, xanchor: "left", yanchor: "bottom",
      text: label, showarrow: false, font: { size: 11, color: colour || C.red },
    };
  }

  function ligandRmsd(host, dynamics) {
    var traces = [{
      x: dynamics.times_ps, y: dynamics.ligand_rmsd_pose1, type: "scatter", mode: "lines",
      name: "to pose 1", line: { color: C.phos, width: 2 },
    }];
    var reference = dynamics.ligand_rmsd_reference;
    if (reference && reference.series && reference.series.length) {
      traces.push({
        x: dynamics.times_ps.slice(0, reference.series.length), y: reference.series,
        type: "scatter", mode: "lines", name: "to the crystal ligand",
        line: { color: C.amber, width: 2 },
      });
    }
    var l = layout("Ligand RMSD", "ps", "A", {
      showlegend: true,
      legend: { orientation: "h", y: 1.12, font: { size: 11 } },
      shapes: [threshold(2.0)], annotations: [annotation(2.0, "2 A")],
    });
    Plotly.newPlot(host, traces, l, CONFIG);
  }

  function proteinRmsd(host, dynamics) {
    var traces = [
      { x: dynamics.times_ps, y: dynamics.protein_ca_rmsd, type: "scatter", mode: "lines",
        name: "protein Ca", line: { color: C.grey, width: 2 } },
      { x: dynamics.times_ps, y: dynamics.pocket_ca_rmsd, type: "scatter", mode: "lines",
        name: "pocket Ca", line: { color: C.phos, width: 2 } },
    ];
    Plotly.newPlot(host, traces, layout("Backbone RMSD from the first frame", "ps", "A", {
      showlegend: true, legend: { orientation: "h", y: 1.12, font: { size: 11 } },
    }), CONFIG);
  }

  function rmsf(host, dynamics) {
    var data = dynamics.rmsf || {};
    if (!data.residues || !data.values) return;
    var pocket = {};
    (dynamics.pocket_residues || []).forEach(function (r) { pocket[r] = true; });
    var colours = data.residues.map(function (r) { return pocket[r] ? C.phos : C.grey; });
    Plotly.newPlot(host, [{
      x: data.residues, y: data.values, type: "bar",
      marker: { color: colours }, name: "RMSF",
    }], layout("Per-residue fluctuation, pocket residues in phosphor", "residue", "A"), CONFIG);
  }

  function volume(host, dynamics) {
    if (!dynamics.pocket_volume || !dynamics.pocket_volume.length) return;
    var shapes = [];
    var annotations = [];
    if (dynamics.reference_volume) {
      shapes.push(threshold(dynamics.reference_volume, "crystal", C.amber));
      annotations.push(annotation(dynamics.reference_volume, "crystal pocket", C.amber));
    }
    Plotly.newPlot(host, [{
      x: dynamics.times_ps, y: dynamics.pocket_volume, type: "scatter", mode: "lines",
      line: { color: C.phos, width: 2 },
    }], layout("Pocket volume", "ps", "A cubed", { shapes: shapes, annotations: annotations }), CONFIG);
  }

  function contacts(host, dynamics) {
    var data = dynamics.contacts || {};
    if (!data.residues || !data.matrix || !data.matrix.length) return;
    Plotly.newPlot(host, [{
      z: data.matrix, x: dynamics.times_ps, y: data.residues, type: "heatmap",
      colorscale: [[0, "#243044"], [1, C.phos]], showscale: false,
    }], layout("Contact persistence, 4 A heavy atom", "ps", "residue", {
      margin: { l: 62, r: 16, t: 34, b: 44 },
    }), CONFIG);
  }

  window.GobPlots = {
    all: function (dynamics) {
      if (!dynamics || !dynamics.frames) return;
      if (document.getElementById("plot-ligand")) ligandRmsd("plot-ligand", dynamics);
      if (document.getElementById("plot-protein")) proteinRmsd("plot-protein", dynamics);
      if (document.getElementById("plot-rmsf")) rmsf("plot-rmsf", dynamics);
      if (document.getElementById("plot-volume")) volume("plot-volume", dynamics);
      if (document.getElementById("plot-contacts")) contacts("plot-contacts", dynamics);
    },
  };
})();
