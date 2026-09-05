/* The sequence track: one SVG, no library.

   Residues are 6 px cells on a scrolling strip, Pfam domains are bands beneath
   them, UniProt features sit under those, and family-specific positions (KLIFS
   pocket residues, GPCRdb microswitches) are phosphor ticks above. Clicking a
   cell toggles that residue in the pocket selection.

   A 1200-residue protein is 7200 px wide, which is why the strip scrolls: the
   alternative is 1 px per residue, and a 1 px cell cannot be clicked. The view
   scrolls itself to the first Pfam domain on render, because that is where the
   binding site is on every target this app is pointed at. */
(function () {
  "use strict";

  // The cell width is computed to fit the panel rather than fixed at 6 px.
  // EGFR at 6 px is a 7,260 px SVG inside a 1,300 px panel, which works only as
  // long as every browser gets overflow and grid minimums exactly right, and is
  // unreadable even when they do. Below MIN_CELL the track scrolls, but it is
  // then scrolling a strip a few times the panel width, not fifty.
  var MAX_CELL = 6;
  var MIN_CELL = 2;
  var TRACK_Y = 26;
  var TRACK_H = 16;
  var DOMAIN_Y = 46;
  var DOMAIN_H = 8;
  var FEATURE_Y = 58;
  var FEATURE_H = 6;
  var HEIGHT = 78;
  var NS = "http://www.w3.org/2000/svg";

  function el(name, attrs) {
    var node = document.createElementNS(NS, name);
    Object.keys(attrs || {}).forEach(function (k) { node.setAttribute(k, attrs[k]); });
    return node;
  }

  function Track(host, options) {
    this.host = typeof host === "string" ? document.getElementById(host) : host;
    this.options = options || {};
    this.selected = {};        // residue number -> true
    this.length = 0;
    this.cells = {};
    this.onToggle = this.options.onToggle || function () {};
    this.onHover = this.options.onHover || function () {};
  }

  Track.prototype.render = function (data) {
    var self = this;
    var sequence = data.sequence || "";
    this.length = sequence.length;
    this.cells = {};
    this.host.innerHTML = "";
    if (!this.length) return;

    var available = this.host.clientWidth || 900;
    var CELL = Math.max(MIN_CELL, Math.min(MAX_CELL, Math.floor(available / this.length)));
    this.cell = CELL;
    var width = this.length * CELL;
    var svg = el("svg", { width: width, height: HEIGHT, viewBox: "0 0 " + width + " " + HEIGHT });

    // Ruler every 50 residues, labelled every 100.
    for (var pos = 50; pos <= this.length; pos += 50) {
      var x = (pos - 1) * CELL;
      svg.appendChild(el("line", { x1: x, y1: TRACK_Y - 4, x2: x, y2: TRACK_Y,
                                   stroke: "#3a4a63", "stroke-width": 1 }));
      if (pos % 100 === 0) {
        var label = el("text", { x: x + 2, y: TRACK_Y - 6 });
        label.textContent = pos;
        svg.appendChild(label);
      }
    }

    // Pfam domains, drawn under the residues.
    (data.pfam || []).forEach(function (domain) {
      var x = (domain.start - 1) * CELL;
      var w = (domain.end - domain.start + 1) * CELL;
      var band = el("rect", { x: x, y: DOMAIN_Y, width: w, height: DOMAIN_H, class: "domain" });
      band.appendChild(_title(domain.pfam + " " + domain.name + " (" + domain.start + "-" + domain.end + ")"));
      svg.appendChild(band);
    });

    // UniProt features, one thin band each.
    (data.features || []).forEach(function (feature) {
      var x = (feature.start - 1) * CELL;
      var w = Math.max(CELL, (feature.end - feature.start + 1) * CELL);
      var band = el("rect", { x: x, y: FEATURE_Y, width: w, height: FEATURE_H, class: "feature" });
      band.appendChild(_title(feature.type + (feature.description ? ": " + feature.description : "")));
      svg.appendChild(band);
    });

    // The residues themselves.
    for (var i = 0; i < this.length; i++) {
      var number = i + 1;
      var cell = el("rect", {
        x: i * CELL, y: TRACK_Y, width: CELL - 1, height: TRACK_H,
        class: "cell", "data-residue": number,
      });
      cell.appendChild(_title(sequence[i] + number));
      svg.appendChild(cell);
      this.cells[number] = cell;
    }

    // Family positions on top, as ticks.
    (data.positions || []).forEach(function (position) {
      var x = (position.residue - 1) * CELL + CELL / 2;
      var tick = el("line", { x1: x, y1: TRACK_Y - 8, x2: x, y2: TRACK_Y, class: "tick" });
      tick.appendChild(_title(position.label + " " + position.residue));
      svg.appendChild(tick);
    });

    svg.addEventListener("click", function (event) {
      var residue = event.target.getAttribute && event.target.getAttribute("data-residue");
      if (residue) self.toggle(parseInt(residue, 10));
    });
    svg.addEventListener("mouseover", function (event) {
      var residue = event.target.getAttribute && event.target.getAttribute("data-residue");
      if (residue) self.onHover(parseInt(residue, 10));
    });

    this.host.appendChild(svg);
    this.repaint();

    var first = (data.pfam || [])[0];
    if (first) this.host.scrollLeft = Math.max(0, (first.start - 20) * CELL);
    this.watchResize(data);
  };

  function _title(text) {
    var node = document.createElementNS(NS, "title");
    node.textContent = text;
    return node;
  }

  /* Re-render on resize: the cell width is a function of the panel width, and a
     window that changes size otherwise leaves the track at its old scale. */
  Track.prototype.watchResize = function (data) {
    var self = this;
    if (self._resizeBound) return;
    self._resizeBound = true;
    var timer = null;
    window.addEventListener("resize", function () {
      window.clearTimeout(timer);
      timer = window.setTimeout(function () {
        var selected = self.residues();
        self.render(data);
        self.setSelected(selected);
      }, 200);
    });
  };

  Track.prototype.toggle = function (residue) {
    if (this.selected[residue]) delete this.selected[residue];
    else this.selected[residue] = true;
    this.repaint();
    this.onToggle(this.residues());
  };

  Track.prototype.setSelected = function (residues) {
    this.selected = {};
    (residues || []).forEach(function (r) { this.selected[r] = true; }, this);
    this.repaint();
  };

  Track.prototype.residues = function () {
    return Object.keys(this.selected).map(Number).sort(function (a, b) { return a - b; });
  };

  Track.prototype.repaint = function () {
    var self = this;
    Object.keys(this.cells).forEach(function (number) {
      var cell = self.cells[number];
      if (self.selected[number]) cell.classList.add("selected");
      else cell.classList.remove("selected");
    });
  };

  window.SequenceTrack = Track;
})();
