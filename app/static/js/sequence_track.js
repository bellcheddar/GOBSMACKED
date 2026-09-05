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
  // A band above the ruler for the slanted labels. Everything else is offset
  // below it, so the geometry stays readable as one column of constants.
  // Tall enough for three diagonal lanes. A 16-character label at 9 px is about
  // 80 px long, and at 45 degrees that is 57 px of height, plus two lane steps.
  var LABEL_BAND = 106;
  var LABEL_ANGLE = -45;
  // Two labels rotated 45 degrees are separated by (dx + LANE_STEP) / sqrt(2),
  // not by LANE_STEP: a vertical offset partly slides a label along its own
  // diagonal instead of away from its neighbour. The gatekeeper and the hinge
  // are one residue apart, so dx is 2 px, and a 14 px step left them 11 px
  // apart, which is the height of the text. 22 px gives 17 px of clearance.
  var LANE_STEP = 22;
  var LANES = 3;
  // 9 px Sora averages a little under 5 px a character. Only used to decide
  // which lane a label goes in, so an approximation is fine.
  var CHAR_PX = 4.8;

  var TRACK_Y = LABEL_BAND + 26;
  var TRACK_H = 16;
  var DOMAIN_Y = LABEL_BAND + 46;
  var DOMAIN_H = 13;
  var FEATURE_Y = LABEL_BAND + 63;
  var FEATURE_H = 6;
  var HEIGHT = LABEL_BAND + 78;
  var NS = "http://www.w3.org/2000/svg";

  function el(name, attrs) {
    var node = document.createElementNS(NS, name);
    Object.keys(attrs || {}).forEach(function (k) {
      // setAttribute stringifies, so a null would land as the literal "null"
      // and quietly break whatever it was meant to control.
      if (attrs[k] === null || attrs[k] === undefined) return;
      node.setAttribute(k, attrs[k]);
    });
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
      var full = domain.pfam + " " + domain.name + " (" + domain.start + "-" + domain.end + ")";
      band.appendChild(_title(full));
      svg.appendChild(band);
      // Written into the band itself when it is wide enough to hold something
      // legible, and left to the tooltip when it is not: a clipped half-word is
      // worse than no word.
      if (w > 70) {
        var caption = el("text", {
          x: x + 4, y: DOMAIN_Y + DOMAIN_H - 4, class: "domain-label",
        });
        caption.textContent = w > 300 ? full : domain.pfam;
        svg.appendChild(caption);
      }
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

    // The named positions from the Family panel, written onto the graphic.
    //
    // Slanting alone is not enough. EGFR's hinge is residues 791, 792 and 793,
    // and the gatekeeper is 790: at two pixels a residue those four labels
    // start within six pixels of each other, and a 45 degree rotation still
    // stacks them on top of one another. So runs of the same landmark are
    // merged into one label ("hinge 791-793"), and what is left is dealt into
    // diagonal lanes, each with a longer tick so it is clear which position it
    // belongs to.
    layoutPositions(data.positions, CELL).forEach(function (position) {
      var x = position.x;
      var top = LABEL_BAND - position.lane * LANE_STEP;
      var tick = el("line", { x1: x, y1: top, x2: x, y2: TRACK_Y, class: "tick" });
      tick.appendChild(_title(position.text));
      svg.appendChild(tick);

      var label = el("text", {
        x: x, y: top - 3, class: "position-label", "text-anchor": "end",
        transform: "rotate(" + LABEL_ANGLE + " " + x + " " + (top - 3) + ")",
      });
      // No <title> here: the label is already the words, and a title inside a
      // <text> makes its textContent read back doubled.
      label.textContent = position.text;
      svg.appendChild(label);
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

    // Open on the landmarks, not on residue 1 and not on whichever Pfam domain
    // happens to come first. EGFR's first domain starts at 57 and its kinase
    // landmarks are at 745 to 856: scrolling to the former put every label this
    // track exists to show off the right-hand edge.
    var focusOn = ((data.positions || [])[0] || {}).residue
      || ((data.pfam || [])[0] || {}).start;
    if (focusOn) {
      var target = (focusOn - 1) * CELL - this.host.clientWidth * 0.25;
      this.host.scrollLeft = Math.max(0, target);
    }
    this.watchResize(data);
  };

  /* Merge runs of the same landmark, then deal what is left into lanes.

     Merging first is what makes the lanes enough: three hinge residues become
     one label rather than three competing for the same six pixels. */
  function layoutPositions(positions, cell) {
    var sorted = (positions || []).slice().sort(function (a, b) {
      return a.residue - b.residue;
    });

    var merged = [];
    sorted.forEach(function (position) {
      // "hinge 1", "hinge 2" -> "hinge". A trailing index is a part number,
      // not a different landmark.
      var base = String(position.label).replace(/\s+\d+$/, "");
      var last = merged[merged.length - 1];
      if (last && last.base === base && position.residue - last.last <= 3) {
        last.last = position.residue;
        return;
      }
      merged.push({ base: base, first: position.residue, last: position.residue });
    });

    var laneEnds = [];
    return merged.map(function (item) {
      var text = item.base + " " + item.first +
        (item.last !== item.first ? "-" + item.last : "");
      var x = (item.first - 1) * cell + cell / 2;
      // Rotated 45 degrees, a label of length L occupies L * cos(45) of
      // horizontal room to the left of its anchor.
      var footprint = text.length * CHAR_PX * 0.71;
      var lane = 0;
      while (lane < LANES && laneEnds[lane] !== undefined && x - footprint < laneEnds[lane]) {
        lane++;
      }
      if (lane === LANES) lane = LANES - 1;   // out of lanes: stack rather than drop
      laneEnds[lane] = x;
      return { text: text, x: x, lane: lane };
    });
  }

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
