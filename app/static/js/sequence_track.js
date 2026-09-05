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

  var TRACK_H = 16;
  var DOMAIN_H = 13;
  var FEATURE_H = 6;
  var BELOW = 78;          // ruler, residues, domains and features, under the band
  var BAND_MIN = 14;       // when there is nothing to label
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

  /* Choose what to show. A 1,210-residue precursor at the two-pixel floor is a
     2,420 px strip in a 900 px panel, so whatever is not scrolled to is simply
     not there: the labels this track exists to show were off the right-hand
     edge every time. Showing the domain that holds the landmarks instead means
     they are on screen without anybody having to hunt for them, and the whole
     chain is one chip away. */
  Track.prototype.defaultRange = function (data) {
    var domains = data.pfam || [];
    var first = (data.positions || [])[0];
    if (first) {
      var holding = domains.filter(function (d) {
        return d.start <= first.residue && first.residue <= d.end;
      })[0];
      if (holding) return { from: holding.start, to: holding.end, label: holding.pfam };
    }
    var biggest = domains.slice().sort(function (a, b) {
      return (b.end - b.start) - (a.end - a.start);
    })[0];
    if (biggest) return { from: biggest.start, to: biggest.end, label: biggest.pfam };
    return { from: 1, to: this.length || (data.sequence || "").length, label: "whole chain" };
  };

  Track.prototype.render = function (data, range) {
    var self = this;
    var sequence = data.sequence || "";
    this.length = sequence.length;
    this.data = data;
    this.cells = {};
    this.host.innerHTML = "";
    if (!this.length) return;

    range = range || this.range || this.defaultRange(data);
    range = { from: Math.max(1, range.from), to: Math.min(this.length, range.to),
              label: range.label };
    this.range = range;
    var span = range.to - range.from + 1;

    this.host.appendChild(this.viewChips(data));

    var scroller = el2("div", "track-scroll");
    var available = this.host.clientWidth || 900;
    var CELL = Math.max(MIN_CELL, Math.min(MAX_CELL, Math.floor(available / span)));
    this.cell = CELL;
    var width = span * CELL;

    // The label band is sized to the labels rather than fixed, so a target with
    // two landmarks does not carry the empty headroom a target with eight
    // needs, and one with none carries no band at all.
    var placed = layoutPositions((data.positions || []).filter(function (p) {
      return p.residue >= range.from && p.residue <= range.to;
    }), CELL, range.from);
    var LABEL_BAND = BAND_MIN;
    placed.forEach(function (position) {
      var reach = position.lane * LANE_STEP + position.text.length * CHAR_PX * 0.71 + 8;
      if (reach > LABEL_BAND) LABEL_BAND = reach;
    });
    LABEL_BAND = Math.round(Math.min(LABEL_BAND, 150));

    var TRACK_Y = LABEL_BAND + 26;
    var DOMAIN_Y = LABEL_BAND + 46;
    var FEATURE_Y = LABEL_BAND + 63;
    var HEIGHT = LABEL_BAND + BELOW;

    var svg = el("svg", { width: width, height: HEIGHT, viewBox: "0 0 " + width + " " + HEIGHT });
    var xOf = function (residue) { return (residue - range.from) * CELL; };

    // Ruler: every 50 residues, numbered every 100, or every 10/50 on a short
    // range where every-100 would draw one tick and no numbers at all.
    var step = span > 400 ? 50 : span > 120 ? 25 : 10;
    var numberEvery = step * 2;
    for (var pos = Math.ceil(range.from / step) * step; pos <= range.to; pos += step) {
      var x = xOf(pos);
      svg.appendChild(el("line", { x1: x, y1: TRACK_Y - 4, x2: x, y2: TRACK_Y,
                                   stroke: "#3a4a63", "stroke-width": 1 }));
      if (pos % numberEvery === 0) {
        var label = el("text", { x: x + 2, y: TRACK_Y - 6 });
        label.textContent = pos;
        svg.appendChild(label);
      }
    }

    // Pfam domains, clipped to the range.
    (data.pfam || []).forEach(function (domain) {
      if (domain.end < range.from || domain.start > range.to) return;
      var from = Math.max(domain.start, range.from);
      var to = Math.min(domain.end, range.to);
      var x = xOf(from);
      var w = (to - from + 1) * CELL;
      var band = el("rect", { x: x, y: DOMAIN_Y, width: w, height: DOMAIN_H, class: "domain" });
      var full = domain.pfam + " " + domain.name + " (" + domain.start + "-" + domain.end + ")";
      band.appendChild(_title(full));
      svg.appendChild(band);
      // Written into the band when it is wide enough to hold something
      // legible, and left to the tooltip when it is not: a clipped half-word
      // is worse than no word.
      if (w > 70) {
        var caption = el("text", { x: x + 4, y: DOMAIN_Y + DOMAIN_H - 4, class: "domain-label" });
        caption.textContent = w > 300 ? full : domain.pfam;
        svg.appendChild(caption);
      }
    });

    // UniProt features, one thin band each.
    (data.features || []).forEach(function (feature) {
      if (feature.end < range.from || feature.start > range.to) return;
      var from = Math.max(feature.start, range.from);
      var to = Math.min(feature.end, range.to);
      var band = el("rect", { x: xOf(from), y: FEATURE_Y,
                              width: Math.max(CELL, (to - from + 1) * CELL),
                              height: FEATURE_H, class: "feature" });
      band.appendChild(_title(feature.type + (feature.description ? ": " + feature.description : "")));
      svg.appendChild(band);
    });

    // The residues themselves, keyed by their real number.
    for (var residue = range.from; residue <= range.to; residue++) {
      var cell = el("rect", {
        x: xOf(residue), y: TRACK_Y, width: Math.max(1, CELL - 1), height: TRACK_H,
        class: "cell", "data-residue": residue,
      });
      cell.appendChild(_title(sequence[residue - 1] + residue));
      svg.appendChild(cell);
      this.cells[residue] = cell;
    }

    // The named positions from the Family panel, written onto the graphic.
    //
    // Slanting alone is not enough. EGFR's hinge is residues 791, 792 and 793,
    // and the gatekeeper is 790: at three pixels a residue those four labels
    // start within nine pixels of each other, and a 45 degree rotation still
    // stacks them. So runs of the same landmark are merged into one label
    // ("hinge 791-793"), and what is left is dealt into diagonal lanes, each
    // with a longer tick so it is clear which position it belongs to.
    placed.forEach(function (position) {
      var x = position.x;
      var top = LABEL_BAND - position.lane * LANE_STEP;
      var tick = el("line", { x1: x, y1: top, x2: x, y2: TRACK_Y, class: "tick" });
      tick.appendChild(_title(position.text));
      svg.appendChild(tick);

      var text = el("text", {
        x: x, y: top - 3, class: "position-label", "text-anchor": "end",
        transform: "rotate(" + LABEL_ANGLE + " " + x + " " + (top - 3) + ")",
      });
      // No <title> here: the label is already the words, and a title inside a
      // <text> makes its textContent read back doubled.
      text.textContent = position.text;
      svg.appendChild(text);
    });

    svg.addEventListener("click", function (event) {
      var residue = event.target.getAttribute && event.target.getAttribute("data-residue");
      if (residue) self.toggle(parseInt(residue, 10));
    });
    svg.addEventListener("mouseover", function (event) {
      var residue = event.target.getAttribute && event.target.getAttribute("data-residue");
      if (residue) self.onHover(parseInt(residue, 10));
    });

    scroller.appendChild(svg);
    this.host.appendChild(scroller);
    this.scroller = scroller;
    this.repaint();
    this.watchResize(data);
  };

  /* One chip per Pfam domain, plus the whole chain. */
  Track.prototype.viewChips = function (data) {
    var self = this;
    var bar = el2("div", "track-views");
    var options = (data.pfam || []).map(function (d) {
      return { from: d.start, to: d.end, label: d.pfam, title: d.name };
    });
    options.push({ from: 1, to: this.length, label: "whole chain",
                   title: "all " + this.length + " residues" });
    options.forEach(function (option) {
      var chip = document.createElement("button");
      chip.className = "chip";
      chip.type = "button";
      chip.textContent = option.label + " " + option.from + "-" + option.to;
      chip.title = option.title || "";
      var active = self.range && self.range.from === option.from && self.range.to === option.to;
      chip.setAttribute("aria-pressed", active ? "true" : "false");
      chip.addEventListener("click", function () {
        var selected = self.residues();
        self.render(data, option);
        self.setSelected(selected);
      });
      bar.appendChild(chip);
    });
    return bar;
  };

  /* Merge runs of the same landmark, then deal what is left into lanes.

     Merging first is what makes the lanes enough: three hinge residues become
     one label rather than three competing for the same few pixels. */
  function layoutPositions(positions, cell, from) {
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
      var x = (item.first - from) * cell + cell / 2;
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

  function el2(name, className) {
    var node = document.createElement(name);
    node.className = className;
    return node;
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
        self.render(data, self.range);
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
