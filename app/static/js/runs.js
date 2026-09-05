/* The Runs table. Private runs are fetched with the owner keys this browser
   holds; without a key they are absent from the response entirely, so there is
   nothing here to hide. */
(function () {
  "use strict";
  var rows = [];

  function keys() {
    try { return JSON.parse(localStorage.getItem("gobsmacked_keys") || "{}"); }
    catch (err) { return {}; }
  }

  function anyKey() {
    var all = keys();
    var names = Object.keys(all);
    return names.length ? all[names[names.length - 1]] : "";
  }

  function load() {
    // One key per request: the API matches a single token hash, so several
    // owned runs need several requests. They are merged by job ID here.
    var all = keys();
    var tokens = Object.keys(all).map(function (k) { return all[k]; });
    var unique = tokens.filter(function (t, i) { return tokens.indexOf(t) === i; });
    var requests = [fetch("/api/runs").then(function (r) { return r.json(); })];
    unique.forEach(function (token) {
      requests.push(fetch("/api/runs", { headers: { "X-Owner-Token": token } })
        .then(function (r) { return r.json(); }));
    });
    Promise.all(requests).then(function (results) {
      var seen = {};
      rows = [];
      results.forEach(function (result) {
        (result.runs || []).forEach(function (run) {
          if (seen[run.job_id]) return;
          seen[run.job_id] = true;
          rows.push(run);
        });
      });
      rows.sort(function (a, b) { return a.created < b.created ? 1 : -1; });
      render();
    });
  }

  function render() {
    var q = document.getElementById("search").value.trim().toLowerCase();
    var family = document.getElementById("family-filter").value;
    var status = document.getElementById("status-filter").value;
    var grade = document.getElementById("grade-filter").value;
    var body = document.getElementById("runs-body");
    var shown = rows.filter(function (run) {
      if (family && run.family !== family) return false;
      if (status && run.status !== status) return false;
      if (grade && run.grade !== grade) return false;
      if (!q) return true;
      return [run.job_id, run.uniprot, run.protein_name, run.ligand_name, run.reference_pdb,
              run.title].join(" ").toLowerCase().indexOf(q) >= 0;
    });
    if (!shown.length) {
      body.innerHTML = "<tr><td colspan='10' class='muted'>No runs match.</td></tr>";
      return;
    }
    body.innerHTML = shown.map(function (run) {
      return "<tr>" +
        "<td><a href='" + run.url + "'>" + run.job_id.slice(-12) + "</a>" +
        (run.visibility === "private" ? " <span class='pill private'>private</span>" : "") + "</td>" +
        "<td class='muted'>" + (run.created || "").slice(0, 10) + "</td>" +
        "<td>" + (run.uniprot || "&mdash;") + "<br><span class='muted'>" +
        (run.protein_name || "").slice(0, 40) + "</span></td>" +
        "<td>" + (run.ligand_name || "&mdash;") + "</td>" +
        "<td>" + (run.family || "") + "</td>" +
        "<td>" + (run.reference_pdb || "<span class='muted'>none</span>") + "</td>" +
        "<td><span class='ministrip'>" + (run.stages || []).map(function (s) {
          return "<i class='" + s + "'></i>";
        }).join("") + "</span></td>" +
        "<td class='num'>" + (run.gobsmack_score !== null && run.gobsmack_score !== undefined
          ? run.gobsmack_score.toFixed(1) : "&mdash;") + "</td>" +
        "<td><span class='grade-tile " + (run.grade || "none") + "'>" +
        (run.grade || "&mdash;") + "</span></td>" +
        "<td>" + modeLed(run) + "</td></tr>";
    }).join("");
  }

  function modeLed(run) {
    if (run.mode_match === 1) return "<i class='led'></i>match";
    if (run.mode_match === 0) return "<i class='led amber'></i>differs";
    return "<i class='led off'></i>&mdash;";
  }

  document.addEventListener("DOMContentLoaded", function () {
    ["search", "family-filter", "status-filter", "grade-filter"].forEach(function (id) {
      document.getElementById(id).addEventListener("input", render);
    });
    document.getElementById("owner-key-btn").addEventListener("click", function () {
      var key = document.getElementById("owner-key").value.trim();
      if (!key) return;
      try {
        var all = keys();
        all["pasted_" + key.slice(0, 6)] = key;
        localStorage.setItem("gobsmacked_keys", JSON.stringify(all));
      } catch (err) { /* nothing to remember in private browsing */ }
      load();
    });
    load();
  });
})();
