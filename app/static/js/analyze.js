/* Upload a results archive and go where the answer is. */
(function () {
  "use strict";
  document.addEventListener("DOMContentLoaded", function () {
    var button = document.getElementById("upload-btn");
    var out = document.getElementById("upload-result");
    button.addEventListener("click", function () {
      var file = document.getElementById("archive").files[0];
      if (!file) { out.innerHTML = "<p class='note' style='color:var(--red)'>Choose the archive first.</p>"; return; }
      var form = new FormData();
      form.append("file", file);
      button.disabled = true;
      button.textContent = "Analysing...";
      out.innerHTML = "<p class='note'>Validating, superposing, running PLIP, classifying. " +
        "This takes seconds.</p>";
      fetch("/api/upload", { method: "POST", body: form })
        .then(function (r) { return r.json().then(function (d) {
          if (!r.ok) throw new Error(d.error || "upload failed"); return d; }); })
        .then(function (data) { window.location = data.url; })
        .catch(function (err) {
          out.innerHTML = "<p class='note' style='color:var(--red)'>" + err.message + "</p>";
          button.disabled = false;
          button.textContent = "Analyse";
        });
    });
  });
})();
