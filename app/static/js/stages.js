/* The one animation in the app: when a results page loads, the stage strip's
   completed cells come up green from left to right. Everything else changes
   instantly, and prefers-reduced-motion turns even this off (in CSS). */
(function () {
  "use strict";
  document.addEventListener("DOMContentLoaded", function () {
    var strip = document.querySelector(".stages.animate");
    if (!strip) return;
    // The animation is declared in CSS with per-cell delays; adding the class
    // after paint is what starts it, so a cached page still animates.
    requestAnimationFrame(function () { strip.classList.add("animate"); });
  });
})();
