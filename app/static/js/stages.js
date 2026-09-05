/* The one animation in the app: when a results page loads, the rail's completed
   stages come up green from top to bottom. Everything else changes instantly,
   and prefers-reduced-motion turns even this off (in CSS). */
(function () {
  "use strict";
  document.addEventListener("DOMContentLoaded", function () {
    var rail = document.querySelector(".rail-list.animate");
    if (!rail) return;
    // The animation is declared in CSS with per-cell delays; re-adding the
    // class after paint is what starts it, so a cached page still animates.
    rail.classList.remove("animate");
    requestAnimationFrame(function () { rail.classList.add("animate"); });
  });
})();
