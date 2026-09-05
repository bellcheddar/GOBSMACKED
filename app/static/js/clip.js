/* Attaches the motion clip's source, last and only when it is wanted.
 *
 * The <video> ships with a poster and no src, so the browser paints the first
 * frame immediately and fetches nothing. The real file is attached only after
 * window.load, which keeps half a megabyte of video out of the way while the
 * page is still pulling Mol* and Plotly, both of which are larger and both of
 * which the reader is more likely to interact with.
 *
 * Two ways this deliberately does nothing. A visitor who has asked for reduced
 * motion never downloads the clip, and with JavaScript off the poster simply
 * stays. Both leave the same still image the clip opens on, and the panel's
 * caption says what the animation would have shown, so nothing is lost that
 * appears nowhere else on the page.
 */
(function () {
  "use strict";

  var still = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  function load() {
    var videos = document.querySelectorAll("video.clip[data-src]");
    for (var i = 0; i < videos.length; i++) {
      var video = videos[i];
      video.src = video.dataset.src;
      video.removeAttribute("data-src");
      // Safari rejects the promise rather than throwing when it refuses to
      // autoplay, which surfaces as an unhandled rejection in the console on a
      // page where nothing has actually gone wrong.
      var started = video.play();
      if (started && started.catch) started.catch(function () {});
    }
  }

  if (still) return;
  if (document.readyState === "complete") load();
  else window.addEventListener("load", load);
})();
