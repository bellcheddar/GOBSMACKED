# Vendored browser libraries

Both are served from this directory rather than a CDN so the app has no runtime
dependency on a third party, and so a deploy pins the exact file that was
tested.

| File | Version | Licence | Source |
|---|---|---|---|
| `molstar.js`, `molstar.css` | see `molstar.version` | MIT | https://github.com/molstar/molstar |
| `plotly.min.js` | 2.35.2 | MIT | https://github.com/plotly/plotly.js |

Both are cached hard by nginx (`immutable`), so every reference from a template
carries a `?v=<mtime>` stamp. Replacing a file therefore changes its URL.
