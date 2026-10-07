# Vendored libraries for the run report

`build_report.py` inlines every `<script src>` of the template into the page, so each
`report.html` carries these and renders with no network.

| File | Package | Source file in the npm tarball | Licence |
|------|---------|--------------------------------|---------|
| `alpine.min.js` | `alpinejs@3.17.4` | `dist/cdn.min.js` | MIT (`LICENSE-alpine`, from the repository's `v3.17.4` tag: the npm tarball ships none) |
| `chart.umd.min.js` | `chart.js@4.5.1` | `dist/chart.umd.min.js` | MIT (`LICENSE-chartjs`); bundles `@kurkle/color@0.3.2`, MIT (`LICENSE-kurkle-color`) |
| `chartjs-plugin-annotation.min.js` | `chartjs-plugin-annotation@3.1.0` | `dist/chartjs-plugin-annotation.min.js` | MIT (`LICENSE-chartjs-plugin-annotation`) |
| `chartjs-plugin-zoom.min.js` | `chartjs-plugin-zoom@2.2.0` | `dist/chartjs-plugin-zoom.min.js` | MIT (`LICENSE-chartjs-plugin-zoom`) |
| `hammer.min.js` | `hammerjs@2.0.8` | `hammer.min.js` | MIT (`LICENSE-hammerjs`) |

Each file is the published build unchanged except for a trailing `//# sourceMappingURL=`
comment removed, if present (the map is not shipped, so a browser would only report it
missing); the builds keep their own licence banners. Alpine's has one prepended.

All are classic scripts. Load order in the template matters: Chart.js defines
`window.Chart`; hammer.js defines `window.Hammer`, which the zoom plugin reads for
drag-to-pan; both plugins register themselves with `window.Chart` as they load; the
report's own scripts follow; Alpine comes last, because it starts on the next microtask
and its `alpine:init` handlers live in `report/app.js`.

To update: `npm pack <package>@<version>` into an empty directory, extract it, copy
the file above, strip the sourceMappingURL line, take the licence from the package (or,
for Alpine, from the matching repository tag), update the versions here, and check the
page in a browser (`file://` and the HTTP dev loop).
