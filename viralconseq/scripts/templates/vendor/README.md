# Vendored libraries for the run report

The report's components (`../report/*.js`) are ES modules that import these by the
bare names in the template's import map (`preact`, `preact/hooks`, `htm`).
`build_report.py` replaces every import-map path with a `data:` URL holding the file,
so each `report.html` carries them and renders with no network.

| File | Package | Source file in the npm tarball | Licence |
|------|---------|--------------------------------|---------|
| `preact.mjs` | `preact@11.0.0` | `dist/preact.mjs` | MIT (`LICENSE-preact`) |
| `preact-hooks.mjs` | `preact@11.0.0` | `hooks/dist/hooks.mjs` | MIT (`LICENSE-preact`) |
| `htm.mjs` | `htm@3.1.1` | `dist/htm.module.js` | Apache-2.0 (`LICENSE-htm`) |

Each file is the published build unchanged except for a one-line licence banner
prepended and the trailing `//# sourceMappingURL=` comment removed (the map is not
shipped, so a browser would only report it missing).

A vendored module must be a single file whose imports are all bare names listed in the
import map: a module loaded from a `data:` URL cannot resolve a relative import.

To update: `npm pack preact@11 htm@3` into an empty directory, extract each tarball,
copy the files above, re-apply the banner, strip the sourceMappingURL line, update the
versions here, and check the page in a browser (`file://` and the HTTP dev loop).
