# Vendored libraries for the run report

`build_report.py` inlines every `<script src>` of the template into the page, so each
`report.html` carries these and renders with no network.

| File | Package | Source file in the npm tarball | Licence |
|------|---------|--------------------------------|---------|
| `alpine.min.js` | `alpinejs@3.17.4` | `dist/cdn.min.js` | MIT (`LICENSE-alpine`, from the repository's `v3.17.4` tag: the npm tarball ships none) |

The file is the published build unchanged except for a one-line licence banner
prepended and a trailing `//# sourceMappingURL=` comment removed, if present (the map
is not shipped, so a browser would only report it missing).

`cdn.min.js` is a classic script: it defines `window.Alpine` and starts it on the
next microtask, so the template loads it after the report's own scripts, which
register their stores and components on `alpine:init`.

To update: `npm pack alpinejs@3` into an empty directory, extract it, copy the file
above, re-apply the banner, take `LICENSE.md` from the matching repository tag,
update the version here, and check the page in a browser (`file://` and the HTTP
dev loop).
