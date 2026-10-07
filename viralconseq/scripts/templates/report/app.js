"use strict";
/* The state the page shares, as an Alpine store: report.html reads it as $store.ui.
   Registered on alpine:init, which Alpine fires as it starts, after these scripts. */

document.addEventListener("alpine:init", () => {
  Alpine.store("ui", {
    /* the table's view */
    filter: "all",
    query: "",
    sortKey: "cov_min",
    sortDir: -1,
    /* the row in the detail panel */
    selected: (ranked[0] ?? S[0])?.key ?? null,
    /* bumped when the theme's colours change, so text that names the theme updates */
    colours: 0,

    get rows() { return visibleRows(this); },
    get current() { return S.find((d) => d.key === this.selected) || null; },
    /* the selected row's place in the table, which Previous/Next walk; -1 if filtered out */
    get index() { return this.rows.findIndex((r) => r.key === this.selected); },

    /* A new column sorts descending if numeric, ascending otherwise; the same one flips. */
    sortBy(k) {
      if (k === this.sortKey) this.sortDir = -this.sortDir;
      else { this.sortKey = k; this.sortDir = COLS.find((c) => c.k === k).num ? -1 : 1; }
    },
    /* Selecting from the table or a chart scrolls the detail into view; Previous/Next do not. */
    select(key, scroll = true) {
      this.selected = key;
      if (scroll) Alpine.nextTick(() =>
        document.getElementById("detail").scrollIntoView({ behavior: "smooth", block: "nearest" }));
    },
    step(by) { this.select(this.rows[this.index + by].key, false); },

    /* OS preference by default, explicit once clicked. Charts read colours as they draw. */
    toggleTheme() {
      document.documentElement.dataset.theme = isDark() ? "light" : "dark";
      this.coloursChanged();
    },
    coloursChanged() { this.colours++; redrawAll(); },
  });
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => Alpine.store("ui").coloursChanged());
});

/* Chart marks select rows through this. */
function select(key) { Alpine.store("ui").select(key); }

document.title = `${DATA.run.label} — viralconseq run report`;
