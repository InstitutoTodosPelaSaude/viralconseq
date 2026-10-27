"use strict";

document.title = `${DATA.run.label} — viralconseq run report`;

// The viewer's theme and depth-axis choices, kept in this browser across visits. Storage can be
// missing or refuse access (private windows, blocked site data); the page then uses its defaults.
const PREFERENCE_KEYS = { theme: "viralconseq.report.theme", logDepth: "viralconseq.report.logDepth" };

function readPreference(key) {
  try {
    return localStorage.getItem(key);
  } catch {
    return null;
  }
}

function writePreference(key, value) {
  try {
    localStorage.setItem(key, value);
  } catch {
    // not persisted; the choice still applies to this page
  }
}

// Applied before Alpine starts, so the first paint of every chart already uses the stored theme.
const STORED_THEME = readPreference(PREFERENCE_KEYS.theme);
if (STORED_THEME === "light" || STORED_THEME === "dark") document.documentElement.dataset.theme = STORED_THEME;

document.addEventListener("alpine:init", () => {
  Alpine.store("ui", createUiStore());
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => Alpine.store("ui").coloursChanged());
});

function createUiStore() {
  return {
    filter: "all",
    query: "",
    sortKey: "cov_min",
    sortDirection: -1,
    selectedKey: (RANKED_ROWS[0] ?? ROWS[0])?.key ?? null,
    dark: isDark(),
    logDepth: readPreference(PREFERENCE_KEYS.logDepth) !== "false",

    get visibleRows() {
      return rowsInView(this);
    },

    get selectedRow() {
      return ROWS.find((row) => row.key === this.selectedKey) || null;
    },

    get selectedIndex() {
      return this.visibleRows.findIndex((row) => row.key === this.selectedKey);
    },

    get hasPrevious() {
      return this.selectedIndex > 0;
    },

    get hasNext() {
      return this.selectedIndex >= 0 && this.selectedIndex < this.visibleRows.length - 1;
    },

    sortBy(key) {
      if (key === this.sortKey) {
        this.sortDirection = -this.sortDirection;
        return;
      }
      const column = TABLE_COLUMNS.find((candidate) => candidate.key === key);
      this.sortKey = key;
      this.sortDirection = column.numeric ? -1 : 1;
    },

    select(key, { scrollIntoView = true } = {}) {
      this.selectedKey = key;
      if (scrollIntoView) Alpine.nextTick(scrollDetailIntoView);
    },

    selectNeighbour(offset) {
      this.select(this.visibleRows[this.selectedIndex + offset].key, { scrollIntoView: false });
    },

    toggleLogDepth() {
      this.logDepth = !this.logDepth;
      writePreference(PREFERENCE_KEYS.logDepth, String(this.logDepth));
      rebuildCharts();
    },

    setTheme(theme) {
      if ((theme === "dark") === isDark()) return;
      document.documentElement.dataset.theme = theme;
      writePreference(PREFERENCE_KEYS.theme, theme);
      this.coloursChanged();
    },

    toggleTheme() {
      this.setTheme(isDark() ? "light" : "dark");
    },

    coloursChanged() {
      this.dark = isDark();
      rebuildCharts();
    },
  };
}

function scrollDetailIntoView() {
  document.getElementById("detail").scrollIntoView({ behavior: "smooth", block: "nearest" });
}

function selectRow(key) {
  Alpine.store("ui").select(key);
}
