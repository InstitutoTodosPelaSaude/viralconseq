"use strict";

document.title = `${DATA.run.label} — viralconseq run report`;

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

    toggleTheme() {
      document.documentElement.dataset.theme = isDark() ? "light" : "dark";
      this.coloursChanged();
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
