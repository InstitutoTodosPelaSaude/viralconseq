"use strict";

document.title = `${DATA.run.label} — viralconseq run report`;

// The viewer's theme and depth-axis choices, kept in this browser across visits. Storage can be
// missing or refuse access (private windows, blocked site data); the page then uses its defaults.
const PREFERENCE_KEYS = { theme: "viralconseq.report.theme", logDepth: "viralconseq.report.logDepth", modal: "viralconseq.report.modal" };

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

// The opening theme (?theme, else the stored choice, else the system's) and the embed mode are set by
// report/boot.js in <head>, before the first paint.

deferAlpineStartPastFirstPaint();

// Alpine (vendor/alpine.min.js, loaded right after this file) starts itself on a microtask as soon as its script
// has run, before the browser has painted anything, and rendering the report takes it ~70 ms; meanwhile the
// window shows the browser's own blank page, in whatever colour the browser picks. This holds that one start
// until just after the first paint, so the page shows its themed background (its content still x-cloak'ed) at
// once. Only the start is held: any other microtask passes straight through. The timeout covers a page whose
// frames are not running yet (an iframe out of view).
function deferAlpineStartPastFirstPaint() {
  const queueMicrotask = window.queueMicrotask.bind(window);
  window.queueMicrotask = (callback) => {
    if (!/\.start\(\)/.test(String(callback))) return queueMicrotask(callback);
    window.queueMicrotask = queueMicrotask;
    let started = false;
    const start = () => { if (!started) { started = true; callback(); } };
    requestAnimationFrame(() => setTimeout(start, 0));
    setTimeout(start, 100);
  };
}

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
    selectedKey: IS_EMBED ? (EMBED_ROW?.key ?? null) : ((RANKED_ROWS[0] ?? ROWS[0])?.key ?? null),
    dark: isDark(),
    logDepth: readPreference(PREFERENCE_KEYS.logDepth) !== "false",
    // "Open in modal": a sample picked from the table or a chart opens in a modal over the page. Never when
    // embedded, where the detail card is the whole page.
    modalMode: !IS_EMBED && readPreference(PREFERENCE_KEYS.modal) === "true",
    modalOpen: false,

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

    // A pick from the table or a chart (scrollIntoView) shows the detail: in the modal, or by scrolling to it.
    // Picks made inside the detail (Prev/Next, a segment) keep it where it is.
    select(key, { scrollIntoView = true } = {}) {
      this.selectedKey = key;
      if (IS_EMBED || !scrollIntoView) return;
      if (this.modalMode) this.openModal();
      else Alpine.nextTick(scrollDetailIntoView);
    },

    openModal() {
      if (!this.modalOpen) modalReturnFocus = document.activeElement;
      this.modalOpen = true;
      document.documentElement.classList.add("modal-open");
      Alpine.nextTick(() => document.getElementById("detail-close")?.focus());
    },

    closeModal() {
      if (!this.modalOpen) return;
      this.modalOpen = false;
      document.documentElement.classList.remove("modal-open");
      if (modalReturnFocus && modalReturnFocus.isConnected) modalReturnFocus.focus();
      modalReturnFocus = null;
    },

    toggleModalMode() {
      this.modalMode = !this.modalMode;
      writePreference(PREFERENCE_KEYS.modal, String(this.modalMode));
      if (!this.modalMode) this.closeModal();
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

// Embedded in an iframe, the page tells its host how tall it is, so the frame can fit it without scrolling:
// window.addEventListener("message", (e) => { if (e.data?.type === "viralconseq:resize") frame.style.height = e.data.height + "px"; });
function reportHeightToHost() {
  if (!IS_EMBED || window.parent === window) return;
  let last = 0;
  const post = () => {
    const height = Math.ceil(document.documentElement.getBoundingClientRect().height);
    if (height === last) return;
    last = height;
    window.parent.postMessage({ type: "viralconseq:resize", sample: EMBED_SAMPLE, height }, "*");
  };
  new ResizeObserver(post).observe(document.body);
  post();
}

document.addEventListener("alpine:initialized", reportHeightToHost);

// Kept outside the store: Alpine would wrap a stored element in a proxy, which cannot be focused.
let modalReturnFocus = null;

function scrollDetailIntoView() {
  document.getElementById("detail").scrollIntoView({ behavior: "smooth", block: "nearest" });
}

function selectRow(key) {
  Alpine.store("ui").select(key);
}
