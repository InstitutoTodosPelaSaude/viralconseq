/* The run document, the figures derived from it, and the formatting every section shares. */
import { h } from "preact";
import htm from "htm";

/* html`...` builds Preact elements from a tagged template: JSX without a build step.
   htm drops whitespace that contains a line break at the edge of a text, so keep a
   space that matters on the same line as the word, tag or ${} it separates. */
export const html = htm.bind(h);

/* A run's page carries its document inline. The raw template, served over HTTP from
   templates/ for work on the page, still holds the marker and fetches data.json (copied
   there from a run) instead. */
async function loadDocument() {
  const inline = document.getElementById("run-data").textContent;
  if (!inline.startsWith("/*")) return JSON.parse(inline);
  const response = await fetch("data.json");
  if (!response.ok) throw new Error("data.json: not found beside the template (copy one from a run)");
  return response.json();
}
export const DATA = await loadDocument();

/* ---------- run settings ---------- */
export const P = DATA.run.params;
export const T = P.consensus_coverage_threshold;
export const FACETS = DATA.run.facets || [];
export const FK = DATA.run.facet_kind || "none";
export const IS_ILLUMINA = DATA.run.data_type === "illumina";
export const PLATFORM = IS_ILLUMINA ? "Illumina" : "Nanopore";
export const FACET_LABEL = FK === "segment" ? "Segment" : FK === "virus" ? "Virus" : null;
export const FACET_PLURAL = FK === "segment" ? "Segments" : "Viruses";
export const ROWS_NOUN = FK === "segment" ? "sample–segment rows" : "samples";
export const version = DATA.run.version || "—";

/* ---------- colours: read from the stylesheet, so they follow the theme ---------- */
export const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
export const isDark = () => css("--is-dark") === "1";
/* facets take the categorical slots in run order; an unknown facet is muted */
const facetSlot = (f) => { const i = FACETS.indexOf(f); return i < 0 ? null : `--cat-${(i % 8) + 1}`; };
export const facetColor = (f) => { const s = facetSlot(f); return s ? `var(${s})` : "var(--muted)"; };

/* ---------- formatting ---------- */
export const fmt = (n) => n === null || n === undefined ? "—" : n.toLocaleString("en-US");
export const fmt1 = (n) => n === null || n === undefined ? "—" : n.toLocaleString("en-US", { maximumFractionDigits: 1 });
export const pct = (n) => n === null || n === undefined ? "—" : n.toFixed(1) + "%";
/* a 0–1 fraction as a percentage */
export const frac = (n) => n === null || n === undefined ? "—" : pct(100 * n);
export const pctOf = (a, b) => a === null || a === undefined || !b ? "—" : (100 * a / b).toFixed(1) + "%";
export const compact = (n) => n >= 1e6 ? (n / 1e6).toFixed(n >= 1e7 ? 0 : 1) + "M"
                            : n >= 1e4 ? Math.round(n / 1e3) + "K" : n.toLocaleString("en-US");
export const niceStep = (span) => { const p = 10 ** Math.floor(Math.log10(Math.max(1, span))); const r = span / p;
  return (r >= 5 ? 5 : r >= 2 ? 2 : 1) * p; };
export const short = (s, n = 18) => s && s.length > n ? s.slice(0, n - 1) + "…" : s;

/* ---------- derived run figures ---------- */
export const S = DATA.samples;
export const IDS = [...new Set(S.map((d) => d.id))];
const rowsOf = (id) => S.filter((d) => d.id === id);
export const anyRow = (id, pred) => rowsOf(id).some(pred);
const firstRow = (id) => rowsOf(id)[0];
export const isComplete = (d) => d.cov_min !== null && d.cov_min >= T;
export const gradeOf = (d) => d.viralqc && d.viralqc.grade ? String(d.viralqc.grade).toUpperCase() : null;
export const isGoodGrade = (d) => ["A", "B"].includes(gradeOf(d));
export const HAS_VQC = S.some((d) => d.viralqc !== null);
export const withCov = S.filter((d) => d.cov_min !== null);
export const complete = S.filter(isComplete);
export const goodGrade = S.filter(isGoodGrade);
export const totalReads = IDS.reduce((a, id) => a + (firstRow(id).total_reads || 0), 0);
const median = (a) => { if (!a.length) return null; const b = [...a].sort((x, y) => x - y), m = b.length >> 1;
  return b.length % 2 ? b[m] : (b[m - 1] + b[m]) / 2; };
export const medDepth = median(S.map((d) => d.median_depth).filter((v) => v !== null));
export const medMapped = median(S.map((d) => d.pct_mapped).filter((v) => v !== null));
export const facetCount = FACETS.map((f) => [f, S.filter((d) => d.facet === f).length]).filter(([, n]) => n);
/* rows with coverage, best first: the ranked chart's order and the first row selected */
export const ranked = withCov.slice().sort((a, b) => b.cov_min - a.cov_min);
/* tooltip rows naming a row's facet, when the run has facets */
export const facetRows = (d) => FACET_LABEL && d.facet ? [[FACET_LABEL, d.facet]] : [];

/* ---------- status: pill class + human wording; unknown statuses show verbatim ---------- */
const STATUS = {
  ok:               { cls: "good",     text: "Assembled" },
  no_mapped_reads:  { cls: "warning",  text: "No mapped reads" },
  empty_consensus:  { cls: "critical", text: "Empty consensus" },
  viralqc_failed:   { cls: "warning",  text: "viralQC failed" },
  viralqc_partial:  { cls: "warning",  text: "viralQC partial" },
  viralqc_skipped:  { cls: "neutral",  text: "viralQC skipped" },
  viralqc_missing:  { cls: "neutral",  text: "viralQC missing" },
  missing_stats:    { cls: "critical", text: "Missing stats" },
};
export const statusMeta = (k) => STATUS[k] || { cls: "neutral", text: k };
