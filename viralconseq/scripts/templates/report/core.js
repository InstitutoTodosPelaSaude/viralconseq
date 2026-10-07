"use strict";
/* The run document, the figures derived from it, and the formatting every part of the
   page shares. A classic script: what it declares is global, so report.html's Alpine
   expressions and the other report/*.js read it by name. */

/* A run's page carries its document inline. The raw template, served over HTTP from
   templates/ for work on the page, still holds the marker and reads data.json (copied
   there from a run) instead; synchronously, because everything below reads DATA at load. */
function readSibling(name) {
  const req = new XMLHttpRequest();
  req.open("GET", name, false);
  req.send();
  if (req.status !== 200) throw new Error(`${name}: not found beside the template (copy one from a run)`);
  return req.responseText;
}
const INLINE = document.getElementById("run-data").textContent;
const DATA = JSON.parse(INLINE.startsWith("/*") ? readSibling("data.json") : INLINE);

/* ---------- run settings ---------- */
const P = DATA.run.params;
const T = P.consensus_coverage_threshold;
const FACETS = DATA.run.facets || [];
const FK = DATA.run.facet_kind || "none";
const IS_ILLUMINA = DATA.run.data_type === "illumina";
const PLATFORM = IS_ILLUMINA ? "Illumina" : "Nanopore";
const FACET_LABEL = FK === "segment" ? "Segment" : FK === "virus" ? "Virus" : null;
const FACET_PLURAL = FK === "segment" ? "Segments" : "Viruses";
const ROWS_NOUN = FK === "segment" ? "sample–segment rows" : "samples";
const version = DATA.run.version || "—";

/* ---------- colours: read from the stylesheet, so they follow the theme ---------- */
const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const isDark = () => css("--is-dark") === "1";
/* facets take the categorical slots in run order; an unknown facet is muted */
const facetSlot = (f) => { const i = FACETS.indexOf(f); return i < 0 ? null : `--cat-${(i % 8) + 1}`; };
const facetColor = (f) => { const s = facetSlot(f); return s ? `var(${s})` : "var(--muted)"; };

/* ---------- formatting ---------- */
const fmt = (n) => n === null || n === undefined ? "—" : n.toLocaleString("en-US");
const fmt1 = (n) => n === null || n === undefined ? "—" : n.toLocaleString("en-US", { maximumFractionDigits: 1 });
const pct = (n) => n === null || n === undefined ? "—" : n.toFixed(1) + "%";
/* a 0–1 fraction as a percentage */
const frac = (n) => n === null || n === undefined ? "—" : pct(100 * n);
const pctOf = (a, b) => a === null || a === undefined || !b ? "—" : (100 * a / b).toFixed(1) + "%";
const compact = (n) => n >= 1e6 ? (n / 1e6).toFixed(n >= 1e7 ? 0 : 1) + "M"
                     : n >= 1e4 ? Math.round(n / 1e3) + "K" : n.toLocaleString("en-US");
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const niceStep = (span) => { const p = 10 ** Math.floor(Math.log10(Math.max(1, span))); const r = span / p;
  return (r >= 5 ? 5 : r >= 2 ? 2 : 1) * p; };
const short = (s, n = 18) => s && s.length > n ? s.slice(0, n - 1) + "…" : s;

/* ---------- derived run figures ---------- */
const S = DATA.samples;
const IDS = [...new Set(S.map((d) => d.id))];
const rowsOf = (id) => S.filter((d) => d.id === id);
const anyRow = (id, pred) => rowsOf(id).some(pred);
const firstRow = (id) => rowsOf(id)[0];
const isComplete = (d) => d.cov_min !== null && d.cov_min >= T;
const gradeOf = (d) => d.viralqc && d.viralqc.grade ? String(d.viralqc.grade).toUpperCase() : null;
const isGoodGrade = (d) => ["A", "B"].includes(gradeOf(d));
const HAS_VQC = S.some((d) => d.viralqc !== null);
const withCov = S.filter((d) => d.cov_min !== null);
const complete = S.filter(isComplete);
const goodGrade = S.filter(isGoodGrade);
const totalReads = IDS.reduce((a, id) => a + (firstRow(id).total_reads || 0), 0);
const median = (a) => { if (!a.length) return null; const b = [...a].sort((x, y) => x - y), m = b.length >> 1;
  return b.length % 2 ? b[m] : (b[m - 1] + b[m]) / 2; };
const medDepth = median(S.map((d) => d.median_depth).filter((v) => v !== null));
const medMapped = median(S.map((d) => d.pct_mapped).filter((v) => v !== null));
const facetCount = FACETS.map((f) => [f, S.filter((d) => d.facet === f).length]).filter(([, n]) => n);
/* rows with coverage, best first: the ranked chart's order and the first row selected */
const ranked = withCov.slice().sort((a, b) => b.cov_min - a.cov_min);
/* tooltip rows naming a row's facet, when the run has facets */
const facetRows = (d) => FACET_LABEL && d.facet ? [[FACET_LABEL, esc(d.facet)]] : [];

/* ---------- status, pills and chips ---------- */
/* status -> pill class + human wording; unknown statuses show verbatim */
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
const statusMeta = (k) => STATUS[k] || { cls: "neutral", text: k };

/* Small HTML snippets the template drops in with x-html; every value goes through esc(). */
const pill = (k) => { const m = statusMeta(k); return `<span class="pill ${m.cls}"><i class="dot"></i>${esc(m.text)}</span>`; };
const gradePill = (g, extra = "") => g
  ? `<span class="pill grade grade-${esc(String(g).toLowerCase())}"><i class="dot"></i>${esc(g)}${extra}</span>`
  : `<span class="dim">—</span>`;
/* viralQC's overall and per-check verdicts */
function qcPill(status) {
  const s = (status || "").toLowerCase();
  const cls = s === "good" || s === "pass" ? "good" : s === "mediocre" || s === "warn" ? "warning"
            : s === "bad" || s === "fail" ? "critical" : "neutral";
  return `<span class="pill ${cls}"><i class="dot"></i>${esc(status || "—")}</span>`;
}
const facetChip = (f) => f
  ? `<span class="chip"><i class="sw" style="background:${facetColor(f)}"></i>${esc(f)}</span>` : "";
