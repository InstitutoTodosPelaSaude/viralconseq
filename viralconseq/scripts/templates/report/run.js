"use strict";

const DATA = loadRunDocument();

function loadRunDocument() {
  const inline = document.getElementById("run-data").textContent;
  const isTemplateWithoutData = inline.startsWith("/*");
  return JSON.parse(isTemplateWithoutData ? readFileBesideTemplate("data.json") : inline);
}

// Synchronous on purpose: every script after this one reads DATA as it loads.
function readFileBesideTemplate(name) {
  const request = new XMLHttpRequest();
  request.open("GET", name, false);
  request.send();
  if (request.status !== 200) throw new Error(`${name}: not found beside the template (copy one from a run)`);
  return request.responseText;
}

const PARAMS = DATA.run.params;
const MIN_DEPTH = PARAMS.min_depth;
const COVERAGE_THRESHOLD = PARAMS.consensus_coverage_threshold;
const IS_ILLUMINA = DATA.run.data_type === "illumina";
const PLATFORM = IS_ILLUMINA ? "Illumina" : "Nanopore";
const PIPELINE_VERSION = DATA.run.version || "—";

const FACETS = DATA.run.facets || [];
const FACET_KIND = DATA.run.facet_kind || "none";
const IS_SEGMENTED = FACET_KIND === "segment";
const FACET_LABEL = { segment: "Segment", virus: "Virus" }[FACET_KIND] || null;
const FACET_PLURAL = IS_SEGMENTED ? "Segments" : "Viruses";
const ROWS_NOUN = IS_SEGMENTED ? "sample–segment rows" : "samples";

const ROWS = DATA.samples;
const SAMPLE_IDS = [...new Set(ROWS.map((row) => row.id))];
const HAS_VIRALQC = ROWS.some((row) => row.viralqc !== null);

function isComplete(row) {
  return row.cov_min !== null && row.cov_min >= COVERAGE_THRESHOLD;
}

function gradeOf(row) {
  return row.viralqc && row.viralqc.grade ? String(row.viralqc.grade).toUpperCase() : null;
}

function isGoodGrade(row) {
  return ["A", "B"].includes(gradeOf(row));
}

function countSamplesWithAnyRow(predicate) {
  return SAMPLE_IDS.filter((id) => ROWS.some((row) => row.id === id && predicate(row))).length;
}

function countSamples(rows) {
  return new Set(rows.map((row) => row.id)).size;
}

function median(values) {
  const sorted = values.filter((value) => value !== null).sort((a, b) => a - b);
  if (sorted.length === 0) return null;
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
}

const COMPLETE_ROWS = ROWS.filter(isComplete);
const GOOD_GRADE_ROWS = ROWS.filter(isGoodGrade);
const RANKED_ROWS = ROWS.filter((row) => row.cov_min !== null).sort((a, b) => b.cov_min - a.cov_min);
const TOTAL_READS = SAMPLE_IDS.reduce((sum, id) => sum + (ROWS.find((row) => row.id === id).total_reads || 0), 0);
const MEDIAN_DEPTH = median(ROWS.map((row) => row.median_depth));
const MEDIAN_PCT_MAPPED = median(ROWS.map((row) => row.pct_mapped));
const FACET_COUNTS = FACETS
  .map((facet) => [facet, ROWS.filter((row) => row.facet === facet).length])
  .filter(([, count]) => count > 0);

function coverageOf(row) {
  return DATA.coverage[row.key];
}

function ampliconDepthsOf(row) {
  return DATA.amplicons[row.key];
}

function primerExtentOf(row) {
  return DATA.run.extents[row.facet || "*"] || null;
}

function cssVar(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

function isDark() {
  return cssVar("--is-dark") === "1";
}

function facetColor(facet) {
  if (FACET_KIND === "none") return "var(--cat-1)";
  const slot = FACETS.indexOf(facet);
  return slot < 0 ? "var(--muted)" : `var(--cat-${(slot % 8) + 1})`;
}

const STATUS_DISPLAY = {
  ok:              { tone: "good",     text: "Assembled" },
  no_mapped_reads: { tone: "warning",  text: "No mapped reads" },
  empty_consensus: { tone: "critical", text: "Empty consensus" },
  viralqc_failed:  { tone: "warning",  text: "viralQC failed" },
  viralqc_partial: { tone: "warning",  text: "viralQC partial" },
  viralqc_skipped: { tone: "neutral",  text: "viralQC skipped" },
  viralqc_missing: { tone: "neutral",  text: "viralQC missing" },
  missing_stats:   { tone: "critical", text: "Missing stats" },
};

function statusDisplay(status) {
  return STATUS_DISPLAY[status] || { tone: "neutral", text: status };
}
