"use strict";

const KPI_TILES = [
  {
    label: "Samples processed",
    value: SAMPLE_IDS.length,
    foot: `${formatNumber(TOTAL_READS)} reads sequenced`,
  },
  {
    label: `≥${COVERAGE_THRESHOLD}% covered at ≥${MIN_DEPTH}×`,
    value: COMPLETE_ROWS.length,
    foot: `${(100 * COMPLETE_ROWS.length / (ROWS.length || 1)).toFixed(0)}% of ${ROWS.length} ${ROWS_NOUN}`,
  },
  HAS_VIRALQC && {
    label: "viralQC grade A or B",
    value: GOOD_GRADE_ROWS.length,
    foot: `${ROWS.filter(gradeOf).length} graded genomes`,
  },
  {
    label: "Median depth",
    value: formatDecimal(MEDIAN_DEPTH),
    unit: "×",
    foot: "across rows with depth",
  },
  FACET_KIND !== "none" && {
    label: `${FACET_PLURAL} detected`,
    value: FACET_COUNTS.length,
    foot: FACET_COUNTS.map(([facet, count]) => `${facet} ${count}`).join(" · ") || "none",
  },
  {
    label: "Median % mapped",
    value: formatDecimal(MEDIAN_PCT_MAPPED),
    unit: "%",
    foot: "of reads, per row",
  },
].filter(Boolean);

const FUNNEL_STAGES = [
  {
    label: "Samples with reads",
    count: countSamplesWithAnyRow((row) => row.total_reads > 0),
    note: `${formatNumber(TOTAL_READS)} reads`,
    showsShare: false,
  },
  IS_ILLUMINA && {
    label: "QC-passed",
    count: countSamplesWithAnyRow((row) => row.qc_passed_reads > 0),
    note: "reads surviving fastp",
    showsShare: true,
  },
  {
    label: "Mapped",
    count: countSamplesWithAnyRow((row) => row.mapped_reads > 0),
    note: "at least one read on the reference",
    showsShare: true,
  },
  {
    label: `Covered ≥${COVERAGE_THRESHOLD}% at ≥${MIN_DEPTH}×`,
    count: countSamples(COMPLETE_ROWS),
    note: IS_SEGMENTED ? "samples with at least one qualifying segment" : "in the filtered multi-FASTA",
    showsShare: true,
  },
  HAS_VIRALQC && {
    label: "viralQC grade A or B",
    count: countSamples(GOOD_GRADE_ROWS),
    note: "samples with at least one A/B genome",
    showsShare: true,
  },
].filter(Boolean);

const FUNNEL_NOTE = IS_SEGMENTED
  ? "Each stage counts unique samples. A sample with several segments contributes once here and one row per segment in the table."
  : "Each stage counts unique samples.";

const DIM_DASH_HTML = `<span class="dim">${MISSING}</span>`;

function pillHtml(tone, text) {
  return `<span class="pill ${tone}"><i class="dot"></i>${escapeHtml(text)}</span>`;
}

function statusPillHtml(status) {
  const { tone, text } = statusDisplay(status);
  return pillHtml(tone, text);
}

function gradePillHtml(grade, suffix = "") {
  if (!grade) return DIM_DASH_HTML;
  return `<span class="pill grade grade-${escapeHtml(grade.toLowerCase())}"><i class="dot"></i>${escapeHtml(grade)}${suffix}</span>`;
}

const VERDICT_TONES = { good: "good", pass: "good", mediocre: "warning", warn: "warning", bad: "critical", fail: "critical" };

function verdictPillHtml(verdict) {
  const tone = VERDICT_TONES[(verdict || "").toLowerCase()] || "neutral";
  return pillHtml(tone, verdict || MISSING);
}

function facetChipHtml(facet) {
  if (!facet) return "";
  return `<span class="chip"><i class="sw" style="background:${facetColor(facet)}"></i>${escapeHtml(facet)}</span>`;
}

function cladeText(row) {
  if (!row.viralqc) return "";
  return [row.viralqc.clade, row.viralqc.lineage].filter(Boolean).join(" · ");
}

function numericCell(value, text) {
  return { className: value === null ? "num dim" : "num", html: escapeHtml(text) };
}

function coverageCell(row) {
  if (row.cov_min === null) return { className: "num", html: DIM_DASH_HTML };
  const width = Math.max(2, row.cov_min);
  const colour = isComplete(row) ? "var(--ord-3)" : "var(--seq-mid)";
  const meter = `<span class="track"><i class="fill" style="width:${width}%;background:${colour}"></i></span>`;
  return { className: "num", html: `<span class="meter">${meter}<span class="val">${row.cov_min.toFixed(1)}%</span></span>` };
}

const TABLE_COLUMNS = [
  {
    key: "id", title: "Sample", numeric: false,
    sortValue: (row) => row.id,
    cell: (row) => ({ className: "id", html: escapeHtml(row.id) }),
  },
  {
    key: "status", title: "Status", numeric: false,
    sortValue: (row) => row.status,
    cell: (row) => ({ className: "", html: statusPillHtml(row.status) }),
  },
  FACET_LABEL && {
    key: "facet", title: FACET_LABEL, numeric: false,
    sortValue: (row) => row.facet || "",
    cell: (row) => ({ className: "", html: row.facet ? facetChipHtml(row.facet) : DIM_DASH_HTML }),
  },
  {
    key: "total_reads", title: "Reads", numeric: true,
    sortValue: (row) => row.total_reads,
    cell: (row) => numericCell(row.total_reads, formatNumber(row.total_reads)),
  },
  IS_ILLUMINA && {
    key: "qc_passed_reads", title: "QC-passed", numeric: true,
    sortValue: (row) => row.qc_passed_reads,
    cell: (row) => numericCell(row.qc_passed_reads, formatNumber(row.qc_passed_reads)),
  },
  {
    key: "mapped_reads", title: "Mapped", numeric: true,
    sortValue: (row) => row.mapped_reads,
    cell: (row) => numericCell(row.mapped_reads, formatNumber(row.mapped_reads)),
  },
  {
    key: "pct_mapped", title: "% mapped", numeric: true,
    sortValue: (row) => row.pct_mapped,
    cell: (row) => numericCell(row.pct_mapped, formatPercent(row.pct_mapped)),
  },
  {
    key: "median_depth", title: "Median depth", numeric: true,
    sortValue: (row) => row.median_depth,
    cell: (row) => numericCell(row.median_depth, row.median_depth === null ? MISSING : formatDecimal(row.median_depth) + "×"),
  },
  {
    key: "cov_min", title: `Cov ≥${MIN_DEPTH}×`, numeric: true,
    sortValue: (row) => row.cov_min,
    cell: coverageCell,
  },
  {
    key: "n_pct", title: "N %", numeric: true,
    sortValue: (row) => (row.consensus ? row.consensus.n_pct : null),
    cell: (row) => numericCell(row.consensus ? row.consensus.n_pct : null,
                               row.consensus ? formatPercent(row.consensus.n_pct) : MISSING),
  },
  HAS_VIRALQC && {
    key: "grade", title: "Grade", numeric: false,
    sortValue: (row) => gradeOf(row) || "",
    cell: (row) => ({ className: "", html: gradePillHtml(gradeOf(row)) }),
  },
  HAS_VIRALQC && {
    key: "clade", title: "Clade / lineage", numeric: false,
    sortValue: cladeText,
    cell: (row) => ({ className: cladeText(row) ? "" : "dim", html: escapeHtml(cladeText(row) || MISSING) }),
  },
].filter(Boolean);

const TABLE_FILTERS = [
  { key: "all", title: "All", matches: () => true },
  { key: "complete", title: `≥${COVERAGE_THRESHOLD}% covered`, matches: isComplete },
  HAS_VIRALQC && { key: "grade", title: "Grade A/B", matches: isGoodGrade },
  { key: "flagged", title: "Flagged", matches: (row) => row.status !== "ok" },
].filter(Boolean).map((filter) => ({ ...filter, count: ROWS.filter(filter.matches).length }));

const SEARCH_PLACEHOLDER =
  `Filter by sample, ${FACET_LABEL ? FACET_LABEL.toLowerCase() + ", " : ""}virus, clade or status`;

function rowsInView({ filter, query, sortKey, sortDirection }) {
  const { matches } = TABLE_FILTERS.find((candidate) => candidate.key === filter);
  const search = query.trim().toLowerCase();
  const column = TABLE_COLUMNS.find((candidate) => candidate.key === sortKey) || TABLE_COLUMNS[0];
  return ROWS
    .filter((row) => matches(row) && matchesSearch(row, search))
    .sort((a, b) => compareRows(a, b, column, sortDirection));
}

function matchesSearch(row, search) {
  if (!search) return true;
  const viralqc = row.viralqc || {};
  const searchable = [row.id, row.facet, viralqc.virus, viralqc.clade, viralqc.lineage, row.status, statusDisplay(row.status).text];
  return searchable.some((value) => value && String(value).toLowerCase().includes(search));
}

const SORTS_LAST_AS_TEXT = "￿";

function compareRows(a, b, column, direction) {
  const first = sortableValue(column.sortValue(a), column);
  const second = sortableValue(column.sortValue(b), column);
  if (first < second) return -direction;
  if (first > second) return direction;
  return a.key < b.key ? -1 : 1;
}

function sortableValue(value, column) {
  const isEmpty = value === null || value === undefined || value === "";
  if (!isEmpty) return value;
  return column.numeric ? -Infinity : SORTS_LAST_AS_TEXT;
}

const VIRALQC_FAILURE_STATUSES = ["viralqc_failed", "viralqc_partial", "viralqc_skipped", "viralqc_missing"];
const EXPLAINED_STATUSES = ["ok", "no_mapped_reads", "empty_consensus", "missing_stats", ...VIRALQC_FAILURE_STATUSES];

function gradeScoreSuffix(viralqc) {
  return viralqc && viralqc.score !== null ? ` · ${formatDecimal(viralqc.score)}/24` : "";
}

function headlineFigures(row) {
  return [
    { label: "Reads", value: formatNumber(row.total_reads) },
    IS_ILLUMINA && {
      label: "QC-passed",
      value: formatNumber(row.qc_passed_reads),
      note: row.total_reads ? `${formatShare(row.qc_passed_reads, row.total_reads)} of reads` : "",
    },
    { label: "Mapped", value: formatNumber(row.mapped_reads), note: `${formatPercent(row.pct_mapped)} of reads` },
    {
      label: "Depth",
      value: formatDecimal(row.median_depth) + "<small>× median</small>",
      note: `mean ${formatDecimal(row.mean_depth)}×, zeros included`,
    },
    {
      label: `Covered at ≥${MIN_DEPTH}×`,
      value: formatPercent(row.cov_min),
      note: `10×: ${formatPercent(row.cov10)} · 100×: ${formatPercent(row.cov100)}`,
    },
    row.consensus && {
      label: "Consensus",
      value: formatNumber(row.consensus.length) + "<small> bp</small>",
      note: `${formatNumber(row.consensus.n)} N (${formatPercent(row.consensus.n_pct)})`,
    },
    row.variants && { label: "Variant records", value: formatNumber(row.variants.records), note: "lines in the final VCF" },
    row.isnv_count !== null && {
      label: "iSNVs",
      value: formatNumber(row.isnv_count),
      note: PARAMS.af_threshold !== null ? `AF ≥ ${PARAMS.af_threshold}` : "",
    },
    row.clair3_model && { label: "Clair3 model", value: `<small>${escapeHtml(row.clair3_model)}</small>` },
  ].filter(Boolean);
}

function readQcFigures(fastp) {
  return [
    { label: "Reads in", value: formatNumber(fastp.total_reads_before) },
    {
      label: "Reads out",
      value: formatNumber(fastp.total_reads_after),
      note: `${formatShare(fastp.total_reads_after, fastp.total_reads_before)} kept`,
    },
    {
      label: "Q30 rate",
      value: formatFraction(fastp.q30_after),
      note: isMissing(fastp.q30_before) ? "" : `before: ${formatFraction(fastp.q30_before)}`,
    },
    { label: "Duplication", value: formatFraction(fastp.duplication) },
    {
      label: "Mean length",
      value: `${formatNumber(fastp.read1_mean_length)}<small> / ${formatNumber(fastp.read2_mean_length)} bp</small>`,
      note: "R1 / R2 after trimming",
    },
    { label: "Low quality", value: formatNumber(fastp.low_quality), note: "reads removed" },
    { label: "Too short", value: formatNumber(fastp.too_short), note: "reads removed" },
    { label: "Adapter trimmed", value: formatNumber(fastp.adapter_trimmed), note: "reads" },
  ];
}

function variantCallingFigures(row) {
  return [
    { label: "Clair3 model", value: `<small>${escapeHtml(row.clair3_model || PARAMS.clair3_model || MISSING)}</small>` },
    { label: "MAPQ floor", value: formatNumber(PARAMS.minimum_map_quality) },
    { label: "Allele fraction", value: PARAMS.af_threshold === null ? MISSING : String(PARAMS.af_threshold) },
    {
      label: "Min read length",
      value: PARAMS.minimum_length === null ? MISSING : formatNumber(PARAMS.minimum_length) + "<small> bp</small>",
    },
    { label: "Variant records", value: row.variants ? formatNumber(row.variants.records) : MISSING },
  ];
}

function mappedReadsFloorNote(row) {
  const floor = row.gate ? row.gate.minimum_mapped_reads : null;
  return isMissing(floor) ? "" : `, below the ${formatNumber(floor)}-read floor`;
}

function depthTraceNotes(row) {
  const contigs = coverageOf(row).contigs || [];
  const extent = primerExtentOf(row);
  const contigNote = contigs.length > 1 ? `; dashed lines separate the ${contigs.length} reference contigs` : "";
  const extentNote = extent
    ? `; the pale bands at each end lie outside the primer extent (${formatNumber(extent.outside_bp)} bp) `
      + "and describe scheme geometry, not a measured coverage limit"
    : "";
  return contigNote + extentNote;
}

function ampliconCaption(row) {
  const depths = ampliconDepthsOf(row);
  const below = depths.filter((depth) => depth < MIN_DEPTH).length;
  return `${below} of ${depths.length} amplicons below ${MIN_DEPTH}×`;
}

function alignmentCoverage(viralqc) {
  return viralqc.coverage !== null ? (100 * viralqc.coverage).toFixed(1) + "%" : MISSING;
}

function viralqcCheckVerdicts(viralqc) {
  const humanName = (check) => check.replace(/([A-Z])/g, " $1").toLowerCase();
  return Object.entries(viralqc.checks || {})
    .filter(([, verdict]) => verdict)
    .map(([check, verdict]) => [humanName(check), verdict]);
}

const MAX_GRID_VERSION_LENGTH = 26;
const MANIFEST = DATA.run.manifest || {};

function viralqcRunSummary(status) {
  if (!status) return null;
  const parts = [
    status.status,
    status.exit_code !== undefined ? `exit ${status.exit_code}` : null,
    status.result_rows !== undefined && status.input_records !== undefined
      ? `${status.result_rows}/${status.input_records} records`
      : null,
  ];
  return parts.filter((part) => !isMissing(part)).join(", ");
}

const PROVENANCE_ROWS = [
  ...DATA.run.versions.filter(([, version]) => version.length <= MAX_GRID_VERSION_LENGTH),
  ["platform", PLATFORM],
  ["minimum depth", `${MIN_DEPTH}×`],
  ["consensus coverage threshold", `${COVERAGE_THRESHOLD}% at ≥${MIN_DEPTH}×`],
  ["primer scheme", PARAMS.primer_scheme],
  [IS_ILLUMINA ? "allele-frequency threshold" : "allele-fraction threshold (ALT / (REF + ALT))", PARAMS.af_threshold],
  ["minimum read length", PARAMS.minimum_length === null ? null : `${formatNumber(PARAMS.minimum_length)} bp`],
  ["MAPQ floor", PARAMS.minimum_map_quality],
  ["Clair3 model", IS_ILLUMINA ? null : PARAMS.clair3_model],
  ["viralQC", viralqcRunSummary(DATA.run.viralqc_status)],
  ["run created (UTC)", MANIFEST.created_utc],
  ["run finished (UTC)", MANIFEST.finished_utc],
  ["run status", MANIFEST.status],
  ["config sha256", MANIFEST.config_sha256 ? MANIFEST.config_sha256.slice(0, 12) : null],
]
  .filter(([, value]) => !isMissing(value) && value !== "")
  .map(([name, value]) => [name, String(value)]);

const LONG_VERSIONS = DATA.run.versions.filter(([, version]) => version.length > MAX_GRID_VERSION_LENGTH);

const INPUT_FILES = (MANIFEST.inputs || []).map((input) => ({
  sample: input.sample,
  file: input.path ? String(input.path).split("/").pop() : MISSING,
  details: ` · sha256 ${input.sha256 || MISSING}${isMissing(input.size_bytes) ? "" : ` · ${formatCompact(input.size_bytes)} B`}`,
}));
