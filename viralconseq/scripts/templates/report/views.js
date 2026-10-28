"use strict";

function fact(label, value) {
  return { label, value };
}

const GRADED_ROWS = ROWS.filter(gradeOf);

// The run at a glance, in the cards the sample detail uses: label/value rows, the explanation in each tooltip.
const RUN_FACT_GROUPS = [
  {
    label: "Samples",
    tip: "Samples in the sample sheet and the reads sequenced across them."
      + (FACET_KIND === "none" ? "" : ` ${FACET_PLURAL} detected, with the rows each one has: `
        + `${FACET_COUNTS.map(([facet, count]) => `${facet} ${count}`).join(", ") || "none"}.`),
    facts: [
      fact("Processed", formatNumber(SAMPLE_IDS.length)),
      fact("Reads sequenced", formatNumber(TOTAL_READS)),
      FACET_KIND !== "none" && fact(`${FACET_PLURAL} detected`, formatNumber(FACET_COUNTS.length)),
    ].filter(Boolean),
  },
  {
    label: "Genomes",
    tip: `Of the ${ROWS.length} ${ROWS_NOUN}: those with ≥${COVERAGE_THRESHOLD}% of reference positions at ≥${MIN_DEPTH}×, which enter the filtered multi-FASTA`
      + (HAS_VIRALQC ? "; and the genomes viralQC graded A or B, out of those it graded." : "."),
    facts: [
      fact(`Covered ≥${COVERAGE_THRESHOLD}% at ≥${MIN_DEPTH}×`,
        `${formatNumber(COMPLETE_ROWS.length)} of ${formatNumber(ROWS.length)} · ${(100 * COMPLETE_ROWS.length / (ROWS.length || 1)).toFixed(0)}%`),
      HAS_VIRALQC && fact("viralQC grade A or B", `${formatNumber(GOOD_GRADE_ROWS.length)} of ${formatNumber(GRADED_ROWS.length)}`),
    ].filter(Boolean),
  },
  {
    label: "Depth & mapping",
    tip: "Medians across rows: the median depth over the rows that have depth, and the share of each row's reads that mapped to the reference.",
    facts: [
      fact("Median depth", isMissing(MEDIAN_DEPTH) ? MISSING : `${formatDecimal(MEDIAN_DEPTH)}×`),
      fact("Median % mapped", isMissing(MEDIAN_PCT_MAPPED) ? MISSING : `${formatDecimal(MEDIAN_PCT_MAPPED)}%`),
    ],
  },
];

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

// Icons from the Nimbus design (vuesax linear arrows, the filled info badge), drawn in currentColor so they follow the theme.
const ICON_ARROW_LEFT = '<svg class="icon" width="16" height="16" viewBox="0 0 16 16" fill="none" aria-hidden="true"><path d="M10 13.28L5.65333 8.93333C5.14 8.42 5.14 7.58 5.65333 7.06667L10 2.72" stroke-width="1.5" stroke-miterlimit="10" stroke-linecap="round" stroke-linejoin="round" stroke="currentColor"/></svg>';
const ICON_ARROW_RIGHT = '<svg class="icon" width="16" height="16" viewBox="0 0 16 16" fill="none" aria-hidden="true"><path d="M5.94 13.28L10.2867 8.93333C10.8 8.42 10.8 7.58 10.2867 7.06667L5.94 2.72" stroke-width="1.5" stroke-miterlimit="10" stroke-linecap="round" stroke-linejoin="round" stroke="currentColor"/></svg>';
const ICON_CLOSE = '<svg class="icon" width="16" height="16" viewBox="0 0 16 16" fill="none" aria-hidden="true"><path d="M4 12L12 4M12 12L4 4" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/></svg>';
const ICON_SUN = '<svg class="icon" width="16" height="16" viewBox="0 0 16 16" fill="none" aria-hidden="true"><circle cx="8" cy="8" r="3" stroke="currentColor" stroke-width="1.5"/><path d="M8 1.5V2.5M8 13.5V14.5M1.5 8H2.5M13.5 8H14.5M3.4 3.4L4.1 4.1M11.9 11.9L12.6 12.6M3.4 12.6L4.1 11.9M11.9 4.1L12.6 3.4" stroke="currentColor" stroke-width="1.5" stroke-linecap="round"/></svg>';
const ICON_MOON = '<svg class="icon" width="16" height="16" viewBox="0 0 16 16" fill="none" aria-hidden="true"><path d="M13.6 10.1A5.9 5.9 0 0 1 5.9 2.4 5.9 5.9 0 1 0 13.6 10.1Z" stroke="currentColor" stroke-width="1.5" stroke-linejoin="round"/></svg>';
const ICON_SEARCH = '<svg class="icon" width="18" height="18" viewBox="0 0 16 16" fill="none" aria-hidden="true"><path d="M7.66665 14C11.1644 14 14 11.1644 14 7.66665C14 4.16884 11.1644 1.33331 7.66665 1.33331C4.16884 1.33331 1.33331 4.16884 1.33331 7.66665C1.33331 11.1644 4.16884 14 7.66665 14Z" stroke="currentColor" stroke-linecap="round" stroke-linejoin="round"/><path d="M14.6666 14.6666L13.3333 13.3333" stroke="currentColor" stroke-linecap="round" stroke-linejoin="round"/></svg>';
const ICON_INFO = '<svg class="icon" width="18" height="18" viewBox="0 0 18 18" aria-hidden="true"><circle cx="9" cy="9" r="9" fill="currentColor"/><circle cx="9" cy="5.45" r="1.1" fill="#fff"/><path d="M9 8.15V12.75" stroke="#fff" stroke-width="1.8" stroke-linecap="round"/></svg>';
const ICON_CHEVRON_DOWN = '<svg class="icon" width="16" height="16" viewBox="0 0 16 16" fill="none" aria-hidden="true"><path d="M13.28 5.96667L8.93333 10.3133C8 11.2467 8 11.2467 7.06667 10.3133L2.72 5.96667" stroke-miterlimit="10" stroke-linecap="round" stroke-linejoin="round" stroke="currentColor"/></svg>';
const ICON_ARROW_DOWN = '<svg class="icon" width="14" height="14" viewBox="0 0 14 14" fill="none" aria-hidden="true"><path d="M7 2.33333L7 11.6667M3.5 8.16667L7 11.6667L10.5 8.16667" stroke-linecap="round" stroke-linejoin="round" stroke="currentColor"/></svg>';

// Centred below its icon, slid sideways to stay inside the window; the arrow keeps pointing at the icon.
const TOOLTIP_MARGIN = 16;

function placeTooltip(info) {
  const tip = info.querySelector(".tip");
  const icon = info.getBoundingClientRect();
  const iconCentre = icon.left + icon.width / 2;
  const viewportWidth = document.documentElement.clientWidth;
  const left = Math.min(Math.max(iconCentre - tip.offsetWidth / 2, TOOLTIP_MARGIN), viewportWidth - TOOLTIP_MARGIN - tip.offsetWidth);
  tip.style.left = `${left - icon.left}px`;
  tip.style.setProperty("--arrow-x", `${iconCentre - left}px`);
}

// A figure cut short by its card shows its full text in a tooltip; a figure below its threshold always does.
function revealCutValue(figure) {
  const value = figure.querySelector(".value");
  figure.dataset.tip = value.scrollWidth > value.clientWidth || figure.classList.contains("low") ? "on" : "off";
  placeTooltip(figure);
}

// The tabs' sliding pill sits under the pressed tab, on its row when the tabs wrap; it only animates
// once it has been placed.
function placeTabIndicator(tabs) {
  const pressed = tabs.querySelector('.tab[aria-pressed="true"]');
  if (!pressed) return;
  tabs.style.setProperty("--indicator-x", `${pressed.offsetLeft}px`);
  tabs.style.setProperty("--indicator-y", `${pressed.offsetTop}px`);
  tabs.style.setProperty("--indicator-width", `${pressed.offsetWidth}px`);
  requestAnimationFrame(() => { tabs.dataset.ready = "true"; });
}

function pillHtml(tone, text, icon = "") {
  return `<span class="pill ${tone}">${icon}${escapeHtml(text)}</span>`;
}

function lowCoveragePillHtml(row) {
  return row.cov_min !== null && !isComplete(row) ? pillHtml("warning", "Low coverage", ICON_ARROW_DOWN) : "";
}

function statusPillHtml(status) {
  const { tone, text } = statusDisplay(status);
  return pillHtml(tone, text);
}

function gradePillHtml(grade, suffix = "") {
  if (!grade) return DIM_DASH_HTML;
  const badge = `<span class="badge grade-${escapeHtml(grade.toLowerCase())}">${escapeHtml(grade)}</span>`;
  return suffix ? `<span class="badge-wrap">${badge}<span class="score">${escapeHtml(suffix)}</span></span>` : badge;
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
  const value = `${row.cov_min.toFixed(1)}%`;
  if (isComplete(row)) return { className: "num", html: `<span class="cov">${value}</span>` };
  const title = `below the ${COVERAGE_THRESHOLD}% consensus threshold`;
  return { className: "num", html: `<span class="cov low" title="${title}">${value}${ICON_ARROW_DOWN}</span>` };
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

// The sample's figures as cards of label/value rows, one card per stage of the pipeline, the short Reads card
// last; each card's tooltip says what its figures mean. Run-wide settings (the read QC and variant-calling floors) sit with the stage they govern.
function factGroups(row) {
  const withShare = (count, share) => (isMissing(count) ? MISSING : `${formatNumber(count)} · ${share}`);
  const atLeast = (value, unit = "") => (isMissing(value) ? MISSING : `≥ ${typeof value === "number" ? formatNumber(value) : value}${unit}`);
  const fastp = row.fastp;
  return [
    IS_ILLUMINA && {
      label: "Read QC",
      tip: "fastp on the read pairs, before and after adapter and quality trimming; the mapped and depth figures are computed from the reads that passed.",
      empty: fastp ? null : "No fastp report for this sample.",
      facts: fastp ? [
        fact("Reads in", formatNumber(fastp.total_reads_before)),
        fact("Reads out", withShare(fastp.total_reads_after, `${formatShare(fastp.total_reads_after, fastp.total_reads_before)} kept`)),
        fact("Q30", isMissing(fastp.q30_before) ? formatFraction(fastp.q30_after) : `${formatFraction(fastp.q30_after)} · before ${formatFraction(fastp.q30_before)}`),
        fact("Duplication", formatFraction(fastp.duplication)),
        fact("Mean length R1 / R2", `${formatNumber(fastp.read1_mean_length)} / ${formatNumber(fastp.read2_mean_length)} bp`),
        fact("Low quality removed", formatNumber(fastp.low_quality)),
        fact("Too short removed", formatNumber(fastp.too_short)),
        fact("Adapter trimmed", formatNumber(fastp.adapter_trimmed)),
      ] : [],
    },
    {
      label: "Depth & coverage",
      tip: `Depth is read depth over every reference position, zeros included. Coverage is the share of reference positions at or above each depth; positions below ${MIN_DEPTH}× are masked to N in the consensus, and a genome enters the filtered multi-FASTA at ≥${COVERAGE_THRESHOLD}% covered at ${MIN_DEPTH}×.`,
      facts: [
        fact("Median depth", isMissing(row.median_depth) ? MISSING : `${formatDecimal(row.median_depth)}×`),
        fact("Mean depth", isMissing(row.mean_depth) ? MISSING : `${formatDecimal(row.mean_depth)}×`),
        {
          ...fact(`Covered ≥${MIN_DEPTH}×`, formatPercent(row.cov_min)),
          low: row.cov_min !== null && !isComplete(row),
          lowNote: `below the ${COVERAGE_THRESHOLD}% consensus threshold`,
        },
        fact("Covered ≥10×", formatPercent(row.cov10)),
        fact("Covered ≥100×", formatPercent(row.cov100)),
      ],
    },
    {
      label: "Consensus",
      tip: "The consensus sequence and its N count; variant records are the lines of the final VCF"
        + (row.isnv_count !== null && PARAMS.af_threshold !== null ? `; iSNVs are intra-host variants at allele frequency ≥ ${PARAMS.af_threshold}.` : "."),
      facts: [
        row.consensus && fact("Length", `${formatNumber(row.consensus.length)} bp`),
        row.consensus && fact("N", withShare(row.consensus.n, formatPercent(row.consensus.n_pct))),
        row.variants && fact("Variant records", formatNumber(row.variants.records)),
        row.isnv_count !== null && fact("iSNVs", formatNumber(row.isnv_count)),
      ].filter(Boolean),
      empty: row.consensus || row.variants || row.isnv_count !== null ? null : "No consensus figures for this row.",
    },
    !IS_ILLUMINA && {
      label: "Variant calling",
      tip: `Clair3, with the model shown, calls on reads passing the length and mapping-quality floors; the consensus applies variants whose allele fraction, ALT over REF plus ALT reads, is at or above the threshold shown, and masks positions below ${MIN_DEPTH}×. Every setting but the model is run-wide and also listed under Provenance.`,
      facts: [
        fact("Model", row.clair3_model || PARAMS.clair3_model || MISSING),
        fact("MAPQ", atLeast(PARAMS.minimum_map_quality)),
        fact("Allele fraction", atLeast(PARAMS.af_threshold === null ? null : String(PARAMS.af_threshold))),
        fact("Read length", atLeast(PARAMS.minimum_length, " bp")),
      ],
    },
    {
      label: "Reads",
      tip: IS_ILLUMINA
        ? "Reads in the sample's FASTQ files. QC-passed reads survived fastp trimming; mapped reads are QC-passed reads aligned to the reference, as a share of all reads."
        : "Reads in the sample's FASTQ file; mapped reads aligned to the reference, as a share of all reads.",
      facts: [
        fact("Total", formatNumber(row.total_reads)),
        IS_ILLUMINA && fact("QC-passed", withShare(row.qc_passed_reads, formatShare(row.qc_passed_reads, row.total_reads))),
        fact("Mapped", withShare(row.mapped_reads, formatPercent(row.pct_mapped))),
      ].filter(Boolean),
    },
  ].filter(Boolean);
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

// The segment table in the detail: every row of the selected sample, one per segment.
const SEGMENT_COLUMNS = [
  { title: "Segment", cell: (row) => ({ className: "strong", html: row.facet ? facetChipHtml(row.facet) : DIM_DASH_HTML }) },
  { title: "Status", cell: (row) => ({ className: "", html: statusPillHtml(row.status) }) },
  { title: "Mapped", cell: (row) => numericCell(row.mapped_reads, formatNumber(row.mapped_reads)) },
  {
    title: "Median depth",
    cell: (row) => numericCell(row.median_depth, row.median_depth === null ? MISSING : formatDecimal(row.median_depth) + "×"),
  },
  { title: `Coverage ≥${MIN_DEPTH}×`, cell: coverageCell },
  {
    title: "N %",
    cell: (row) => numericCell(row.consensus ? row.consensus.n_pct : null,
                               row.consensus ? formatPercent(row.consensus.n_pct) : MISSING),
  },
  HAS_VIRALQC && { title: "Quality", cell: (row) => ({ className: "", html: gradePillHtml(gradeOf(row)) }) },
].filter(Boolean);

function segmentRowsOf(row) {
  return IS_SEGMENTED ? ROWS.filter((candidate) => candidate.id === row.id) : [];
}

// viralQC's figures for the genome, in the detail's card style; its per-check statuses follow as a table.
function viralqcFactGroups(row) {
  const qc = row.viralqc;
  const grade = gradeOf(row);
  const score = gradeScoreSuffix(qc);
  return [
    {
      label: "Genome quality",
      tip: "viralQC's genome grade (A best, D worst) with its score out of 24, the overall QC status, the share of the reference the consensus aligned to, and the input sequence's status.",
      facts: [
        { label: "Grade", value: grade ? `${grade}${score}` : MISSING, html: grade ? gradePillHtml(grade, score) : null },
        { label: "Overall status", value: qc.overall || MISSING, html: qc.overall ? verdictPillHtml(qc.overall) : null },
        fact("Alignment coverage", alignmentCoverage(qc)),
        qc.input_status && fact("Input sequence", qc.input_status),
      ].filter(Boolean),
    },
    {
      label: "Classification",
      tip: "The virus, clade and lineage viralQC assigned, and the dataset it compared the consensus against.",
      facts: [
        fact("Virus", qc.virus || MISSING),
        qc.species && fact("Species", qc.species),
        fact("Clade", qc.clade || MISSING),
        fact("Lineage", qc.lineage || MISSING),
        fact("Dataset", qc.dataset || MISSING),
        qc.dataset_version && fact("Dataset version", qc.dataset_version),
      ].filter(Boolean),
    },
    {
      label: "Differences",
      tip: "What viralQC found in the consensus against the dataset's reference: substitutions, deletions and missing sites.",
      facts: [
        fact("Substitutions", formatNumber(qc.substitutions)),
        fact("Deletions", formatNumber(qc.deletions)),
        fact("Missing sites", formatNumber(qc.missing)),
      ],
    },
  ];
}

function viralqcCheckVerdicts(viralqc) {
  const humanName = (check) => {
    const words = check.replace(/([A-Z])/g, " $1").toLowerCase();
    return (words.charAt(0).toUpperCase() + words.slice(1)).replace(/\bsnp\b/i, "SNP");
  };
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
