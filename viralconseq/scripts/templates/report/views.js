"use strict";
/* What each part of report.html shows, worked out once from the document: the template
   reads these by name and keeps only the markup. */

/* ---------- KPI tiles ---------- */
const KPI_TILES = [
  { label: "Samples processed", value: IDS.length, foot: `${fmt(totalReads)} reads sequenced` },
  { label: `≥${T}% covered at ≥${P.min_depth}×`, value: complete.length,
    foot: `${(100 * complete.length / (S.length || 1)).toFixed(0)}% of ${S.length} ${ROWS_NOUN}` },
  HAS_VQC ? { label: "viralQC grade A or B", value: goodGrade.length,
    foot: `${S.filter((d) => gradeOf(d)).length} graded genomes` } : null,
  { label: "Median depth", value: fmt1(medDepth), unit: "×", foot: "across rows with depth" },
  FK !== "none" ? { label: `${FACET_PLURAL} detected`, value: facetCount.length,
    foot: facetCount.map(([f, n]) => `${f} ${n}`).join(" · ") || "none" } : null,
  { label: "Median % mapped", value: fmt1(medMapped), unit: "%", foot: "of reads, per row" },
].filter(Boolean);

/* ---------- funnel stages: unique samples surviving each step ---------- */
const FUNNEL_STAGES = [
  { label: "Samples with reads", n: IDS.filter((id) => anyRow(id, (d) => (d.total_reads || 0) > 0)).length,
    note: `${fmt(totalReads)} reads`, share: false },
  IS_ILLUMINA ? { label: "QC-passed", n: IDS.filter((id) => anyRow(id, (d) => (d.qc_passed_reads || 0) > 0)).length,
    note: "reads surviving fastp", share: true } : null,
  { label: "Mapped", n: IDS.filter((id) => anyRow(id, (d) => (d.mapped_reads || 0) > 0)).length,
    note: "at least one read on the reference", share: true },
  { label: `Covered ≥${T}% at ≥${P.min_depth}×`, n: new Set(complete.map((d) => d.id)).size,
    note: FK === "segment" ? "samples with at least one qualifying segment" : "in the filtered multi-FASTA", share: true },
  HAS_VQC ? { label: "viralQC grade A or B", n: new Set(goodGrade.map((d) => d.id)).size,
    note: "samples with at least one A/B genome", share: true } : null,
].filter(Boolean);

/* ---------- sample table ---------- */
const cladeText = (d) => d.viralqc ? [d.viralqc.clade, d.viralqc.lineage].filter(Boolean).join(" · ") : "";
const dash = `<span class="dim">—</span>`;
/* a numeric cell, dimmed when the value is missing */
const numCell = (v, text) => ({ cls: v === null ? "num dim" : "num", html: esc(text) });
function coverageCell(d) {
  if (d.cov_min === null) return { cls: "num", html: dash };
  const fill = `<span class="track"><i class="fill" style="width:${Math.max(2, d.cov_min)}%;background:${
    isComplete(d) ? "var(--ord-3)" : "var(--seq-mid)"}"></i></span>`;
  return { cls: "num", html: `<span class="meter">${fill}<span class="val">${d.cov_min.toFixed(1)}%</span></span>` };
}
/* k: sort key and column id; t: header; num: numeric sort (missing values last); get: the
   value sorted on; cell(d): the <td>'s class and inner HTML */
const COLS = [
  { k: "id", t: "Sample", num: false, get: (d) => d.id, cell: (d) => ({ cls: "id", html: esc(d.id) }) },
  { k: "status", t: "Status", num: false, get: (d) => d.status, cell: (d) => ({ cls: "", html: pill(d.status) }) },
  FACET_LABEL ? { k: "facet", t: FACET_LABEL, num: false, get: (d) => d.facet || "",
    cell: (d) => ({ cls: "", html: d.facet ? facetChip(d.facet) : dash }) } : null,
  { k: "total_reads", t: "Reads", num: true, get: (d) => d.total_reads, cell: (d) => numCell(d.total_reads, fmt(d.total_reads)) },
  IS_ILLUMINA ? { k: "qc_passed_reads", t: "QC-passed", num: true, get: (d) => d.qc_passed_reads,
    cell: (d) => numCell(d.qc_passed_reads, fmt(d.qc_passed_reads)) } : null,
  { k: "mapped_reads", t: "Mapped", num: true, get: (d) => d.mapped_reads, cell: (d) => numCell(d.mapped_reads, fmt(d.mapped_reads)) },
  { k: "pct_mapped", t: "% mapped", num: true, get: (d) => d.pct_mapped, cell: (d) => numCell(d.pct_mapped, pct(d.pct_mapped)) },
  { k: "median_depth", t: "Median depth", num: true, get: (d) => d.median_depth,
    cell: (d) => numCell(d.median_depth, d.median_depth === null ? "—" : fmt1(d.median_depth) + "×") },
  { k: "cov_min", t: `Cov ≥${P.min_depth}×`, num: true, get: (d) => d.cov_min, cell: coverageCell },
  { k: "n_pct", t: "N %", num: true, get: (d) => d.consensus ? d.consensus.n_pct : null,
    cell: (d) => numCell(d.consensus ? d.consensus.n_pct : null, d.consensus ? pct(d.consensus.n_pct) : "—") },
  HAS_VQC ? { k: "grade", t: "Grade", num: false, get: (d) => gradeOf(d) || "",
    cell: (d) => ({ cls: "", html: gradePill(gradeOf(d)) }) } : null,
  HAS_VQC ? { k: "clade", t: "Clade / lineage", num: false, get: (d) => cladeText(d),
    cell: (d) => ({ cls: cladeText(d) ? "" : "dim", html: esc(cladeText(d) || "—") }) } : null,
].filter(Boolean);

const FILTERS = [
  { k: "all", t: "All", match: () => true },
  { k: "complete", t: `≥${T}% covered`, match: isComplete },
  HAS_VQC ? { k: "grade", t: "Grade A/B", match: isGoodGrade } : null,
  { k: "flagged", t: "Flagged", match: (d) => d.status !== "ok" },
].filter(Boolean).map((f) => ({ ...f, n: S.filter(f.match).length }));

const FILTER_PLACEHOLDER = `Filter by sample, ${FACET_LABEL ? FACET_LABEL.toLowerCase() + ", " : ""}virus, clade or status`;

/* What the table shows for a view ({ filter, query, sortKey, sortDir }), in order. The
   detail panel's Previous/Next walk the same list. */
function visibleRows(view) {
  const f = FILTERS.find((x) => x.k === view.filter);
  const q = view.query.trim().toLowerCase();
  const rows = S.filter((d) => f.match(d) &&
    (!q || [d.id, d.facet, d.viralqc && d.viralqc.virus, d.viralqc && d.viralqc.clade,
            d.viralqc && d.viralqc.lineage, d.status, statusMeta(d.status).text]
      .some((v) => v && String(v).toLowerCase().includes(q))));
  const col = COLS.find((c) => c.k === view.sortKey) || COLS[0];
  return rows.sort((a, b) => {
    let x = col.get(a), y = col.get(b);
    if (x === null || x === undefined || x === "") x = col.num ? -Infinity : "￿";
    if (y === null || y === undefined || y === "") y = col.num ? -Infinity : "￿";
    if (x < y) return -view.sortDir;
    if (x > y) return view.sortDir;
    return a.key < b.key ? -1 : 1;
  });
}

/* ---------- sample detail ---------- */
const extentOf = (d) => DATA.run.extents[d.facet || "*"] || null;
const scoreNote = (vq) => vq && vq.score !== null ? ` · ${fmt1(vq.score)}/24` : "";
/* statuses the detail panel explains in words; anything else is shown verbatim */
const VIRALQC_STATUSES = ["viralqc_failed", "viralqc_partial", "viralqc_skipped", "viralqc_missing"];
const EXPLAINED_STATUSES = ["ok", "no_mapped_reads", "empty_consensus", "missing_stats", ...VIRALQC_STATUSES];

/* Stat strips: k is the label, v the value (HTML, escaped), n an optional note. */
function headlineStats(d) {
  const items = [{ k: "Reads", v: fmt(d.total_reads) }];
  if (IS_ILLUMINA) items.push({ k: "QC-passed", v: fmt(d.qc_passed_reads),
    n: d.total_reads ? `${pctOf(d.qc_passed_reads, d.total_reads)} of reads` : "" });
  items.push({ k: "Mapped", v: fmt(d.mapped_reads), n: `${pct(d.pct_mapped)} of reads` });
  items.push({ k: "Depth", v: fmt1(d.median_depth) + "<small>× median</small>",
    n: `mean ${fmt1(d.mean_depth)}×, zeros included` });
  items.push({ k: `Covered at ≥${P.min_depth}×`, v: pct(d.cov_min), n: `10×: ${pct(d.cov10)} · 100×: ${pct(d.cov100)}` });
  if (d.consensus) items.push({ k: "Consensus", v: fmt(d.consensus.length) + "<small> bp</small>",
    n: `${fmt(d.consensus.n)} N (${pct(d.consensus.n_pct)})` });
  if (d.variants) items.push({ k: "Variant records", v: fmt(d.variants.records), n: "lines in the final VCF" });
  if (d.isnv_count !== null) items.push({ k: "iSNVs", v: fmt(d.isnv_count),
    n: P.af_threshold !== null ? `AF ≥ ${P.af_threshold}` : "" });
  if (d.clair3_model) items.push({ k: "Clair3 model", v: `<small>${esc(d.clair3_model)}</small>` });
  return items;
}

function fastpStats(f) {
  return [
    { k: "Reads in", v: fmt(f.total_reads_before) },
    { k: "Reads out", v: fmt(f.total_reads_after), n: `${pctOf(f.total_reads_after, f.total_reads_before)} kept` },
    { k: "Q30 rate", v: frac(f.q30_after), n: f.q30_before === null || f.q30_before === undefined ? "" : `before: ${frac(f.q30_before)}` },
    { k: "Duplication", v: frac(f.duplication) },
    { k: "Mean length", v: `${fmt(f.read1_mean_length)}<small> / ${fmt(f.read2_mean_length)} bp</small>`, n: "R1 / R2 after trimming" },
    { k: "Low quality", v: fmt(f.low_quality), n: "reads removed" },
    { k: "Too short", v: fmt(f.too_short), n: "reads removed" },
    { k: "Adapter trimmed", v: fmt(f.adapter_trimmed), n: "reads" },
  ];
}

function variantCallingStats(d) {
  return [
    { k: "Clair3 model", v: `<small>${esc(d.clair3_model || P.clair3_model || "—")}</small>` },
    { k: "MAPQ floor", v: fmt(P.minimum_map_quality) },
    { k: "Allele fraction", v: P.af_threshold === null ? "—" : String(P.af_threshold) },
    { k: "Min read length", v: P.minimum_length === null ? "—" : fmt(P.minimum_length) + "<small> bp</small>" },
    { k: "Variant records", v: d.variants ? fmt(d.variants.records) : "—" },
  ];
}

/* viralQC's per-check verdicts that are present, with readable names */
const viralqcChecks = (vq) => Object.entries(vq.checks || {}).filter(([, v]) => v)
  .map(([k, v]) => [k.replace(/([A-Z])/g, " $1").toLowerCase(), v]);

/* ---------- provenance ---------- */
/* versions longer than this go to the fine print instead of the grid */
const LONG = 26;
const MANIFEST = DATA.run.manifest;
function viralqcLine(vs) {
  if (!vs) return null;
  return [vs.status, vs.exit_code !== undefined ? `exit ${vs.exit_code}` : null,
    vs.result_rows !== undefined && vs.input_records !== undefined ? `${vs.result_rows}/${vs.input_records} records` : null]
    .filter((v) => v !== null && v !== undefined).join(", ");
}
const PROVENANCE_ROWS = [
  ...DATA.run.versions.filter(([, v]) => v.length <= LONG),
  ["platform", PLATFORM],
  ["minimum depth", `${P.min_depth}×`],
  ["consensus coverage threshold", `${T}% at ≥${P.min_depth}×`],
  ["primer scheme", P.primer_scheme],
  [IS_ILLUMINA ? "allele-frequency threshold" : "allele-fraction threshold (ALT / (REF + ALT))", P.af_threshold],
  ["minimum read length", P.minimum_length === null ? null : `${fmt(P.minimum_length)} bp`],
  ["MAPQ floor", P.minimum_map_quality],
  ["Clair3 model", IS_ILLUMINA ? null : P.clair3_model],
  ["viralQC", viralqcLine(DATA.run.viralqc_status)],
  ["run created (UTC)", MANIFEST ? MANIFEST.created_utc : null],
  ["run finished (UTC)", MANIFEST ? MANIFEST.finished_utc : null],
  ["run status", MANIFEST ? MANIFEST.status : null],
  ["config sha256", MANIFEST && MANIFEST.config_sha256 ? MANIFEST.config_sha256.slice(0, 12) : null],
].filter(([, v]) => v !== null && v !== undefined && v !== "").map(([k, v]) => [k, String(v)]);
const LONG_VERSIONS = DATA.run.versions.filter(([, v]) => v.length > LONG);
const INPUT_LINES = (MANIFEST && MANIFEST.inputs || []).map((i) => ({
  sample: i.sample, file: i.path ? String(i.path).split("/").pop() : "—",
  rest: ` · sha256 ${i.sha256 || "—"}${i.size_bytes != null ? ` · ${compact(i.size_bytes)} B` : ""}`,
}));
