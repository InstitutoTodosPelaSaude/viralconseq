/* The sample table: filter buttons, free-text filter, sortable columns. */
import {
  html, FACET_LABEL, HAS_VQC, IS_ILLUMINA, P, ROWS_NOUN, S, T,
  fmt, fmt1, gradeOf, isComplete, isGoodGrade, pct, statusMeta,
} from "report/core";
import { FacetChip, GradePill, Pill } from "report/ui";

const dash = html`<span class="dim">—</span>`;
const cladeText = (d) => d.viralqc ? [d.viralqc.clade, d.viralqc.lineage].filter(Boolean).join(" · ") : "";
/* a numeric cell, dimmed when the value is missing */
const numCell = (v, text) => html`<td class=${v === null ? "num dim" : "num"}>${text}</td>`;

function CoverageCell({ d }) {
  if (d.cov_min === null) return html`<td class="num">${dash}</td>`;
  const fill = `width:${Math.max(2, d.cov_min)}%;background:${isComplete(d) ? "var(--ord-3)" : "var(--seq-mid)"}`;
  return html`<td class="num"><span class="meter">
    <span class="track"><i class="fill" style=${fill}></i></span><span class="val">${d.cov_min.toFixed(1)}%</span>
  </span></td>`;
}

/* k: sort key and column id; t: header; num: numeric sort (missing values last); get: the
   value sorted on; cell: the <td> */
export const COLS = [
  { k: "id", t: "Sample", num: false, get: (d) => d.id,
    cell: (d) => html`<td class="id">${d.id}</td>` },
  { k: "status", t: "Status", num: false, get: (d) => d.status,
    cell: (d) => html`<td><${Pill} status=${d.status} /></td>` },
  FACET_LABEL ? { k: "facet", t: FACET_LABEL, num: false, get: (d) => d.facet || "",
    cell: (d) => html`<td>${d.facet ? html`<${FacetChip} facet=${d.facet} />` : dash}</td>` } : null,
  { k: "total_reads", t: "Reads", num: true, get: (d) => d.total_reads,
    cell: (d) => numCell(d.total_reads, fmt(d.total_reads)) },
  IS_ILLUMINA ? { k: "qc_passed_reads", t: "QC-passed", num: true, get: (d) => d.qc_passed_reads,
    cell: (d) => numCell(d.qc_passed_reads, fmt(d.qc_passed_reads)) } : null,
  { k: "mapped_reads", t: "Mapped", num: true, get: (d) => d.mapped_reads,
    cell: (d) => numCell(d.mapped_reads, fmt(d.mapped_reads)) },
  { k: "pct_mapped", t: "% mapped", num: true, get: (d) => d.pct_mapped,
    cell: (d) => numCell(d.pct_mapped, pct(d.pct_mapped)) },
  { k: "median_depth", t: "Median depth", num: true, get: (d) => d.median_depth,
    cell: (d) => numCell(d.median_depth, d.median_depth === null ? "—" : fmt1(d.median_depth) + "×") },
  { k: "cov_min", t: `Cov ≥${P.min_depth}×`, num: true, get: (d) => d.cov_min,
    cell: (d) => html`<${CoverageCell} d=${d} />` },
  { k: "n_pct", t: "N %", num: true, get: (d) => d.consensus ? d.consensus.n_pct : null,
    cell: (d) => numCell(d.consensus ? d.consensus.n_pct : null, d.consensus ? pct(d.consensus.n_pct) : "—") },
  HAS_VQC ? { k: "grade", t: "Grade", num: false, get: (d) => gradeOf(d) || "",
    cell: (d) => html`<td><${GradePill} grade=${gradeOf(d)} /></td>` } : null,
  HAS_VQC ? { k: "clade", t: "Clade / lineage", num: false, get: (d) => cladeText(d),
    cell: (d) => html`<td class=${cladeText(d) ? null : "dim"}>${cladeText(d) || "—"}</td>` } : null,
].filter(Boolean);

const FILTERS = [
  { k: "all", t: "All", match: () => true },
  { k: "complete", t: `≥${T}% covered`, match: isComplete },
  HAS_VQC ? { k: "grade", t: "Grade A/B", match: isGoodGrade } : null,
  { k: "flagged", t: "Flagged", match: (d) => d.status !== "ok" },
].filter(Boolean);

/* What the table shows for a view ({ filter, query, sortKey, sortDir }), in order. The
   detail panel's Previous/Next walk the same list. */
export function visibleRows(view) {
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

/* A new column sorts descending if numeric, ascending otherwise; the same one flips. */
const sortedBy = (view, k) => k === view.sortKey
  ? { ...view, sortDir: -view.sortDir }
  : { ...view, sortKey: k, sortDir: COLS.find((c) => c.k === k).num ? -1 : 1 };

export function SampleTable({ view, setView, rows, selected, onSelect }) {
  const activate = (e, key) => {
    if (e.key === "Enter" || e.key === " ") { e.preventDefault(); onSelect(key); }
  };
  return html`
    <section>
      <div class="sec-head">
        <h2>Samples</h2>
        <p id="table-caption">${rows.length} of ${S.length} ${ROWS_NOUN}. Click a row for its detail.</p>
      </div>
      <div class="filters" id="filters">
        ${FILTERS.map((f) => html`
          <button class="fbtn" type="button" data-f=${f.k} aria-pressed=${String(f.k === view.filter)}
            onClick=${() => setView({ ...view, filter: f.k })}>
            ${f.t}<span class="n">${S.filter(f.match).length}</span>
          </button>`)}
        <input type="search" id="q" aria-label="Filter samples" value=${view.query}
          placeholder=${`Filter by sample, ${FACET_LABEL ? FACET_LABEL.toLowerCase() + ", " : ""}virus, clade or status`}
          onInput=${(e) => setView({ ...view, query: e.target.value })} />
      </div>
      <div class="tablewrap">
        <table class="data" id="table">
          <thead><tr id="thead-row">
            ${COLS.map((c) => {
              const on = c.k === view.sortKey;
              return html`<th data-k=${c.k} class=${c.num ? "num" : ""} scope="col"
                  aria-sort=${on ? (view.sortDir === 1 ? "ascending" : "descending") : null}
                  onClick=${() => setView(sortedBy(view, c.k))}>
                ${c.t}<span class="arrow">${on && view.sortDir !== 1 ? "▼" : "▲"}</span>
              </th>`;
            })}
          </tr></thead>
          <tbody id="tbody">
            ${rows.map((d) => html`
              <tr key=${d.key} data-id=${d.key} aria-selected=${String(d.key === selected)} tabindex="0"
                  onClick=${() => onSelect(d.key)} onKeyDown=${(e) => activate(e, d.key)}>
                ${COLS.map((c) => c.cell(d))}
              </tr>`)}
          </tbody>
        </table>
      </div>
      <p class="note" style="margin-top:8px" id="table-foot">
        <b>Cov ≥${P.min_depth}×</b> is the share of reference positions at or above the masking depth;
        positions below it are <code>N</code> in the consensus. A row enters the filtered multi-FASTA at
        ≥${T}%.${HAS_VQC ? html` <b>Grade</b> is viralQC's A–D genome quality.` : null}${
        IS_ILLUMINA ? html` <b>QC-passed</b> counts reads surviving fastp.` : null}
      </p>
    </section>`;
}
