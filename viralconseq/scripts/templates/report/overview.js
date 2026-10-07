/* Run overview: the funnel, ranked completeness and depth against completeness. */
import {
  html, FK, HAS_VQC, IDS, IS_ILLUMINA, P, S, T,
  anyRow, compact, complete, facetColor, facetCount, facetRows, fmt, fmt1, goodGrade, pct, ranked, totalReads,
} from "report/core";
import { Chart, ChartSvg, FacetChip, Hover, inkOn, percentGrid, svgText, textWidth } from "report/ui";

/* ---------- funnel: unique samples surviving each stage ---------- */
const STAGES = [
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

function drawFunnel(W) {
  const rowH = 58, padT = 4, h = STAGES.length * rowH + padT;
  const barMax = W - 46;
  const top = STAGES[0].n || 1;
  return html`<${ChartSvg} w=${W} h=${h} label="Samples with reads through coverage filtering">
    ${STAGES.map((st, i) => {
      const y = padT + i * rowH;
      const w = Math.max(3, barMax * st.n / top);
      const ord = `--ord-${STAGES.length > 1 ? Math.round(i * 3 / (STAGES.length - 1)) + 1 : 1}`;
      const val = `${st.n}`;
      const inside = textWidth(val, 11) + 18 < w;
      return html`
        ${svgText(0, y + 12, st.label, "lbl-strong")}
        <${Hover} title=${st.label} rows=${[["Samples", fmt(st.n)],
            ["Share of run", (100 * st.n / top).toFixed(1) + "%"], ["This stage keeps", st.note]]}>
          <rect x="0" y=${y + 19} width=${w} height="16" rx="4" fill=${`var(${ord})`} />
          <!-- square the baseline end, round the data end -->
          <rect x="0" y=${y + 19} width=${Math.min(6, w)} height="16" fill=${`var(${ord})`} />
        <//>
        ${svgText(inside ? w - 8 : w + 8, y + 31, val, inside ? "in-bar" : "lbl-strong",
                  inside ? { "text-anchor": "end", fill: inkOn(ord) } : {})}
        ${svgText(0, y + 49, st.note + (st.share ? ` · ${(100 * st.n / top).toFixed(0)}% of samples` : ""), "tick")}`;
    })}
  <//>`;
}

/* ---------- ranked completeness: one bar per row, best first ---------- */
function drawRanked(W, onSelect) {
  const padL = 30, padR = 8, padT = 8, plotH = 168, axisH = 26;
  const innerW = W - padL - padR;
  const band = innerW / Math.max(1, ranked.length);
  const barW = Math.max(3, Math.min(24, band - 2)); /* 2px surface gap between bars */
  const y = (v) => padT + plotH - (v / 100) * plotH;
  const ty = y(T);
  return html`<${ChartSvg} w=${W} h=${plotH + axisH + padT}
      label=${`Genome completeness at ${P.min_depth}x for each row, ranked`}>
    ${percentGrid(y, padL, W - padR)}
    ${ranked.map((d, i) => {
      const x = padL + i * band + (band - barW) / 2;
      return html`<${Hover} title=${d.key} onActivate=${() => onSelect(d.key)}
          rows=${[...facetRows(d), [`≥${P.min_depth}× coverage`, pct(d.cov_min)],
                  ["10× coverage", pct(d.cov10)], ["Median depth", fmt1(d.median_depth) + "×"]]}>
        <rect x=${x} y=${y(d.cov_min)} width=${barW} height=${plotH + padT - y(d.cov_min)} rx="3" fill=${facetColor(d.facet)} />
        <rect x=${x} y=${plotH + padT - 3} width=${barW} height="3" fill=${facetColor(d.facet)} />
        <rect class="hit" x=${padL + i * band} y=${padT} width=${band} height=${plotH} />
      <//>`;
    })}
    <line class="rule" x1=${padL} x2=${W - padR} y1=${ty} y2=${ty} />
    ${svgText(W - padR, ty - 5, `${T}% threshold`, "rule-lbl", { "text-anchor": "end" })}
    <line class="axis-line" x1=${padL} x2=${W - padR} y1=${padT + plotH} y2=${padT + plotH} />
    ${svgText(padL, padT + plotH + 16, `${ranked.length} rows with coverage, best first`, "tick")}
    ${svgText(W - padR, padT + plotH + 16, `${complete.length} clear the threshold`, "tick", { "text-anchor": "end" })}
  <//>`;
}

/* ---------- depth against completeness: one dot per row, log depth ---------- */
function drawScatter(W, onSelect) {
  const padL = 34, padR = 12, padT = 8, plotH = 168, axisH = 34;
  const pts = S.filter((d) => d.median_depth !== null && d.cov_min !== null);
  const lx = (v) => Math.log10(Math.max(1, v));
  const xMax = lx(Math.max(...pts.map((d) => d.median_depth), 10)) * 1.04;
  const x = (v) => padL + (lx(v) / xMax) * (W - padL - padR);
  const y = (v) => padT + plotH - (v / 100) * plotH;
  const ty = y(T);
  return html`<${ChartSvg} w=${W} h=${plotH + axisH + padT} label="Median depth against genome completeness for each row">
    ${percentGrid(y, padL, W - padR)}
    ${[1, 10, 100, 1000, 10000].filter((v) => lx(v) <= xMax).map((v) => html`
      <line class="grid-line" x1=${x(v)} x2=${x(v)} y1=${padT} y2=${padT + plotH} />
      ${svgText(x(v), padT + plotH + 15, compact(v) + "×", "tick", { "text-anchor": "middle" })}`)}
    <line class="rule" x1=${padL} x2=${W - padR} y1=${ty} y2=${ty} />
    ${svgText(W - padR, ty - 5, `${T}%`, "rule-lbl", { "text-anchor": "end" })}
    ${pts.map((d) => html`<${Hover} title=${d.key} onActivate=${() => onSelect(d.key)}
        rows=${[...facetRows(d), ["Median depth", fmt1(d.median_depth) + "×"],
                [`≥${P.min_depth}× coverage`, pct(d.cov_min)],
                ["Clade", (d.viralqc && (d.viralqc.clade || d.viralqc.lineage)) || "—"]]}>
      <!-- 2px surface ring keeps overlapping dots legible -->
      <circle cx=${x(d.median_depth)} cy=${y(d.cov_min)} r="5.5"
        fill=${facetColor(d.facet)} stroke="var(--surface)" stroke-width="2" />
      <circle class="hit" cx=${x(d.median_depth)} cy=${y(d.cov_min)} r="13" />
    <//>`)}
    <line class="axis-line" x1=${padL} x2=${W - padR} y1=${padT + plotH} y2=${padT + plotH} />
    ${svgText(padL, padT + plotH + 29, "Median depth (log scale)", "tick")}
    ${svgText(0, 0, `% covered at ≥${P.min_depth}×`, "tick",
              { transform: `translate(11 ${padT + plotH / 2}) rotate(-90)`, "text-anchor": "middle" })}
  <//>`;
}

/* one chip per facet with its row count */
const FacetLegend = () => facetCount.map(([f, n]) =>
  html`<${FacetChip} facet=${f}> <span style="color:var(--muted);font-weight:500">${n}</span><//>`);

export function Overview({ onSelect }) {
  return html`
    <section>
      <div class="sec-head">
        <h2>Run overview</h2>
        <p>What the run yielded, where samples were lost, and how depth converts into a genome.</p>
      </div>
      <div class="grid-3">
        <div class="card">
          <h3>From sample to genome</h3>
          <p class="note" id="funnel-note">${FK === "segment"
            ? "Each stage counts unique samples. A sample with several segments contributes once here and one row per segment in the table."
            : "Each stage counts unique samples."}</p>
          <${Chart} id="funnel" draw=${drawFunnel} />
        </div>
        <div class="card">
          <h3>Genome completeness, ranked</h3>
          <p class="note" id="ranked-note">Percent of reference positions with depth ≥${P.min_depth}×. Masking and variant filters also affect the consensus.</p>
          <${Chart} id="ranked" draw=${(W) => drawRanked(W, onSelect)} />
          <div class="legend" id="ranked-legend"><${FacetLegend} /></div>
        </div>
        <div class="card">
          <h3>Depth against completeness</h3>
          <p class="note">Median depth over every reference position, zeros included.
            The relationship shown is descriptive and does not establish the cause of low coverage.</p>
          <${Chart} id="scatter" draw=${(W) => drawScatter(W, onSelect)} />
          <div class="legend" id="scatter-legend"><${FacetLegend} /></div>
        </div>
      </div>
    </section>`;
}
