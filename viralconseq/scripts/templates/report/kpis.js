/* The headline tiles under the masthead. */
import {
  html, FACET_PLURAL, FK, HAS_VQC, IDS, P, ROWS_NOUN, S, T,
  complete, facetCount, fmt, fmt1, goodGrade, gradeOf, medDepth, medMapped, totalReads,
} from "report/core";

const TILES = [
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

export function Kpis() {
  return html`
    <section id="kpi-section">
      <div class="kpis" id="kpis" data-n=${TILES.length}>
        ${TILES.map((k) => html`
          <div class="kpi">
            <div class="label">${k.label}</div>
            <div class="value">${k.value}${k.unit && html`<span class="unit">${k.unit}</span>`}</div>
            <div class="foot">${k.foot}</div>
          </div>`)}
      </div>
    </section>`;
}
