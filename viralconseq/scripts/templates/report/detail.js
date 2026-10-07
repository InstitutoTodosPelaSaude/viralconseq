/* The selected row: headline figures, read QC or variant calling, the depth trace,
   amplicon depths and the viralQC verdict. */
import { useState } from "preact/hooks";
import {
  html, DATA, IS_ILLUMINA, P, S,
  compact, facetColor, fmt, fmt1, frac, gradeOf, niceStep, pct, pctOf, short, statusMeta,
} from "report/core";
import {
  Chart, ChartSvg, FacetChip, GradePill, Hover, Pill, QcPill, StatStrip, hideTip, showTip, svgText,
} from "report/ui";

const small = (text) => html`<small>${text}</small>`;

/* ---------- headline figures ---------- */
function headlineStats(d) {
  const items = [{ k: "Reads", v: fmt(d.total_reads) }];
  if (IS_ILLUMINA) items.push({ k: "QC-passed", v: fmt(d.qc_passed_reads),
    n: d.total_reads ? `${pctOf(d.qc_passed_reads, d.total_reads)} of reads` : "" });
  items.push({ k: "Mapped", v: fmt(d.mapped_reads), n: `${pct(d.pct_mapped)} of reads` });
  items.push({ k: "Depth", v: html`${fmt1(d.median_depth)}${small("× median")}`,
    n: `mean ${fmt1(d.mean_depth)}×, zeros included` });
  items.push({ k: `Covered at ≥${P.min_depth}×`, v: pct(d.cov_min), n: `10×: ${pct(d.cov10)} · 100×: ${pct(d.cov100)}` });
  if (d.consensus) items.push({ k: "Consensus", v: html`${fmt(d.consensus.length)}${small(" bp")}`,
    n: `${fmt(d.consensus.n)} N (${pct(d.consensus.n_pct)})` });
  if (d.variants) items.push({ k: "Variant records", v: fmt(d.variants.records), n: "lines in the final VCF" });
  if (d.isnv_count !== null) items.push({ k: "iSNVs", v: fmt(d.isnv_count),
    n: P.af_threshold !== null ? `AF ≥ ${P.af_threshold}` : "" });
  if (d.clair3_model) items.push({ k: "Clair3 model", v: small(d.clair3_model) });
  return items;
}

/* Why a row is not plain "ok", in words. */
function StatusCallout({ d }) {
  const g = d.gate || {};
  switch (d.status) {
    case "ok": return null;
    case "no_mapped_reads":
      return html`<p class="callout">Assembly was not attempted: <b>${fmt(g.mapped_reads ?? d.mapped_reads)}</b> reads mapped to
        the reference${g.minimum_mapped_reads != null ? `, below the ${fmt(g.minimum_mapped_reads)}-read floor` : ""}.
        This gate alone does not establish absence of the target or explain the low mapped-read count.</p>`;
    case "empty_consensus":
      return html`<p class="callout flag">The consensus for this row contains no called bases; every position fell
        below ${P.min_depth}× or the alignment produced nothing to call.</p>`;
    case "viralqc_failed": case "viralqc_partial": case "viralqc_skipped": case "viralqc_missing":
      return html`<p class="callout">viralQC did not grade this genome (<b>${statusMeta(d.status).text}</b>).
        Coverage and depth below are unaffected; clade, lineage and grade are unavailable.</p>`;
    case "missing_stats":
      return html`<p class="callout flag">Assembly statistics are missing for this row; the figures shown come only from
        the artefacts that were present.</p>`;
    default:
      return html`<p class="callout">Status <b>${d.status}</b>.</p>`;
  }
}

/* ---------- left column: read QC (Illumina) or variant calling (Nanopore) ---------- */
function ReadQc({ d }) {
  const f = d.fastp;
  return html`<div>
    <h3 class="sub-h3">Read QC (fastp)</h3>
    <p class="note">Paired reads before and after adapter and quality trimming; the mapped and
      depth figures above are computed from the reads that passed.</p>
    ${f ? html`<${StatStrip} items=${[
      { k: "Reads in", v: fmt(f.total_reads_before) },
      { k: "Reads out", v: fmt(f.total_reads_after), n: `${pctOf(f.total_reads_after, f.total_reads_before)} kept` },
      { k: "Q30 rate", v: frac(f.q30_after), n: f.q30_before === null || f.q30_before === undefined ? "" : `before: ${frac(f.q30_before)}` },
      { k: "Duplication", v: frac(f.duplication) },
      { k: "Mean length", v: html`${fmt(f.read1_mean_length)}${small(` / ${fmt(f.read2_mean_length)} bp`)}`, n: "R1 / R2 after trimming" },
      { k: "Low quality", v: fmt(f.low_quality), n: "reads removed" },
      { k: "Too short", v: fmt(f.too_short), n: "reads removed" },
      { k: "Adapter trimmed", v: fmt(f.adapter_trimmed), n: "reads" },
    ]} />` : html`<p class="note">No fastp report for this sample.</p>`}
  </div>`;
}

function VariantCalling({ d }) {
  return html`<div>
    <h3 class="sub-h3">Variant calling</h3>
    <p class="note">Clair3 calls on reads passing the length and mapping-quality floors; the consensus
      applies variants whose allele fraction, ALT over REF plus ALT reads, is at or above the
      threshold below, and masks positions below ${P.min_depth}×.</p>
    <${StatStrip} items=${[
      { k: "Clair3 model", v: small(d.clair3_model || P.clair3_model || "—") },
      { k: "MAPQ floor", v: fmt(P.minimum_map_quality) },
      { k: "Allele fraction", v: P.af_threshold === null ? "—" : String(P.af_threshold) },
      { k: "Min read length", v: P.minimum_length === null ? "—" : html`${fmt(P.minimum_length)}${small(" bp")}` },
      { k: "Variant records", v: d.variants ? fmt(d.variants.records) : "—" },
    ]} />
  </div>`;
}

/* ---------- right column: depth along the genome ---------- */
function DepthTrace({ d, cov, ext, color }) {
  const [hover, setHover] = useState(null); /* bin index under the pointer */
  const draw = (W) => {
    const padL = 34, padR = 8, padT = 8, plotH = 150, axisH = 30;
    const L = cov.genome_length, innerW = W - padL - padR;
    const x = (p) => padL + (p / L) * innerW;
    const maxD = Math.max(10, ...cov.bins);
    const ly = (v) => Math.log10(v + 1) / Math.log10(maxD + 1);
    const y = (v) => padT + plotH - ly(v) * plotH;
    const contigs = cov.contigs || [];
    const mask = cov.mask || [];

    let path = "";
    cov.bins.forEach((v, i) => { path += (i ? "L" : "M") + x(i * cov.bin_bp).toFixed(1) + " " + y(v).toFixed(1); });
    const ticks = [];
    const step = niceStep(L / 6);
    for (let p = 0; p <= L; p += step) ticks.push(p);
    const ry = y(P.min_depth);
    const hx = hover === null ? null : x(hover * cov.bin_bp);

    const onMove = (e) => {
      const box = e.currentTarget.ownerSVGElement.getBoundingClientRect();
      const px = (e.clientX - box.left) * (W / box.width);
      const i = Math.max(0, Math.min(cov.bins.length - 1, Math.round(((px - padL) / innerW) * L / cov.bin_bp)));
      const pos = i * cov.bin_bp;
      const masked = mask.some(([a, b]) => pos + 1 >= a && pos + 1 <= b);
      const contig = contigs.length > 1 ? contigs.find((c) => pos >= c.offset && pos < c.offset + c.length) : null;
      setHover(i);
      showTip(e, `Position ${fmt(pos + 1)}–${fmt(Math.min(L, pos + cov.bin_bp))}`,
        [...(contig ? [["Contig", contig.name], ["Contig position", fmt(pos - contig.offset + 1)]] : []),
         ["Mean depth", fmt1(cov.bins[i]) + "×"],
         ["Consensus", masked ? "masked to N" : "called"],
         ...(ext ? [["Primer extent", pos < ext.start || pos > ext.end ? "outside" : "inside"]] : [])]);
    };
    const onLeave = () => { hideTip(); setHover(null); };

    return html`<${ChartSvg} w=${W} h=${plotH + axisH + padT} label=${`Depth along the genome for ${d.key}`}>
      ${ext && html`
        <rect x=${x(0)} y=${padT} width=${Math.max(1, x(ext.start) - x(0))} height=${plotH} fill="var(--surface-sunk)" />
        <rect x=${x(ext.end)} y=${padT} width=${Math.max(1, x(L) - x(ext.end))} height=${plotH} fill="var(--surface-sunk)" />`}
      ${mask.map(([a, b]) => html`
        <rect x=${x(a - 1)} y=${padT} width=${Math.max(1, x(b) - x(a - 1))} height=${plotH} fill="var(--grid)" />`)}
      ${[1, 10, 100, 1000, 10000].filter((v) => v <= maxD).map((v) => html`
        <line class="grid-line" x1=${padL} x2=${W - padR} y1=${y(v)} y2=${y(v)} />
        ${svgText(padL - 6, y(v) + 3.5, compact(v), "tick", { "text-anchor": "end" })}`)}
      <path d=${path + `L${x((cov.bins.length - 1) * cov.bin_bp).toFixed(1)} ${padT + plotH}L${padL} ${padT + plotH}Z`}
        fill=${color} fill-opacity="0.1" />
      <path d=${path} fill="none" stroke=${color} stroke-width="2" stroke-linejoin="round" stroke-linecap="round" />
      ${contigs.length > 1 && contigs.map((c, i) => html`
        ${i ? html`<line class="contig-line" x1=${x(c.offset)} x2=${x(c.offset)} y1=${padT} y2=${padT + plotH} />` : null}
        ${svgText(x(c.offset + c.length / 2), padT + 11,
                  short(c.name, Math.max(4, Math.floor((x(c.offset + c.length) - x(c.offset)) / 6.5))),
                  "tick", { "text-anchor": "middle" })}`)}
      <line class="rule" x1=${padL} x2=${W - padR} y1=${ry} y2=${ry} />
      ${svgText(W - padR, ry - 5, `${P.min_depth}×`, "rule-lbl", { "text-anchor": "end" })}
      <line class="axis-line" x1=${padL} x2=${W - padR} y1=${padT + plotH} y2=${padT + plotH} />
      ${ticks.map((p) => svgText(x(p), padT + plotH + 15, compact(p), "tick", { "text-anchor": p ? "middle" : "start" }))}
      ${svgText(padL, padT + plotH + 28, `Genome position (bp)${contigs.length === 1 ? ` · ${short(contigs[0].name, 40)}`
        : contigs.length > 1 ? ` · ${contigs.length} contigs, concatenated` : ""}`, "tick")}
      <!-- crosshair + tooltip across the trace -->
      <line class="axis-line" x1=${hx} x2=${hx} y1=${padT} y2=${padT + plotH} opacity=${hover === null ? 0 : 1} />
      <circle cx=${hx} cy=${hover === null ? null : y(cov.bins[hover])} r="4.5" fill=${color}
        stroke="var(--surface)" stroke-width="2" opacity=${hover === null ? 0 : 1} />
      <rect class="hit" x=${padL} y=${padT} width=${innerW} height=${plotH} onMouseMove=${onMove} onMouseLeave=${onLeave} />
    <//>`;
  };
  return html`<${Chart} id="covtrace" draw=${draw} />`;
}

function DepthColumn({ d, color }) {
  const cov = DATA.coverage[d.key];
  const ext = DATA.run.extents[d.facet || "*"] || null;
  if (!cov) return html`<div>
    <h3 class="sub-h3">Depth along the genome</h3>
    <p class="note">No per-position depth table for this row.</p>
  </div>`;
  return html`<div>
    <h3 class="sub-h3">Depth along the genome</h3>
    <p class="note">Mean depth in ${fmt(cov.bin_bp)} bp windows from <code>table_cov_basewise.txt</code>.
      Grey bands are ranges below ${P.min_depth}× and masked to <code>N</code>${
      cov.contigs && cov.contigs.length > 1 ? `; dashed lines separate the ${cov.contigs.length} reference contigs` : ""}${
      ext ? `; the pale bands at each end lie outside the primer extent (${fmt(ext.outside_bp)} bp) and describe scheme geometry, not a measured coverage limit` : ""}.</p>
    <${DepthTrace} key=${d.key} d=${d} cov=${cov} ext=${ext} color=${color} />
  </div>`;
}

/* ---------- amplicon depths ---------- */
function drawAmplicons(W, d, depths, color) {
  const padL = 34, padR = 8, padT = 8, plotH = 120, axisH = 26;
  const innerW = W - padL - padR, band = innerW / depths.length;
  const barW = Math.max(3, Math.min(24, band - 2));
  const dMax = Math.max(...depths, P.min_depth) * 1.06;
  const y = (v) => padT + plotH - (v / dMax) * plotH;
  const grid = [];
  const step = niceStep(dMax / 4);
  for (let v = 0; v <= dMax; v += step) grid.push(v);
  const ry = y(P.min_depth);
  const below = depths.filter((v) => v < P.min_depth).length;
  return html`<${ChartSvg} w=${W} h=${plotH + axisH + padT} label=${`Mean depth per amplicon for ${d.key}`}>
    ${grid.map((v) => html`
      <line class="grid-line" x1=${padL} x2=${W - padR} y1=${y(v)} y2=${y(v)} />
      ${svgText(padL - 6, y(v) + 3.5, compact(v), "tick", { "text-anchor": "end" })}`)}
    ${depths.map((v, i) => {
      const xx = padL + i * band + (band - barW) / 2;
      const low = v < P.min_depth;
      return html`
        <${Hover} title=${`Amplicon ${i + 1}`} rows=${[["Mean depth", fmt1(v) + "×"],
            ["Against threshold", low ? `below ${P.min_depth}×` : `at or above ${P.min_depth}×`]]}>
          <rect x=${xx} y=${y(v)} width=${barW} height=${Math.max(2, padT + plotH - y(v))} rx="3"
            fill=${low ? "var(--accent)" : color} />
          <!-- a true zero would otherwise render as nothing, hiding which amplicon failed -->
          ${v < 1 ? html`<rect x=${xx} y=${padT + plotH - 3} width=${barW} height="3" fill="var(--accent)" />` : null}
          <rect class="hit" x=${padL + i * band} y=${padT} width=${band} height=${plotH} />
        <//>
        ${(i + 1) % 5 === 0 || i === 0
          ? svgText(padL + i * band + band / 2, padT + plotH + 14, i + 1, "tick", { "text-anchor": "middle" }) : null}`;
    })}
    <line class="rule" x1=${padL} x2=${W - padR} y1=${ry} y2=${ry} />
    <line class="axis-line" x1=${padL} x2=${W - padR} y1=${padT + plotH} y2=${padT + plotH} />
    ${svgText(padL, padT + plotH + 26, `${below} of ${depths.length} amplicons below ${P.min_depth}×`, "tick")}
  <//>`;
}

function Amplicons({ d, depths, color }) {
  return html`<div style="margin-top:18px">
    <h3 class="sub-h3">Amplicon coverage</h3>
    <p class="note">Mean depth over each amplicon of <code>${P.primer_scheme || "the primer scheme"}</code>,
      from the same per-position table. A mean above ${P.min_depth}× does not guarantee every position
      in the amplicon met the masking threshold; use the trace for that.</p>
    <${Chart} key=${d.key} id="sampleamp" draw=${(W) => drawAmplicons(W, d, depths, color)} />
    <div class="legend">
      <span class="chip"><i class="sw" style=${"background:" + color}></i>at or above ${P.min_depth}×</span>
      <span class="chip"><i class="sw" style="background:var(--accent)"></i>mean below ${P.min_depth}×</span>
    </div>
  </div>`;
}

/* ---------- viralQC ---------- */
const scoreNote = (vq) => vq.score !== null ? ` · ${fmt1(vq.score)}/24` : "";

function Viralqc({ vq, grade }) {
  const checks = Object.entries(vq.checks || {}).filter(([, v]) => v);
  return html`<div style="margin-top:18px">
    <h3 class="sub-h3">viralQC</h3>
    <p class="note">Virus <b>${vq.virus || "—"}</b>${vq.species ? ` (${vq.species})` : ""} ·
      clade <b>${vq.clade || "—"}</b> · lineage <b>${vq.lineage || "—"}</b>.</p>
    <div class="legend" style="margin-top:8px">
      <span class="chip">Grade <${GradePill} grade=${grade} extra=${scoreNote(vq)} /></span>
      <span class="chip">Overall <${QcPill} status=${vq.overall} /></span>
      ${checks.map(([k, v]) => html`<span class="chip">${k.replace(/([A-Z])/g, " $1").toLowerCase()} <${QcPill} status=${v} /></span>`)}
    </div>
    <p class="note" style="margin-top:8px">Alignment coverage ${vq.coverage !== null ? (100 * vq.coverage).toFixed(1) + "%" : "—"} · ${fmt(vq.substitutions)} substitutions
      · ${fmt(vq.deletions)} deletions · ${fmt(vq.missing)} missing sites · dataset <b>${vq.dataset || "—"}</b>${
      vq.dataset_version ? html` <b>${vq.dataset_version}</b>` : null}.${
      vq.input_status ? html` Input sequence status: <b>${vq.input_status}</b>.` : null}</p>
    <p class="note" style="margin-top:8px">Grades and QC statuses are reported from the selected viralQC
      dataset. Inspect the individual checks alongside coverage and missing sites; coverage alone does
      not replace the other quality checks.</p>
  </div>`;
}

/* ---------- the panel ---------- */

/* rows: the table's visible rows, which Previous/Next walk; onSelect(key, scroll) */
export function SampleDetail({ selected, rows, onSelect }) {
  const d = S.find((x) => x.key === selected);
  if (!d) return html`<p class="note">Select a sample from the table above.</p>`;
  const idx = rows.findIndex((r) => r.key === d.key);
  const vq = d.viralqc;
  const grade = gradeOf(d);
  const color = facetColor(d.facet);
  const amplicons = DATA.amplicons[d.key];
  return html`
    <div class="detail-head">
      <h3>${d.id}</h3>
      <${Pill} status=${d.status} />
      <${FacetChip} facet=${d.facet} />
      ${grade ? html`<${GradePill} grade=${grade} extra=${scoreNote(vq)} />` : null}
      <div class="detail-nav">
        <button type="button" id="prev" disabled=${idx <= 0}
          onClick=${() => onSelect(rows[idx - 1].key, false)}>← Previous</button>
        <button type="button" id="next" disabled=${idx < 0 || idx >= rows.length - 1}
          onClick=${() => onSelect(rows[idx + 1].key, false)}>Next →</button>
      </div>
    </div>
    <${StatStrip} items=${headlineStats(d)} />
    <${StatusCallout} d=${d} />
    <div class="grid-2">
      ${IS_ILLUMINA ? html`<${ReadQc} d=${d} />` : html`<${VariantCalling} d=${d} />`}
      <${DepthColumn} d=${d} color=${color} />
    </div>
    ${amplicons ? html`<${Amplicons} d=${d} depths=${amplicons} color=${color} />` : null}
    ${vq ? html`<${Viralqc} vq=${vq} grade=${grade} />` : null}`;
}
