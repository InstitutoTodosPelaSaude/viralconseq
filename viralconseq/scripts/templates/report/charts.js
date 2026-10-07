"use strict";
/* The five charts: hand-drawn SVG, drawn at their host's pixel width and redrawn on
   resize and on a theme change. The template mounts each one on its host element with
   x-init="mount($el, draw…)". Alpine cannot template inside <svg>, so these stay code. */

/* ---------- tooltip: one element for the whole page ---------- */
const tipEl = document.getElementById("tip");
/* rows is a list of [label, value] pairs; both are HTML, so callers esc() data */
function showTip(evt, title, rows) {
  tipEl.innerHTML = `<div class="t">${title}</div>` +
    rows.map(([k, v]) => `<div class="r"><span>${k}</span><b>${v}</b></div>`).join("");
  tipEl.style.opacity = "1";
  tipEl.setAttribute("aria-hidden", "false");
  moveTip(evt);
}
function moveTip(evt) {
  const r = tipEl.getBoundingClientRect();
  let x = evt.clientX + 14, y = evt.clientY + 16;
  if (x + r.width > innerWidth - 8) x = evt.clientX - r.width - 14;
  if (y + r.height > innerHeight - 8) y = evt.clientY - r.height - 16;
  tipEl.style.left = Math.max(8, x) + "px";
  tipEl.style.top = Math.max(8, y) + "px";
}
function hideTip() { tipEl.style.opacity = "0"; tipEl.setAttribute("aria-hidden", "true"); }
/* attach hover + keyboard focus to a mark; both show the same thing */
function hoverable(node, title, rows, onActivate) {
  node.addEventListener("mouseenter", (e) => showTip(e, title, rows));
  node.addEventListener("mousemove", moveTip);
  node.addEventListener("mouseleave", hideTip);
  if (onActivate) {
    node.setAttribute("tabindex", "0");
    node.addEventListener("focus", (e) => {
      const b = node.getBoundingClientRect();
      showTip({ clientX: b.left + b.width / 2, clientY: b.top }, title, rows);
    });
    node.addEventListener("blur", hideTip);
    node.addEventListener("click", onActivate);
    node.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); onActivate(); } });
  }
}

/* ---------- svg helpers ---------- */
const NS = "http://www.w3.org/2000/svg";
function el(name, attrs = {}, parent = null) {
  const n = document.createElementNS(NS, name);
  for (const k in attrs) if (attrs[k] !== null && attrs[k] !== undefined) n.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(n);
  return n;
}
function svgRoot(host, w, h) {
  host.innerHTML = "";
  const s = el("svg", { class: "chart", viewBox: `0 0 ${w} ${h}`, width: w, height: h,
                        preserveAspectRatio: "xMinYMin meet", role: "img" }, host);
  s.style.width = "100%"; s.style.height = "auto";
  return s;
}
function text(parent, x, y, str, cls = "lbl", extra = {}) {
  const t = el("text", Object.assign({ x, y, class: cls }, extra), parent);
  t.textContent = str;
  return t;
}
/* Approximate rendered width — used to decide whether a label fits inside a mark. */
const textWidth = (str, size) => String(str).length * size * 0.56;

/* A label set inside a coloured fill takes white or ink by the fill's luminance,
   so it clears contrast in either theme. */
function inkOn(cssVar) {
  const hex = css(cssVar).replace("#", "");
  if (hex.length !== 6) return "#fff";
  const lum = [0, 2, 4].map((i) => {
    const c = parseInt(hex.slice(i, i + 2), 16) / 255;
    return c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4;
  }).reduce((a, v, i) => a + v * [0.2126, 0.7152, 0.0722][i], 0);
  return (lum + 0.05) / 0.05 > 4.5 ? css("--brand-deep") : "#fff";
}

/* Re-render on resize so type stays one size and nothing is clipped. */
let mounted = [];
function mount(host, render) {
  if (!host) return;
  mounted.push([host, render]);
  render(host, Math.max(240, host.clientWidth || 320));
}
/* Detail charts are replaced with every selection; drop the ones whose host has left
   the document so the resize handler stays bounded. */
const dropDetached = () => { mounted = mounted.filter(([h]) => h.isConnected); };
function redrawAll() {
  dropDetached();
  mounted.forEach(([h, r]) => r(h, Math.max(240, h.clientWidth || 320)));
}
let rAF = null;
addEventListener("resize", () => {
  hideTip(); /* an open tooltip points at where a mark was before the resize */
  if (rAF) cancelAnimationFrame(rAF);
  rAF = requestAnimationFrame(redrawAll);
});

/* ---------- funnel: unique samples surviving each stage ---------- */
function drawFunnel(host, W) {
  const rowH = 58, padT = 4, h = FUNNEL_STAGES.length * rowH + padT;
  const barMax = W - 46;
  const s = svgRoot(host, W, h);
  s.setAttribute("aria-label", "Samples with reads through coverage filtering");
  const top = FUNNEL_STAGES[0].n || 1;
  const nStages = FUNNEL_STAGES.length;
  FUNNEL_STAGES.forEach((st, i) => {
    const y = padT + i * rowH;
    const w = Math.max(3, barMax * st.n / top);
    const ord = `--ord-${nStages > 1 ? Math.round(i * 3 / (nStages - 1)) + 1 : 1}`;
    text(s, 0, y + 12, st.label, "lbl-strong");
    const g = el("g", {}, s);
    el("rect", { x: 0, y: y + 19, width: w, height: 16, rx: 4, fill: `var(${ord})` }, g);
    /* square the baseline end, round the data end */
    el("rect", { x: 0, y: y + 19, width: Math.min(6, w), height: 16, fill: `var(${ord})` }, g);
    const val = `${st.n}`;
    const inside = textWidth(val, 11) + 18 < w;
    text(s, inside ? w - 8 : w + 8, y + 31, val,
         inside ? "in-bar" : "lbl-strong",
         inside ? { "text-anchor": "end", fill: inkOn(ord) } : {});
    text(s, 0, y + 49, st.note + (st.share ? ` · ${(100 * st.n / top).toFixed(0)}% of samples` : ""), "tick");
    hoverable(g, st.label, [["Samples", fmt(st.n)],
      ["Share of run", (100 * st.n / top).toFixed(1) + "%"], ["This stage keeps", st.note]]);
  });
}

/* ---------- ranked completeness: one bar per row, best first ---------- */
function drawRanked(host, W) {
  const padL = 30, padR = 8, padT = 8, plotH = 168, axisH = 26;
  const s = svgRoot(host, W, plotH + axisH + padT);
  s.setAttribute("aria-label", `Genome completeness at ${P.min_depth}x for each row, ranked`);
  const innerW = W - padL - padR;
  const band = innerW / Math.max(1, ranked.length);
  const barW = Math.max(3, Math.min(24, band - 2)); /* 2px surface gap between bars */
  const y = (v) => padT + plotH - (v / 100) * plotH;

  [0, 25, 50, 75, 100].forEach((v) => {
    el("line", { class: "grid-line", x1: padL, x2: W - padR, y1: y(v), y2: y(v) }, s);
    text(s, padL - 6, y(v) + 3.5, v, "tick", { "text-anchor": "end" });
  });
  ranked.forEach((d, i) => {
    const x = padL + i * band + (band - barW) / 2;
    const g = el("g", {}, s);
    el("rect", { x, y: y(d.cov_min), width: barW, height: plotH + padT - y(d.cov_min), rx: 3,
                 fill: facetColor(d.facet) }, g);
    el("rect", { x, y: plotH + padT - 3, width: barW, height: 3, fill: facetColor(d.facet) }, g);
    el("rect", { class: "hit", x: padL + i * band, y: padT, width: band, height: plotH }, g);
    hoverable(g, esc(d.key), [...facetRows(d), [`≥${P.min_depth}× coverage`, pct(d.cov_min)],
                              ["10× coverage", pct(d.cov10)], ["Median depth", fmt1(d.median_depth) + "×"]],
              () => select(d.key));
  });
  const ty = y(T);
  el("line", { class: "rule", x1: padL, x2: W - padR, y1: ty, y2: ty }, s);
  text(s, W - padR, ty - 5, `${T}% threshold`, "rule-lbl", { "text-anchor": "end" });
  el("line", { class: "axis-line", x1: padL, x2: W - padR, y1: padT + plotH, y2: padT + plotH }, s);
  text(s, padL, padT + plotH + 16, `${ranked.length} rows with coverage, best first`, "tick");
  text(s, W - padR, padT + plotH + 16, `${complete.length} clear the threshold`, "tick", { "text-anchor": "end" });
}

/* ---------- depth against completeness: one dot per row, log depth ---------- */
function drawScatter(host, W) {
  const padL = 34, padR = 12, padT = 8, plotH = 168, axisH = 34;
  const s = svgRoot(host, W, plotH + axisH + padT);
  s.setAttribute("aria-label", "Median depth against genome completeness for each row");
  const pts = S.filter((d) => d.median_depth !== null && d.cov_min !== null);
  const lx = (v) => Math.log10(Math.max(1, v));
  const xMax = lx(Math.max(...pts.map((d) => d.median_depth), 10)) * 1.04;
  const x = (v) => padL + (lx(v) / xMax) * (W - padL - padR);
  const y = (v) => padT + plotH - (v / 100) * plotH;

  [0, 25, 50, 75, 100].forEach((v) => {
    el("line", { class: "grid-line", x1: padL, x2: W - padR, y1: y(v), y2: y(v) }, s);
    text(s, padL - 6, y(v) + 3.5, v, "tick", { "text-anchor": "end" });
  });
  [1, 10, 100, 1000, 10000].forEach((v) => {
    if (lx(v) > xMax) return;
    el("line", { class: "grid-line", x1: x(v), x2: x(v), y1: padT, y2: padT + plotH }, s);
    text(s, x(v), padT + plotH + 15, compact(v) + "×", "tick", { "text-anchor": "middle" });
  });
  const ty = y(T);
  el("line", { class: "rule", x1: padL, x2: W - padR, y1: ty, y2: ty }, s);
  text(s, W - padR, ty - 5, `${T}%`, "rule-lbl", { "text-anchor": "end" });

  pts.forEach((d) => {
    const g = el("g", {}, s);
    /* 2px surface ring keeps overlapping dots legible */
    el("circle", { cx: x(d.median_depth), cy: y(d.cov_min), r: 5.5,
                   fill: facetColor(d.facet), stroke: "var(--surface)", "stroke-width": 2 }, g);
    el("circle", { class: "hit", cx: x(d.median_depth), cy: y(d.cov_min), r: 13 }, g);
    hoverable(g, esc(d.key), [...facetRows(d), ["Median depth", fmt1(d.median_depth) + "×"],
                              [`≥${P.min_depth}× coverage`, pct(d.cov_min)],
                              ["Clade", esc((d.viralqc && (d.viralqc.clade || d.viralqc.lineage)) || "—")]],
              () => select(d.key));
  });
  el("line", { class: "axis-line", x1: padL, x2: W - padR, y1: padT + plotH, y2: padT + plotH }, s);
  text(s, padL, padT + plotH + 29, "Median depth (log scale)", "tick");
  text(s, 0, 0, `% covered at ≥${P.min_depth}×`, "tick",
       { transform: `translate(11 ${padT + plotH / 2}) rotate(-90)`, "text-anchor": "middle" });
}

/* ---------- depth along the genome, for row d ---------- */
function drawTrace(host, W, d) {
  const cov = DATA.coverage[d.key];
  const ext = extentOf(d);
  const color = facetColor(d.facet);
  const padL = 34, padR = 8, padT = 8, plotH = 150, axisH = 30;
  const s = svgRoot(host, W, plotH + axisH + padT);
  s.setAttribute("aria-label", `Depth along the genome for ${d.key}`);
  const L = cov.genome_length, innerW = W - padL - padR;
  const x = (p) => padL + (p / L) * innerW;
  const maxD = Math.max(10, ...cov.bins);
  const ly = (v) => Math.log10(v + 1) / Math.log10(maxD + 1);
  const y = (v) => padT + plotH - ly(v) * plotH;
  const contigs = cov.contigs || [];

  if (ext) {
    el("rect", { x: x(0), y: padT, width: Math.max(1, x(ext.start) - x(0)), height: plotH, fill: "var(--surface-sunk)" }, s);
    el("rect", { x: x(ext.end), y: padT, width: Math.max(1, x(L) - x(ext.end)), height: plotH, fill: "var(--surface-sunk)" }, s);
  }
  (cov.mask || []).forEach(([a, b]) => {
    el("rect", { x: x(a - 1), y: padT, width: Math.max(1, x(b) - x(a - 1)), height: plotH, fill: "var(--grid)" }, s);
  });
  [1, 10, 100, 1000, 10000].filter((v) => v <= maxD).forEach((v) => {
    el("line", { class: "grid-line", x1: padL, x2: W - padR, y1: y(v), y2: y(v) }, s);
    text(s, padL - 6, y(v) + 3.5, compact(v), "tick", { "text-anchor": "end" });
  });
  let path = "";
  cov.bins.forEach((v, i) => { path += (i ? "L" : "M") + x(i * cov.bin_bp).toFixed(1) + " " + y(v).toFixed(1); });
  el("path", { d: path + `L${x((cov.bins.length - 1) * cov.bin_bp).toFixed(1)} ${padT + plotH}L${padL} ${padT + plotH}Z`,
               fill: color, "fill-opacity": 0.10 }, s);
  el("path", { d: path, fill: "none", stroke: color, "stroke-width": 2,
               "stroke-linejoin": "round", "stroke-linecap": "round" }, s);
  if (contigs.length > 1) {
    contigs.forEach((c, i) => {
      if (i) el("line", { class: "contig-line", x1: x(c.offset), x2: x(c.offset), y1: padT, y2: padT + plotH }, s);
      const label = short(c.name, Math.max(4, Math.floor((x(c.offset + c.length) - x(c.offset)) / 6.5)));
      text(s, x(c.offset + c.length / 2), padT + 11, label, "tick", { "text-anchor": "middle" });
    });
  }
  const ry = y(P.min_depth);
  el("line", { class: "rule", x1: padL, x2: W - padR, y1: ry, y2: ry }, s);
  text(s, W - padR, ry - 5, `${P.min_depth}×`, "rule-lbl", { "text-anchor": "end" });
  el("line", { class: "axis-line", x1: padL, x2: W - padR, y1: padT + plotH, y2: padT + plotH }, s);
  const step = niceStep(L / 6);
  for (let p = 0; p <= L; p += step)
    text(s, x(p), padT + plotH + 15, compact(p), "tick", { "text-anchor": p ? "middle" : "start" });
  text(s, padL, padT + plotH + 28, `Genome position (bp)${contigs.length === 1 ? ` · ${short(contigs[0].name, 40)}` : contigs.length > 1 ? ` · ${contigs.length} contigs, concatenated` : ""}`, "tick");

  /* crosshair + tooltip across the trace */
  const cross = el("line", { class: "axis-line", y1: padT, y2: padT + plotH, opacity: 0 }, s);
  const dot = el("circle", { r: 4.5, fill: color, stroke: "var(--surface)", "stroke-width": 2, opacity: 0 }, s);
  const hit = el("rect", { class: "hit", x: padL, y: padT, width: innerW, height: plotH }, s);
  hit.addEventListener("mousemove", (e) => {
    const box = s.getBoundingClientRect();
    const px = (e.clientX - box.left) * (W / box.width);
    const i = Math.max(0, Math.min(cov.bins.length - 1, Math.round(((px - padL) / innerW) * L / cov.bin_bp)));
    const pos = i * cov.bin_bp;
    cross.setAttribute("x1", x(pos)); cross.setAttribute("x2", x(pos)); cross.setAttribute("opacity", 1);
    dot.setAttribute("cx", x(pos)); dot.setAttribute("cy", y(cov.bins[i])); dot.setAttribute("opacity", 1);
    const masked = (cov.mask || []).some(([a, b]) => pos + 1 >= a && pos + 1 <= b);
    const contig = contigs.length > 1 ? contigs.find((c) => pos >= c.offset && pos < c.offset + c.length) : null;
    showTip(e, `Position ${fmt(pos + 1)}–${fmt(Math.min(L, pos + cov.bin_bp))}`,
      [...(contig ? [["Contig", esc(contig.name)], ["Contig position", fmt(pos - contig.offset + 1)]] : []),
       ["Mean depth", fmt1(cov.bins[i]) + "×"],
       ["Consensus", masked ? "masked to N" : "called"],
       ...(ext ? [["Primer extent", pos < ext.start || pos > ext.end ? "outside" : "inside"]] : [])]);
  });
  hit.addEventListener("mouseleave", () => { hideTip(); cross.setAttribute("opacity", 0); dot.setAttribute("opacity", 0); });
}

/* ---------- mean depth per amplicon, for row d ---------- */
function drawAmplicons(host, W, d) {
  const depths = DATA.amplicons[d.key];
  const color = facetColor(d.facet);
  const padL = 34, padR = 8, padT = 8, plotH = 120, axisH = 26;
  const s = svgRoot(host, W, plotH + axisH + padT);
  s.setAttribute("aria-label", `Mean depth per amplicon for ${d.key}`);
  const innerW = W - padL - padR, band = innerW / depths.length;
  const barW = Math.max(3, Math.min(24, band - 2));
  const dMax = Math.max(...depths, P.min_depth) * 1.06;
  const y = (v) => padT + plotH - (v / dMax) * plotH;
  const step = niceStep(dMax / 4);
  for (let v = 0; v <= dMax; v += step) {
    el("line", { class: "grid-line", x1: padL, x2: W - padR, y1: y(v), y2: y(v) }, s);
    text(s, padL - 6, y(v) + 3.5, compact(v), "tick", { "text-anchor": "end" });
  }
  depths.forEach((v, i) => {
    const xx = padL + i * band + (band - barW) / 2;
    const g = el("g", {}, s);
    const low = v < P.min_depth;
    el("rect", { x: xx, y: y(v), width: barW, height: Math.max(2, padT + plotH - y(v)), rx: 3,
                 fill: low ? "var(--accent)" : color }, g);
    /* a true zero would otherwise render as nothing, hiding which amplicon failed */
    if (v < 1) el("rect", { x: xx, y: padT + plotH - 3, width: barW, height: 3, fill: "var(--accent)" }, g);
    el("rect", { class: "hit", x: padL + i * band, y: padT, width: band, height: plotH }, g);
    hoverable(g, `Amplicon ${i + 1}`, [["Mean depth", fmt1(v) + "×"],
      ["Against threshold", low ? `below ${P.min_depth}×` : `at or above ${P.min_depth}×`]]);
    if ((i + 1) % 5 === 0 || i === 0) text(s, padL + i * band + band / 2, padT + plotH + 14, i + 1, "tick", { "text-anchor": "middle" });
  });
  const ry = y(P.min_depth);
  el("line", { class: "rule", x1: padL, x2: W - padR, y1: ry, y2: ry }, s);
  el("line", { class: "axis-line", x1: padL, x2: W - padR, y1: padT + plotH, y2: padT + plotH }, s);
  const below = depths.filter((v) => v < P.min_depth).length;
  text(s, padL, padT + plotH + 26, `${below} of ${depths.length} amplicons below ${P.min_depth}×`, "tick");
}
