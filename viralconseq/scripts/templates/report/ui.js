/* Small pieces the sections share: pills and chips, stat strips, the tooltip, chart hosts. */
import { render } from "preact";
import { useLayoutEffect, useRef, useState } from "preact/hooks";
import { html, css, facetColor, statusMeta } from "report/core";

/* ---------- pills, chips, stat strips ---------- */
export function Pill({ status }) {
  const m = statusMeta(status);
  return html`<span class=${"pill " + m.cls}><i class="dot"></i>${m.text}</span>`;
}

export function GradePill({ grade, extra = "" }) {
  if (!grade) return html`<span class="dim">—</span>`;
  return html`<span class=${`pill grade grade-${String(grade).toLowerCase()}`}><i class="dot"></i>${grade}${extra}</span>`;
}

/* viralQC's overall and per-check verdicts */
export function QcPill({ status }) {
  const s = (status || "").toLowerCase();
  const cls = s === "good" || s === "pass" ? "good" : s === "mediocre" || s === "warn" ? "warning"
            : s === "bad" || s === "fail" ? "critical" : "neutral";
  return html`<span class=${"pill " + cls}><i class="dot"></i>${status || "—"}</span>`;
}

export function FacetChip({ facet, children }) {
  if (!facet) return null;
  return html`<span class="chip"><i class="sw" style=${"background:" + facetColor(facet)}></i>${facet}${children}</span>`;
}

/* A row of figures: k is the label, v the value (text or markup), n an optional note. */
export function StatStrip({ items }) {
  return html`<div class="statstrip">
    ${items.map((i) => html`<div class="s">
      <div class="k">${i.k}</div><div class="v">${i.v}</div>${i.n ? html`<div class="n">${i.n}</div>` : null}
    </div>`)}
  </div>`;
}

/* ---------- tooltip: one element for the whole page, outside the app root ---------- */
const tipEl = document.getElementById("tip");

/* rows is a list of [label, value] pairs */
export function showTip(evt, title, rows) {
  render(html`<div class="t">${title}</div>${rows.map(([k, v]) => html`<div class="r"><span>${k}</span><b>${v}</b></div>`)}`,
         tipEl);
  tipEl.style.opacity = "1";
  tipEl.setAttribute("aria-hidden", "false");
  moveTip(evt);
}

export function moveTip(evt) {
  const r = tipEl.getBoundingClientRect();
  let x = evt.clientX + 14, y = evt.clientY + 16;
  if (x + r.width > innerWidth - 8) x = evt.clientX - r.width - 14;
  if (y + r.height > innerHeight - 8) y = evt.clientY - r.height - 16;
  tipEl.style.left = Math.max(8, x) + "px";
  tipEl.style.top = Math.max(8, y) + "px";
}

export function hideTip() {
  tipEl.style.opacity = "0";
  tipEl.setAttribute("aria-hidden", "true");
}

/* A chart mark with a tooltip on hover and on keyboard focus; onActivate makes it
   focusable and clickable (Enter and Space too). */
export function Hover({ title, rows, onActivate, children }) {
  const focus = (e) => {
    const b = e.currentTarget.getBoundingClientRect();
    showTip({ clientX: b.left + b.width / 2, clientY: b.top }, title, rows);
  };
  const key = (e) => {
    if (e.key === "Enter" || e.key === " ") { e.preventDefault(); onActivate(); }
  };
  return html`<g onMouseEnter=${(e) => showTip(e, title, rows)} onMouseMove=${moveTip} onMouseLeave=${hideTip}
      tabindex=${onActivate ? 0 : null} onFocus=${onActivate && focus} onBlur=${onActivate && hideTip}
      onClick=${onActivate} onKeyDown=${onActivate && key}>${children}</g>`;
}

/* ---------- charts ---------- */

/* Approximate rendered width — used to decide whether a label fits inside a mark. */
export const textWidth = (str, size) => String(str).length * size * 0.56;

/* A label set inside a coloured fill takes white or ink by the fill's luminance,
   so it clears contrast in either theme. */
export function inkOn(cssVar) {
  const hex = css(cssVar).replace("#", "");
  if (hex.length !== 6) return "#fff";
  const lum = [0, 2, 4].map((i) => {
    const c = parseInt(hex.slice(i, i + 2), 16) / 255;
    return c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4;
  }).reduce((a, v, i) => a + v * [0.2126, 0.7152, 0.0722][i], 0);
  return (lum + 0.05) / 0.05 > 4.5 ? css("--brand-deep") : "#fff";
}

/* A chart's host: draws once it knows its width, and again on resize, so type stays
   one size and nothing is clipped. draw(W) returns the chart for a W-pixel-wide host. */
export function Chart({ id, draw }) {
  const ref = useRef(null);
  const [width, setWidth] = useState(null);
  useLayoutEffect(() => {
    let frame = null;
    const measure = () => setWidth(Math.max(240, ref.current.clientWidth || 320));
    /* an open tooltip points at where a mark was before the resize */
    const onResize = () => { hideTip(); if (frame) cancelAnimationFrame(frame); frame = requestAnimationFrame(measure); };
    measure();
    addEventListener("resize", onResize);
    return () => { removeEventListener("resize", onResize); if (frame) cancelAnimationFrame(frame); };
  }, []);
  return html`<div id=${id} ref=${ref}>${width === null ? null : draw(width)}</div>`;
}

/* The SVG every chart draws into: W×H user units, scaled to the host's width. */
export function ChartSvg({ w, h, label, children }) {
  return html`<svg class="chart" viewBox=${`0 0 ${w} ${h}`} width=${w} height=${h}
      preserveAspectRatio="xMinYMin meet" role="img" aria-label=${label}
      style="width: 100%; height: auto;">${children}</svg>`;
}

/* An SVG text label; extra carries anchors, fills and transforms. */
export const svgText = (x, y, str, cls = "lbl", extra = {}) =>
  html`<text x=${x} y=${y} class=${cls} ...${extra}>${str}</text>`;

/* The 0–100 % gridlines with their tick labels, shared by the coverage charts. */
export const percentGrid = (y, padL, right) => [0, 25, 50, 75, 100].map((v) => html`
  <line class="grid-line" x1=${padL} x2=${right} y1=${y(v)} y2=${y(v)} />
  ${svgText(padL - 6, y(v) + 3.5, v, "tick", { "text-anchor": "end" })}`);
