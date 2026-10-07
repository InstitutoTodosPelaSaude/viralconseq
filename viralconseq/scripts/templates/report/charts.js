"use strict";
/* The five charts, on Chart.js with its annotation plugin (thresholds, bands, labels) and
   zoom plugin (the depth trace). Each chart is a spec: { height, label, config() } where
   config() returns a Chart.js configuration. The template mounts a spec on its host
   element with x-init="mount($el, …)"; a theme change builds every chart again. */

/* ---------- colours and type: Chart.js paints concrete colours, so the stylesheet's
   custom properties are read each time a chart is built ---------- */
const paint = (c) => c.startsWith("var(") ? css(c.slice(4, -1)) : c;
/* a resolved #rrggbb colour at the given opacity */
function fade(c, alpha) {
  const hex = paint(c).replace("#", "");
  if (hex.length !== 6) return c;
  const [r, g, b] = [0, 2, 4].map((i) => parseInt(hex.slice(i, i + 2), 16));
  return `rgba(${r}, ${g}, ${b}, ${alpha})`;
}
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
/* Approximate rendered width — used to decide whether a label fits inside a bar. */
const textWidth = (str, size) => String(str).length * size * 0.56;

/* What every chart shares: type, tooltip, short animation, no legend (the page has its own). */
function base(extra = {}) {
  const tooltip = {
    backgroundColor: paint("var(--surface)"), borderColor: paint("var(--border)"), borderWidth: 1,
    titleColor: paint("var(--ink)"), bodyColor: paint("var(--ink-2)"), padding: { x: 11, y: 8 },
    cornerRadius: 6, displayColors: false, titleFont: { weight: 700, size: 12 }, bodyFont: { size: 12 },
  };
  Chart.defaults.font.family = getComputedStyle(document.body).fontFamily;
  return {
    responsive: true, maintainAspectRatio: false, animation: { duration: 300 },
    ...extra,
    plugins: { legend: { display: false }, ...extra.plugins, tooltip: { ...tooltip, ...(extra.plugins || {}).tooltip } },
  };
}
/* Axis chrome in the page's colours. */
const axis = (extra = {}) => ({
  border: { display: false },
  grid: { color: paint("var(--grid)"), tickLength: 0 },
  ...extra,
  ticks: { color: paint("var(--muted)"), font: { size: 10.5 }, padding: 6, ...extra.ticks },
  title: extra.title ? { display: true, color: paint("var(--muted)"), font: { size: 10.5 }, text: extra.title } : undefined,
});
/* A threshold rule across the plot, labelled at its right end. */
const threshold = (value, label) => ({
  type: "line", yMin: value, yMax: value, borderColor: paint("var(--accent)"), borderWidth: 1,
  label: label ? { display: true, content: label, position: "end", yAdjust: -9, padding: 0,
    backgroundColor: "transparent", color: paint("var(--accent)"), font: { size: 10.5, weight: 600 } } : undefined,
});

/* ---------- mounting: one chart per host ---------- */
let mounted = [];
function build(entry) {
  entry.chart = new Chart(entry.canvas, entry.spec.config());
}
function mount(host, spec) {
  host.style.position = "relative";
  host.style.height = spec.height + "px";
  const canvas = document.createElement("canvas");
  canvas.setAttribute("role", "img");
  canvas.setAttribute("aria-label", spec.label);
  host.replaceChildren(canvas);
  const entry = { host, canvas, spec };
  build(entry);
  if (spec.mounted) spec.mounted(entry.chart, canvas);
  mounted.push(entry);
}
/* Detail charts are replaced with every selection: destroy the ones whose host has
   left the document. */
function dropDetached() {
  mounted = mounted.filter((e) => e.host.isConnected || (e.chart.destroy(), false));
}
/* Colours are read when a chart is built, so a theme change builds them all again. */
function redrawAll() {
  dropDetached();
  for (const e of mounted) { e.chart.destroy(); build(e); }
}

/* ---------- funnel: unique samples surviving each stage ---------- */
const FUNNEL_CHART = {
  height: FUNNEL_STAGES.length * 58 + 4,
  label: "Samples with reads through coverage filtering",
  config() {
    const top = FUNNEL_STAGES[0].n || 1, n = FUNNEL_STAGES.length;
    const ord = (i) => `--ord-${n > 1 ? Math.round(i * 3 / (n - 1)) + 1 : 1}`;
    const width = document.getElementById("funnel").clientWidth - 46;
    /* stage name above each bar, its note below, its count inside or after it */
    const labels = FUNNEL_STAGES.flatMap((st, i) => {
      const inside = textWidth(st.n, 11) + 18 < Math.max(3, width * st.n / top);
      const note = st.note + (st.share ? ` · ${(100 * st.n / top).toFixed(0)}% of samples` : "");
      const at = (x, content, extra) => ({ type: "label", xValue: x, yValue: i, content, padding: 0, ...extra });
      return [
        at(0, st.label, { position: { x: "start", y: "center" }, yAdjust: -16, color: paint("var(--ink)"), font: { size: 11.5, weight: 600 } }),
        at(0, note, { position: { x: "start", y: "center" }, yAdjust: 19, color: paint("var(--muted)"), font: { size: 10.5 } }),
        at(st.n, String(st.n), inside
          ? { position: { x: "end", y: "center" }, xAdjust: -8, color: inkOn(ord(i)), font: { size: 11, weight: 600 } }
          : { position: { x: "start", y: "center" }, xAdjust: 8, color: paint("var(--ink)"), font: { size: 11, weight: 600 } }),
      ];
    });
    return {
      type: "bar",
      data: { labels: FUNNEL_STAGES.map((st) => st.label), datasets: [{
        data: FUNNEL_STAGES.map((st) => st.n), barThickness: 16,
        backgroundColor: FUNNEL_STAGES.map((_, i) => paint(`var(${ord(i)})`)),
        borderRadius: { topRight: 4, bottomRight: 4 }, borderSkipped: "start",
      }] },
      options: base({
        indexAxis: "y",
        layout: { padding: { right: 40 } },
        scales: { x: { display: false, min: 0, max: top }, y: { display: false } },
        plugins: {
          annotation: { annotations: labels },
          tooltip: { callbacks: {
            title: (items) => FUNNEL_STAGES[items[0].dataIndex].label,
            label: (item) => { const st = FUNNEL_STAGES[item.dataIndex];
              return [`Samples: ${fmt(st.n)}`, `Share of run: ${(100 * st.n / top).toFixed(1)}%`, `This stage keeps: ${st.note}`]; },
          } },
        },
      }),
    };
  },
};

/* ---------- ranked completeness: one bar per row, best first ---------- */
const RANKED_CHART = {
  height: 186,
  label: `Genome completeness at ${P.min_depth}x for each row, ranked`,
  config: () => ({
    type: "bar",
    data: { labels: ranked.map((d) => d.key), datasets: [{
      data: ranked.map((d) => d.cov_min), maxBarThickness: 24, minBarLength: 3,
      categoryPercentage: 0.92, barPercentage: 1,
      backgroundColor: ranked.map((d) => paint(facetColor(d.facet))), borderRadius: { topLeft: 3, topRight: 3 }, borderSkipped: "start",
    }] },
    options: base({
      /* the whole column under a bar is its hover and click target */
      interaction: { mode: "index", intersect: false },
      onClick: (e, items) => items.length && select(ranked[items[0].index].key),
      onHover: (e, items) => { e.native.target.style.cursor = items.length ? "pointer" : "default"; },
      scales: {
        x: axis({ grid: { display: false }, ticks: { display: false } }),
        y: axis({ min: 0, max: 100, ticks: { stepSize: 25 } }),
      },
      plugins: {
        annotation: { annotations: { threshold: threshold(T, `${T}% threshold`) } },
        tooltip: { callbacks: {
          title: (items) => ranked[items[0].dataIndex].key,
          label: (item) => { const d = ranked[item.dataIndex];
            return [...facetRows(d).map(([k, v]) => `${k}: ${d.facet}`), `≥${P.min_depth}× coverage: ${pct(d.cov_min)}`,
                    `10× coverage: ${pct(d.cov10)}`, `Median depth: ${fmt1(d.median_depth)}×`]; },
        } },
      },
    }),
  }),
};

/* ---------- depth against completeness: one dot per row, log depth ---------- */
const SCATTER_POINTS = S.filter((d) => d.median_depth !== null && d.cov_min !== null);
const DECADES = [1, 10, 100, 1000, 10000, 100000];
const SCATTER_CHART = {
  height: 220,
  label: "Median depth against genome completeness for each row",
  config() {
    const xMax = 10 ** (Math.log10(Math.max(...SCATTER_POINTS.map((d) => d.median_depth), 10)) * 1.04);
    return {
      type: "scatter",
      data: { datasets: [{
        data: SCATTER_POINTS.map((d) => ({ x: Math.max(1, d.median_depth), y: d.cov_min })),
        /* 2px surface ring keeps overlapping dots legible */
        backgroundColor: SCATTER_POINTS.map((d) => paint(facetColor(d.facet))),
        borderColor: paint("var(--surface)"), borderWidth: 2, pointRadius: 6.5, pointHoverRadius: 7.5,
      }] },
      options: base({
        onClick: (e, items) => items.length && select(SCATTER_POINTS[items[0].index].key),
        onHover: (e, items) => { e.native.target.style.cursor = items.length ? "pointer" : "default"; },
        scales: {
          x: axis({ type: "logarithmic", min: 1, max: xMax, title: "Median depth (log scale)",
                    afterBuildTicks: (scale) => { scale.ticks = DECADES.filter((v) => v <= xMax).map((value) => ({ value })); },
                    ticks: { callback: (v) => compact(v) + "×" } }),
          y: axis({ min: 0, max: 100, ticks: { stepSize: 25 }, title: `% covered at ≥${P.min_depth}×` }),
        },
        plugins: {
          annotation: { annotations: { threshold: threshold(T, `${T}%`) } },
          tooltip: { callbacks: {
            title: (items) => SCATTER_POINTS[items[0].dataIndex].key,
            label: (item) => { const d = SCATTER_POINTS[item.dataIndex];
              return [...facetRows(d).map(([k]) => `${k}: ${d.facet}`), `Median depth: ${fmt1(d.median_depth)}×`,
                      `≥${P.min_depth}× coverage: ${pct(d.cov_min)}`,
                      `Clade: ${(d.viralqc && (d.viralqc.clade || d.viralqc.lineage)) || "—"}`]; },
          } },
        },
      }),
    };
  },
};

/* ---------- depth along the genome, for row d: wheel to zoom, drag to pan ---------- */
function traceChart(d) {
  const cov = DATA.coverage[d.key];
  const ext = extentOf(d);
  const L = cov.genome_length;
  const contigs = cov.contigs || [];
  const mask = cov.mask || [];
  const maxD = Math.max(10, ...cov.bins);
  /* a log axis has no zero: plot depth + 1, as the y = log10(depth + 1) trace always did */
  const points = cov.bins.map((v, i) => ({ x: i * cov.bin_bp, y: v + 1 }));
  const band = (from, to, colour) => ({ type: "box", xMin: from, xMax: to, backgroundColor: colour, borderWidth: 0,
                                        drawTime: "beforeDatasetsDraw" });
  return {
    height: 200,
    label: `Depth along the genome for ${d.key}`,
    config() {
      const color = paint(facetColor(d.facet));
      return {
        type: "line",
        data: { datasets: [{
          data: points, borderColor: color, borderWidth: 2, backgroundColor: fade(color, 0.1), fill: "origin",
          pointRadius: 0, pointHoverRadius: 4.5, pointHoverBackgroundColor: color,
          pointHoverBorderColor: paint("var(--surface)"), pointHoverBorderWidth: 2,
        }] },
        options: base({
          parsing: false, normalized: true,
          interaction: { mode: "index", intersect: false, axis: "x" },
          scales: {
            x: axis({ type: "linear", min: 0, max: L, grid: { display: false },
                      title: `Genome position (bp)${contigs.length === 1 ? ` · ${short(contigs[0].name, 40)}`
                        : contigs.length > 1 ? ` · ${contigs.length} contigs, concatenated` : ""}`,
                      ticks: { callback: (v) => compact(v), maxRotation: 0, includeBounds: false } }),
            y: axis({ type: "logarithmic", min: 1, max: maxD + 1,
                      afterBuildTicks: (scale) => { scale.ticks = DECADES.filter((v) => v <= maxD + 1).map((value) => ({ value })); },
                      ticks: { callback: (v) => v === 1 ? "0" : compact(v) } }),
          },
          plugins: {
            annotation: { annotations: [
              /* outside the primer extent, then the masked ranges */
              ...(ext ? [band(0, ext.start, paint("var(--surface-sunk)")), band(ext.end, L, paint("var(--surface-sunk)"))] : []),
              ...mask.map(([a, b]) => band(a - 1, b, paint("var(--grid)"))),
              /* contig boundaries and names */
              ...contigs.slice(1).map((c) => ({ type: "line", xMin: c.offset, xMax: c.offset, borderColor: paint("var(--muted)"),
                                                 borderWidth: 1, borderDash: [3, 3] })),
              ...(contigs.length > 1 ? contigs.map((c) => ({ type: "label", xValue: c.offset + c.length / 2, yValue: maxD + 1,
                content: short(c.name, Math.max(4, Math.floor(c.length / L * 600 / 6.5))), position: { x: "center", y: "start" },
                yAdjust: 3, padding: { x: 4, y: 1 }, backgroundColor: fade(paint("var(--surface)"), 0.85),
                color: paint("var(--muted)"), font: { size: 10.5 } })) : []),
              threshold(P.min_depth + 1, `${P.min_depth}×`),
            ] },
            zoom: {
              zoom: { wheel: { enabled: true }, pinch: { enabled: true }, mode: "x" },
              pan: { enabled: true, mode: "x" },
              limits: { x: { min: 0, max: L, minRange: Math.max(10, cov.bin_bp * 10) } },
            },
            tooltip: { callbacks: {
              title: (items) => { const pos = items[0].dataIndex * cov.bin_bp;
                return `Position ${fmt(pos + 1)}–${fmt(Math.min(L, pos + cov.bin_bp))}`; },
              label: (item) => {
                const i = item.dataIndex, pos = i * cov.bin_bp;
                const masked = mask.some(([a, b]) => pos + 1 >= a && pos + 1 <= b);
                const contig = contigs.length > 1 ? contigs.find((c) => pos >= c.offset && pos < c.offset + c.length) : null;
                return [...(contig ? [`Contig: ${contig.name}`, `Contig position: ${fmt(pos - contig.offset + 1)}`] : []),
                        `Mean depth: ${fmt1(cov.bins[i])}×`, `Consensus: ${masked ? "masked to N" : "called"}`,
                        ...(ext ? [`Primer extent: ${pos < ext.start || pos > ext.end ? "outside" : "inside"}`] : [])];
              },
            } },
          },
        }),
      };
    },
    /* double-click returns to the whole genome */
    mounted: (chart, canvas) => canvas.addEventListener("dblclick", () => {
      const live = mounted.find((e) => e.canvas === canvas);
      if (live) live.chart.resetZoom();
    }),
  };
}

/* ---------- mean depth per amplicon, for row d ---------- */
function ampliconChart(d) {
  const depths = DATA.amplicons[d.key];
  return {
    height: 150,
    label: `Mean depth per amplicon for ${d.key}`,
    config() {
      const color = paint(facetColor(d.facet)), accent = paint("var(--accent)");
      return {
        type: "bar",
        data: { labels: depths.map((_, i) => i + 1), datasets: [{
          /* minBarLength keeps a true zero visible, so the failed amplicon still shows */
          data: depths, maxBarThickness: 24, minBarLength: 3, categoryPercentage: 0.92, barPercentage: 1,
          backgroundColor: depths.map((v) => v < P.min_depth ? accent : color),
          borderRadius: { topLeft: 3, topRight: 3 }, borderSkipped: "start",
        }] },
        options: base({
          interaction: { mode: "index", intersect: false },
          scales: {
            x: axis({ grid: { display: false }, ticks: { autoSkip: false, maxRotation: 0,
                      callback: (v, i) => i === 0 || (i + 1) % 5 === 0 ? i + 1 : "" } }),
            y: axis({ min: 0, ticks: { maxTicksLimit: 5, callback: (v) => compact(v) } }),
          },
          plugins: {
            annotation: { annotations: { threshold: threshold(P.min_depth) } },
            tooltip: { callbacks: {
              title: (items) => `Amplicon ${items[0].dataIndex + 1}`,
              label: (item) => { const v = depths[item.dataIndex];
                return [`Mean depth: ${fmt1(v)}×`, `Against threshold: ${v < P.min_depth ? "below" : "at or above"} ${P.min_depth}×`]; },
            } },
          },
        }),
      };
    },
  };
}
