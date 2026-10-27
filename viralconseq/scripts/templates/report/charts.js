"use strict";

Chart.defaults.font.family = getComputedStyle(document.body).fontFamily;

const LOG_AXIS_DECADES = [1, 10, 100, 1000, 10000, 100000];

function resolveColor(color) {
  return color.startsWith("var(") ? cssVar(color.slice(4, -1)) : color;
}

function hexChannels(color) {
  const hex = resolveColor(color).replace("#", "");
  if (hex.length !== 6) return null;
  return [0, 2, 4].map((offset) => parseInt(hex.slice(offset, offset + 2), 16));
}

function withOpacity(color, opacity) {
  const channels = hexChannels(color);
  return channels ? `rgba(${channels.join(", ")}, ${opacity})` : color;
}

function labelColorOn(fillVariable) {
  const channels = hexChannels(`var(${fillVariable})`);
  if (!channels) return "#fff";
  const linear = channels.map((channel) => {
    const value = channel / 255;
    return value <= 0.04045 ? value / 12.92 : ((value + 0.055) / 1.055) ** 2.4;
  });
  const luminance = 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2];
  const contrastWithBlack = (luminance + 0.05) / 0.05;
  return contrastWithBlack > 4.5 ? cssVar("--brand-deep") : "#fff";
}

function estimateTextWidth(text, fontSize) {
  return String(text).length * fontSize * 0.56;
}

function chartOptions(options) {
  const tooltip = {
    backgroundColor: resolveColor("var(--surface)"),
    borderColor: resolveColor("var(--border)"),
    borderWidth: 1,
    titleColor: resolveColor("var(--ink)"),
    bodyColor: resolveColor("var(--ink-2)"),
    padding: { x: 11, y: 8 },
    cornerRadius: 6,
    displayColors: false,
    titleFont: { weight: 700, size: 12 },
    bodyFont: { size: 12 },
  };
  const plugins = options.plugins || {};
  return {
    responsive: true,
    maintainAspectRatio: false,
    animation: { duration: 300 },
    ...options,
    plugins: { legend: { display: false }, ...plugins, tooltip: { ...tooltip, ...plugins.tooltip } },
  };
}

function axisOptions({ title, ticks, ...options } = {}) {
  return {
    border: { display: false },
    grid: { color: resolveColor("var(--grid)"), tickLength: 0 },
    ...options,
    ticks: { color: resolveColor("var(--muted)"), font: { size: 10.5 }, padding: 6, ...ticks },
    title: title ? { display: true, text: title, color: resolveColor("var(--muted)"), font: { size: 10.5 } } : undefined,
  };
}

function decadeTicks(max) {
  return (scale) => {
    scale.ticks = LOG_AXIS_DECADES.filter((value) => value <= max).map((value) => ({ value }));
  };
}

function thresholdLine(value, label) {
  const line = {
    type: "line", yMin: value, yMax: value, borderColor: resolveColor("var(--accent)"), borderWidth: 1, borderDash: [4, 3],
  };
  if (!label) return line;
  return {
    ...line,
    label: {
      display: true,
      content: label,
      position: "end",
      yAdjust: -9,
      padding: 0,
      backgroundColor: "transparent",
      color: resolveColor("var(--accent)"),
      font: { size: 10.5, weight: 600 },
    },
  };
}

function tooltipFor(titleAt, linesAt) {
  return {
    callbacks: {
      title: (items) => titleAt(items[0].dataIndex),
      label: (item) => linesAt(item.dataIndex),
    },
  };
}

function selectsRowOnClick(rows) {
  return {
    onClick: (event, items) => items.length && selectRow(rows[items[0].index].key),
    onHover: (event, items) => { event.native.target.style.cursor = items.length ? "pointer" : "default"; },
  };
}

function facetLine(row) {
  return FACET_LABEL && row.facet ? [`${FACET_LABEL}: ${row.facet}`] : [];
}

let mountedCharts = [];

function mountChart(host, spec) {
  host.style.position = "relative";
  host.style.height = spec.height + "px";
  const canvas = document.createElement("canvas");
  canvas.setAttribute("role", "img");
  canvas.setAttribute("aria-label", spec.label);
  host.replaceChildren(canvas);

  const mounted = { host, canvas, spec, chart: new Chart(canvas, spec.config()) };
  if (spec.resetsZoomOnDoubleClick) canvas.addEventListener("dblclick", () => mounted.chart.resetZoom());
  mountedCharts.push(mounted);
}

function rebuildCharts() {
  destroyDetachedCharts();
  for (const mounted of mountedCharts) {
    mounted.chart.destroy();
    mounted.chart = new Chart(mounted.canvas, mounted.spec.config());
  }
}

function destroyDetachedCharts() {
  const detached = mountedCharts.filter((mounted) => !mounted.host.isConnected);
  detached.forEach((mounted) => mounted.chart.destroy());
  mountedCharts = mountedCharts.filter((mounted) => mounted.host.isConnected);
}

const FUNNEL_CHART = {
  height: FUNNEL_STAGES.length * 58 + 4,
  label: "Samples with reads through coverage filtering",
  config: funnelConfig,
};

function funnelConfig() {
  const firstStageCount = FUNNEL_STAGES[0].count || 1;
  return {
    type: "bar",
    data: {
      labels: FUNNEL_STAGES.map((stage) => stage.label),
      datasets: [{
        data: FUNNEL_STAGES.map((stage) => stage.count),
        backgroundColor: FUNNEL_STAGES.map((_, index) => resolveColor(`var(${funnelColorVariable(index)})`)),
        barThickness: 16,
        borderRadius: { topRight: 4, bottomRight: 4 },
        borderSkipped: "start",
      }],
    },
    options: chartOptions({
      indexAxis: "y",
      layout: { padding: { right: 40 } },
      scales: { x: { display: false, min: 0, max: firstStageCount }, y: { display: false } },
      plugins: {
        annotation: { annotations: FUNNEL_STAGES.flatMap((stage, index) => funnelLabels(stage, index, firstStageCount)) },
        tooltip: tooltipFor(
          (index) => FUNNEL_STAGES[index].label,
          (index) => funnelTooltipLines(FUNNEL_STAGES[index], firstStageCount),
        ),
      },
    }),
  };
}

function funnelColorVariable(index) {
  const lastIndex = FUNNEL_STAGES.length - 1;
  const step = lastIndex > 0 ? Math.round(index * 3 / lastIndex) + 1 : 1;
  return `--ord-${step}`;
}

function funnelLabels(stage, index, firstStageCount) {
  const barWidth = document.getElementById("funnel").clientWidth - 46;
  const countFitsInBar = estimateTextWidth(stage.count, 11) + 18 < Math.max(3, barWidth * stage.count / firstStageCount);
  const share = stage.showsShare ? ` · ${(100 * stage.count / firstStageCount).toFixed(0)}% of samples` : "";
  const at = (xValue, content, options) => ({ type: "label", xValue, yValue: index, content, padding: 0, ...options });
  const leftAligned = { x: "start", y: "center" };
  const countFont = { size: 11, weight: 600 };
  return [
    at(0, stage.label, { position: leftAligned, yAdjust: -16, color: resolveColor("var(--ink)"), font: { size: 11.5, weight: 600 } }),
    at(0, stage.note + share, { position: leftAligned, yAdjust: 19, color: resolveColor("var(--muted)"), font: { size: 10.5 } }),
    countFitsInBar
      ? at(stage.count, String(stage.count), {
        position: { x: "end", y: "center" }, xAdjust: -8, color: labelColorOn(funnelColorVariable(index)), font: countFont,
      })
      : at(stage.count, String(stage.count), {
        position: leftAligned, xAdjust: 8, color: resolveColor("var(--ink)"), font: countFont,
      }),
  ];
}

function funnelTooltipLines(stage, firstStageCount) {
  return [
    `Samples: ${formatNumber(stage.count)}`,
    `Share of run: ${(100 * stage.count / firstStageCount).toFixed(1)}%`,
    `This stage keeps: ${stage.note}`,
  ];
}

const RANKED_CHART = {
  height: 186,
  label: `Genome completeness at ${MIN_DEPTH}x for each row, ranked`,
  config: rankedConfig,
};

function rankedConfig() {
  return {
    type: "bar",
    data: {
      labels: RANKED_ROWS.map((row) => row.key),
      datasets: [{
        data: RANKED_ROWS.map((row) => row.cov_min),
        backgroundColor: RANKED_ROWS.map((row) => resolveColor(facetColor(row.facet))),
        ...COLUMN_BARS,
      }],
    },
    options: chartOptions({
      interaction: { mode: "index", intersect: false },
      ...selectsRowOnClick(RANKED_ROWS),
      scales: {
        x: axisOptions({ grid: { display: false }, ticks: { display: false } }),
        y: axisOptions({ min: 0, max: 100, ticks: { stepSize: 25 } }),
      },
      plugins: {
        annotation: { annotations: { threshold: thresholdLine(COVERAGE_THRESHOLD, `${COVERAGE_THRESHOLD}% threshold`) } },
        tooltip: tooltipFor((index) => RANKED_ROWS[index].key, (index) => rankedTooltipLines(RANKED_ROWS[index])),
      },
    }),
  };
}

// minBarLength keeps a zero visible, so the row that failed still shows.
const COLUMN_BARS = {
  maxBarThickness: 24,
  minBarLength: 3,
  categoryPercentage: 0.92,
  barPercentage: 1,
  borderRadius: { topLeft: 3, topRight: 3 },
  borderSkipped: "start",
};

function rankedTooltipLines(row) {
  return [
    ...facetLine(row),
    `≥${MIN_DEPTH}× coverage: ${formatPercent(row.cov_min)}`,
    `10× coverage: ${formatPercent(row.cov10)}`,
    `Median depth: ${formatDecimal(row.median_depth)}×`,
  ];
}

const SCATTER_ROWS = ROWS.filter((row) => row.median_depth !== null && row.cov_min !== null);

const SCATTER_CHART = {
  height: 220,
  label: "Median depth against genome completeness for each row",
  config: scatterConfig,
};

function scatterConfig() {
  const deepest = Math.max(...SCATTER_ROWS.map((row) => row.median_depth), 10);
  const xMax = 10 ** (Math.log10(deepest) * 1.04);
  return {
    type: "scatter",
    data: {
      datasets: [{
        data: SCATTER_ROWS.map((row) => ({ x: Math.max(1, row.median_depth), y: row.cov_min })),
        backgroundColor: SCATTER_ROWS.map((row) => resolveColor(facetColor(row.facet))),
        borderColor: resolveColor("var(--surface)"),
        borderWidth: 2,
        pointRadius: 6.5,
        pointHoverRadius: 7.5,
      }],
    },
    options: chartOptions({
      ...selectsRowOnClick(SCATTER_ROWS),
      scales: {
        x: axisOptions({
          type: "logarithmic", min: 1, max: xMax, title: "Median depth (log scale)",
          afterBuildTicks: decadeTicks(xMax),
          ticks: { callback: (value) => formatCompact(value) + "×" },
        }),
        y: axisOptions({ min: 0, max: 100, ticks: { stepSize: 25 }, title: `% covered at ≥${MIN_DEPTH}×` }),
      },
      plugins: {
        annotation: { annotations: { threshold: thresholdLine(COVERAGE_THRESHOLD, `${COVERAGE_THRESHOLD}%`) } },
        tooltip: tooltipFor((index) => SCATTER_ROWS[index].key, (index) => scatterTooltipLines(SCATTER_ROWS[index])),
      },
    }),
  };
}

function scatterTooltipLines(row) {
  const clade = (row.viralqc && (row.viralqc.clade || row.viralqc.lineage)) || MISSING;
  return [
    ...facetLine(row),
    `Median depth: ${formatDecimal(row.median_depth)}×`,
    `≥${MIN_DEPTH}× coverage: ${formatPercent(row.cov_min)}`,
    `Clade: ${clade}`,
  ];
}

function depthTraceChart(row) {
  return {
    height: 200,
    label: `Depth along the genome for ${row.key}`,
    config: () => depthTraceConfig(row),
    resetsZoomOnDoubleClick: true,
  };
}

// A log axis has no zero, so on it depth is plotted as depth + 1 and its ticks read back as depth.
function depthTraceConfig(row) {
  const coverage = coverageOf(row);
  const genomeLength = coverage.genome_length;
  const isLog = Alpine.store("ui").logDepth;
  const offset = isLog ? 1 : 0;
  const plottedMax = Math.max(10, ...coverage.bins) + offset;
  const color = resolveColor(facetColor(row.facet));
  const yAxis = isLog
    ? axisOptions({
      type: "logarithmic", min: 1, max: plottedMax,
      afterBuildTicks: decadeTicks(plottedMax),
      ticks: { callback: (value) => (value === 1 ? "0" : formatCompact(value)) },
    })
    : axisOptions({ type: "linear", min: 0, ticks: { maxTicksLimit: 5, callback: formatCompact } });
  return {
    type: "line",
    data: {
      datasets: [{
        data: coverage.bins.map((depth, index) => ({ x: index * coverage.bin_bp, y: depth + offset })),
        borderColor: color,
        borderWidth: 2,
        backgroundColor: withOpacity(color, 0.1),
        fill: "origin",
        pointRadius: 0,
        pointHoverRadius: 4.5,
        pointHoverBackgroundColor: color,
        pointHoverBorderColor: resolveColor("var(--surface)"),
        pointHoverBorderWidth: 2,
      }],
    },
    options: chartOptions({
      parsing: false,
      normalized: true,
      interaction: { mode: "index", intersect: false, axis: "x" },
      scales: {
        x: axisOptions({
          type: "linear", min: 0, max: genomeLength, grid: { display: false },
          title: `Genome position (bp)${contigSummary(coverage.contigs || [])}`,
          ticks: { callback: formatCompact, maxRotation: 0, includeBounds: false },
        }),
        y: yAxis,
      },
      plugins: {
        annotation: { annotations: depthTraceAnnotations(row, coverage, plottedMax, offset) },
        zoom: {
          zoom: { wheel: { enabled: true }, pinch: { enabled: true }, mode: "x" },
          pan: { enabled: true, mode: "x" },
          limits: { x: { min: 0, max: genomeLength, minRange: Math.max(10, coverage.bin_bp * 10) } },
        },
        tooltip: tooltipFor(
          (index) => positionRange(index, coverage),
          (index) => depthTraceTooltipLines(row, coverage, index),
        ),
      },
    }),
  };
}

function contigSummary(contigs) {
  if (contigs.length === 1) return ` · ${truncate(contigs[0].name, 40)}`;
  if (contigs.length > 1) return ` · ${contigs.length} contigs, concatenated`;
  return "";
}

function depthTraceAnnotations(row, coverage, plottedMax, offset) {
  const genomeLength = coverage.genome_length;
  const extent = primerExtentOf(row);
  const contigs = coverage.contigs || [];
  const band = (from, to, colorVariable) => ({
    type: "box", xMin: from, xMax: to, backgroundColor: resolveColor(colorVariable), borderWidth: 0,
    drawTime: "beforeDatasetsDraw",
  });
  const outsidePrimers = extent
    ? [band(0, extent.start, "var(--outside-band)"), band(extent.end, genomeLength, "var(--outside-band)")]
    : [];
  const masked = (coverage.mask || []).map(([first, last]) => band(first - 1, last, "var(--mask-band)"));
  const contigMarkers = contigs.length > 1 ? contigBoundaries(contigs).concat(contigNames(contigs, genomeLength, plottedMax)) : [];
  return [...outsidePrimers, ...masked, ...contigMarkers, thresholdLine(MIN_DEPTH + offset, `${MIN_DEPTH}× min`)];
}

function contigBoundaries(contigs) {
  return contigs.slice(1).map((contig) => ({
    type: "line", xMin: contig.offset, xMax: contig.offset,
    borderColor: resolveColor("var(--muted)"), borderWidth: 1, borderDash: [3, 3],
  }));
}

function contigNames(contigs, genomeLength, plottedMax) {
  return contigs.map((contig) => ({
    type: "label",
    xValue: contig.offset + contig.length / 2,
    yValue: plottedMax,
    content: truncate(contig.name, Math.max(4, Math.floor(contig.length / genomeLength * 600 / 6.5))),
    position: { x: "center", y: "start" },
    yAdjust: 3,
    padding: { x: 4, y: 1 },
    backgroundColor: withOpacity("var(--surface)", 0.85),
    color: resolveColor("var(--muted)"),
    font: { size: 10.5 },
  }));
}

function positionRange(binIndex, coverage) {
  const start = binIndex * coverage.bin_bp;
  return `Position ${formatNumber(start + 1)}–${formatNumber(Math.min(coverage.genome_length, start + coverage.bin_bp))}`;
}

function depthTraceTooltipLines(row, coverage, binIndex) {
  const position = binIndex * coverage.bin_bp;
  const contigs = coverage.contigs || [];
  const extent = primerExtentOf(row);
  const isMasked = (coverage.mask || []).some(([first, last]) => position + 1 >= first && position + 1 <= last);
  const contig = contigs.length > 1
    ? contigs.find((candidate) => position >= candidate.offset && position < candidate.offset + candidate.length)
    : null;
  const isOutsidePrimers = extent && (position < extent.start || position > extent.end);
  return [
    ...(contig ? [`Contig: ${contig.name}`, `Contig position: ${formatNumber(position - contig.offset + 1)}`] : []),
    `Mean depth: ${formatDecimal(coverage.bins[binIndex])}×`,
    `Consensus: ${isMasked ? "masked to N" : "called"}`,
    ...(extent ? [`Primer extent: ${isOutsidePrimers ? "outside" : "inside"}`] : []),
  ];
}

function ampliconDepthChart(row) {
  return {
    height: 150,
    label: `Mean depth per amplicon for ${row.key}`,
    config: () => ampliconDepthConfig(row),
  };
}

function ampliconDepthConfig(row) {
  const depths = ampliconDepthsOf(row);
  const color = resolveColor(facetColor(row.facet));
  const belowColor = resolveColor("var(--accent)");
  const everyFifthAmplicon = (value, index) => (index === 0 || (index + 1) % 5 === 0 ? index + 1 : "");
  return {
    type: "bar",
    data: {
      labels: depths.map((_, index) => index + 1),
      datasets: [{
        data: depths,
        backgroundColor: depths.map((depth) => (depth < MIN_DEPTH ? belowColor : color)),
        ...COLUMN_BARS,
      }],
    },
    options: chartOptions({
      interaction: { mode: "index", intersect: false },
      scales: {
        x: axisOptions({ grid: { display: false }, ticks: { autoSkip: false, maxRotation: 0, callback: everyFifthAmplicon } }),
        y: axisOptions({ min: 0, ticks: { maxTicksLimit: 5, callback: formatCompact } }),
      },
      plugins: {
        annotation: { annotations: { threshold: thresholdLine(MIN_DEPTH) } },
        tooltip: tooltipFor((index) => `Amplicon ${index + 1}`, (index) => ampliconTooltipLines(depths[index])),
      },
    }),
  };
}

function ampliconTooltipLines(depth) {
  return [
    `Mean depth: ${formatDecimal(depth)}×`,
    `Against threshold: ${depth < MIN_DEPTH ? "below" : "at or above"} ${MIN_DEPTH}×`,
  ];
}
