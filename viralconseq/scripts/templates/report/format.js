"use strict";

const MISSING = "—";

function isMissing(value) {
  return value === null || value === undefined;
}

function formatNumber(value) {
  return isMissing(value) ? MISSING : value.toLocaleString("en-US");
}

function formatDecimal(value) {
  return isMissing(value) ? MISSING : value.toLocaleString("en-US", { maximumFractionDigits: 1 });
}

function formatPercent(value) {
  return isMissing(value) ? MISSING : value.toFixed(1) + "%";
}

function formatFraction(fraction) {
  return isMissing(fraction) ? MISSING : formatPercent(100 * fraction);
}

function formatShare(part, whole) {
  return isMissing(part) || !whole ? MISSING : formatPercent(100 * part / whole);
}

function formatCompact(value) {
  if (value >= 1e7) return (value / 1e6).toFixed(0) + "M";
  if (value >= 1e6) return (value / 1e6).toFixed(1) + "M";
  if (value >= 1e4) return Math.round(value / 1e3) + "K";
  return value.toLocaleString("en-US");
}

function truncate(text, maxLength = 18) {
  return text && text.length > maxLength ? text.slice(0, maxLength - 1) + "…" : text;
}

const HTML_ESCAPES = { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" };

function escapeHtml(text) {
  return String(text).replace(/[&<>"]/g, (character) => HTML_ESCAPES[character]);
}
