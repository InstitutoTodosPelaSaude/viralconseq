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

// An ISO 8601 timestamp, as the run manifest and versions.tsv write them.
const ISO_TIMESTAMP = /^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}(:\d{2}(\.\d+)?)?(Z|[+-]\d{2}:?\d{2})?$/;

function isTimestamp(value) {
  return typeof value === "string" && ISO_TIMESTAMP.test(value);
}

// A timestamp in the viewer's locale and time zone, named so the hour is unambiguous. Fractions of a second
// are cut to milliseconds first: Python writes microseconds, which not every browser parses.
function formatTimestamp(value) {
  const date = new Date(value.replace(/(\.\d{3})\d+/, "$1"));
  if (Number.isNaN(date.getTime())) return value;
  return date.toLocaleString(undefined, {
    year: "numeric", month: "2-digit", day: "2-digit", hour: "2-digit", minute: "2-digit", timeZoneName: "short",
  });
}

function truncate(text, maxLength = 18) {
  return text && text.length > maxLength ? text.slice(0, maxLength - 1) + "…" : text;
}

const HTML_ESCAPES = { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" };

function escapeHtml(text) {
  return String(text).replace(/[&<>"]/g, (character) => HTML_ESCAPES[character]);
}
