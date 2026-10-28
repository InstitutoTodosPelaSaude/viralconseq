"use strict";

// Runs in <head>, before the page is first painted: the theme the page opens in (?theme=light|dark, else the
// viewer's stored choice, else the system's) and the embed mode (?sample=<id>), so neither changes after paint.
(() => {
  const params = new URLSearchParams(location.search);
  const urlTheme = params.get("theme");
  let storedTheme = null;
  try {
    storedTheme = localStorage.getItem("viralconseq.report.theme");
  } catch {
    // storage unavailable: the system's theme
  }
  const theme = ["light", "dark"].includes(urlTheme) ? urlTheme : storedTheme;
  if (theme === "light" || theme === "dark") document.documentElement.dataset.theme = theme;
  if (params.get("sample") !== null) document.documentElement.dataset.embed = "sample";
})();
