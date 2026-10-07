/* The page header: run name, headline figures and the theme toggle. */
import { useEffect } from "preact/hooks";
import { html, DATA, IDS, PLATFORM, T, complete, isDark, version } from "report/core";

document.title = `${DATA.run.label} — viralconseq run report`;

/* onTheme is called whenever the colours change, so the charts can redraw with them:
   on a click here, and when the OS switches between light and dark. */
export function Masthead({ onTheme }) {
  useEffect(() => {
    const os = matchMedia("(prefers-color-scheme: dark)");
    os.addEventListener("change", onTheme);
    return () => os.removeEventListener("change", onTheme);
  }, []);
  /* OS preference by default, explicit once clicked */
  const toggle = () => {
    document.documentElement.dataset.theme = isDark() ? "light" : "dark";
    onTheme();
  };
  const meta = [
    ["Platform", html`<span class="platform">${PLATFORM}</span>`],
    ["Samples", IDS.length],
    [`Genomes ≥${T}%`, complete.length],
    ["Pipeline", "viralconseq " + version],
  ];
  return html`
    <header class="masthead">
      <div class="wrap">
        <div>
          <h1 id="run-title">${DATA.run.label} <span>run report</span></h1>
          <p class="sub" id="run-sub">viralconseq consensus genomes · ${PLATFORM} · generated ${DATA.run.generated}</p>
        </div>
        <div class="mast-meta" id="mast-meta">${meta.map(([k, v]) => html`<div>${k}<b>${v}</b></div>`)}</div>
        <button class="theme-toggle" type="button" id="theme" aria-label="Switch colour theme"
          onClick=${toggle}>${isDark() ? "Light theme" : "Dark theme"}</button>
      </div>
    </header>`;
}
