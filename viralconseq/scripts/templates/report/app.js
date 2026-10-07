/* The page: the state the sections share, and their order. */
import { render } from "preact";
import { useEffect, useRef, useState } from "preact/hooks";
import { html, S, ranked } from "report/core";
import { Masthead } from "report/masthead";
import { Kpis } from "report/kpis";
import { Overview } from "report/overview";
import { SampleTable, visibleRows } from "report/table";
import { SampleDetail } from "report/detail";
import { Provenance } from "report/provenance";

function App() {
  /* Charts read the theme's colours as they draw: bumping this redraws them. */
  const [, setColours] = useState(0);
  const [view, setView] = useState({ filter: "all", query: "", sortKey: "cov_min", sortDir: -1 });
  const [selected, setSelected] = useState((ranked[0] ?? S[0])?.key);
  /* Selecting from the table or a chart scrolls the detail into view; Previous/Next do not. */
  const [scrolls, setScrolls] = useState(0);
  const detail = useRef(null);
  useEffect(() => {
    if (scrolls) detail.current.scrollIntoView({ behavior: "smooth", block: "nearest" });
  }, [scrolls]);

  const rows = visibleRows(view);
  const select = (key, scroll = true) => {
    setSelected(key);
    if (scroll) setScrolls((n) => n + 1);
  };

  return html`
    <${Masthead} onTheme=${() => setColours((n) => n + 1)} />
    <div class="wrap">
      <${Kpis} />
      <${Overview} onSelect=${select} />
      <${SampleTable} view=${view} setView=${setView} rows=${rows} selected=${selected} onSelect=${select} />
      <section>
        <div class="sec-head">
          <h2>Sample detail</h2>
          <p>Select any row above.</p>
        </div>
        <div class="card" id="detail" ref=${detail}>
          <${SampleDetail} selected=${selected} rows=${rows} onSelect=${select} />
        </div>
      </section>
      <${Provenance} />
    </div>`;
}

/* The root holds a note for browsers that cannot run the page; render() adds to what is
   already there, so clear it first. */
const root = document.getElementById("app");
root.replaceChildren();
render(html`<${App} />`, root);
