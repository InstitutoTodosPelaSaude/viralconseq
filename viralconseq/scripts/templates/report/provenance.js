/* The footer: tool versions, run settings and where every figure came from. */
import { html, DATA, IS_ILLUMINA, P, PLATFORM, T, compact, fmt } from "report/core";

/* versions longer than this go to the fine print instead of the grid */
const LONG = 26;

const provRow = ([k, v]) => v === null || v === undefined || v === "" ? null
  : html`<div class="row"><span>${k}</span><b>${String(v)}</b></div>`;

function viralqcLine(vs) {
  if (!vs) return null;
  return [vs.status, vs.exit_code !== undefined ? `exit ${vs.exit_code}` : null,
    vs.result_rows !== undefined && vs.input_records !== undefined ? `${vs.result_rows}/${vs.input_records} records` : null]
    .filter((v) => v !== null && v !== undefined).join(", ");
}

export function Provenance() {
  const M = DATA.run.manifest;
  const versions = DATA.run.versions;
  const settings = [
    ["platform", PLATFORM],
    ["minimum depth", `${P.min_depth}×`],
    ["consensus coverage threshold", `${T}% at ≥${P.min_depth}×`],
    ["primer scheme", P.primer_scheme],
    [IS_ILLUMINA ? "allele-frequency threshold" : "allele-fraction threshold (ALT / (REF + ALT))", P.af_threshold],
    ["minimum read length", P.minimum_length === null ? null : `${fmt(P.minimum_length)} bp`],
    ["MAPQ floor", P.minimum_map_quality],
    ["Clair3 model", IS_ILLUMINA ? null : P.clair3_model],
    ["viralQC", viralqcLine(DATA.run.viralqc_status)],
    ["run created (UTC)", M ? M.created_utc : null],
    ["run finished (UTC)", M ? M.finished_utc : null],
    ["run status", M ? M.status : null],
    ["config sha256", M && M.config_sha256 ? M.config_sha256.slice(0, 12) : null],
  ];
  const inputs = M && M.inputs && M.inputs.length ? M.inputs : [];
  return html`
    <footer>
      <div class="sec-head"><h2>Provenance</h2><p>Everything the page reports comes from these versions and settings.</p></div>
      <div class="prov" id="prov">
        ${versions.filter(([, v]) => v.length <= LONG).map(provRow)}
        ${settings.map(provRow)}
      </div>
      <p class="fine" id="fine">
        Read counts and coverage come from <code>summary.tsv</code>; the depth trace from each
        row's <code>table_cov_basewise.txt</code>; variant counts are VCF records; the N mask is derived from
        depth ${"<"} ${P.min_depth}×. Clade, lineage and grade are viralQC's. Source: <code>${DATA.run.output_dir}</code>.
        ${versions.filter(([, v]) => v.length > LONG).map(([k, v]) => html`<br /><code>${k}</code>: ${v}`)}
        ${inputs.length ? html`<br />Inputs:` : null}
        ${inputs.map((i) => html`<br /><code>${i.sample}</code> ${i.path ? String(i.path).split("/").pop() : "—"} · sha256 ${
          i.sha256 || "—"}${i.size_bytes != null ? ` · ${compact(i.size_bytes)} B` : ""}`)}
      </p>
    </footer>`;
}
