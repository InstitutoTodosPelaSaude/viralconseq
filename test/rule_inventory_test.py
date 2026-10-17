"""Static conventions over the Snakemake rule files.

Pure text checks, no Snakemake import: the four entry workflows under
``viralconseq/scripts/`` and everything they ``include:`` are scanned rule by
rule. The conventions enforced here are the ones that, when broken, hide
failures rather than cause them (a rule without a log, a shell body that never
writes to its log) or move resources the operator never named (a rule reading
another rule's ``_cpus`` key).
"""

import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from typing import Dict, Iterator, List, Tuple

import yaml

import viralconseq
from viralconseq import integrity
from viralconseq.constants import ResourceDefaults
from viralconseq.scripts.python import build_report

SCRIPTS_DIR = Path(viralconseq.__file__).resolve().parent / "scripts"
ENTRY_FILES = sorted(SCRIPTS_DIR.glob("consensus_*.smk"))

RULE_RE = re.compile(r"^rule (\w+):", re.MULTILINE)
INCLUDE_RE = re.compile(r'^include:\s*"([^"]+\.smk)"', re.MULTILINE)
# A rule body ends at the next line that starts in column 0.
TOP_LEVEL_RE = re.compile(r"^(?=\S)", re.MULTILINE)
# cpus("<rule>") / ram_mb("<rule>") from rules/common.smk, or a raw config read.
RESOURCE_KEY_RE = re.compile(r'(?:cpus|ram_mb)\("(\w+)"\)|config(?:\.get\(|\[)"(\w+)_(?:cpus|ram)"')
DIRECTIVE_RE = {
    name: re.compile(rf"^\s+{name}:", re.MULTILINE)
    for name in ("log", "benchmark", "shell", "script", "run")
}

# Target-only rules never execute anything, so a log would be empty.
LOG_EXEMPT = {"all"}

# Analysis parameters the CLI always writes. Rules must read them as
# config["key"] (required), never config.get("key", <literal>): a literal
# fallback is a second, silently drifting default. Mirrors _REQUIRED_KEYS in
# rules/common.smk; keep both lists in step.
REQUIRED_COMMON = [
    "samples",
    "data",
    "output",
    "threads",
    "reference",
    "scheme",
    "minimum_depth",
    "normalise_depth",
    "minimum_length",
    "af_threshold",
    "consensus_coverage_threshold",
]
# plus run_viralqc_ram when run_viralqc is on, and infer_consensus_sequence_ram
# on nanopore (ConfigGenerator always writes them for the memory rules).
REQUIRED_ILLUMINA = [
    "adapters",
    "trim_head",
    "trim_tail",
    "cut_front_mean_quality",
    "cut_tail_mean_quality",
    "cut_right_window_size",
    "cut_right_mean_quality",
    "af_isnv_threshold",
]
REQUIRED_NANOPORE = [
    "infer_consensus_sequence_ram",
    "minimum_mapped_reads",
    "chunk_size",
    "clair3_model",
    "variant_quality",
    "frameshift_variant_quality",
    "variant_depth",
    "minimum_map_quality",
]
ANALYSIS_KEYS = set(REQUIRED_COMMON + REQUIRED_ILLUMINA + REQUIRED_NANOPORE) - {
    # sentinel-valued keys are legitimately read with .get
    "scheme",
    "adapters",
}
FALLBACK_RE = re.compile(r'config\.get\("(\w+)",')
DRYRUN_CONFIG_DIR = Path(__file__).resolve().parent / "dryrun_configs"


def workflow_files(entry: Path) -> List[Path]:
    """Return the entry file plus every transitively included ``.smk``."""
    seen: List[Path] = []
    queue = [entry]
    while queue:
        path = queue.pop(0)
        if path in seen:
            continue
        seen.append(path)
        text = path.read_text()
        for rel in INCLUDE_RE.findall(text):
            queue.append((path.parent / rel).resolve())
    return seen


def rule_blocks(text: str) -> Iterator[Tuple[str, str]]:
    """Yield ``(rule name, indented body)`` for every rule in ``text``."""
    matches = list(RULE_RE.finditer(text))
    for match in matches:
        body_start = match.end()
        nxt = TOP_LEVEL_RE.search(text, body_start + 1)
        body_end = nxt.start() if nxt else len(text)
        yield match.group(1), text[body_start:body_end]


def rules_for_entry(entry: Path) -> Dict[str, Tuple[Path, str]]:
    rules: Dict[str, Tuple[Path, str]] = {}
    for path in workflow_files(entry):
        for name, body in rule_blocks(path.read_text()):
            rules[name] = (path, body)
    return rules


def all_rules() -> Dict[Tuple[str, str], Tuple[Path, str]]:
    """``{(entry name, rule name): (defining file, body)}`` over all workflows."""
    out: Dict[Tuple[str, str], Tuple[Path, str]] = {}
    for entry in ENTRY_FILES:
        for name, located in rules_for_entry(entry).items():
            out[(entry.name, name)] = located
    return out


def has(directive: str, body: str) -> bool:
    return bool(DIRECTIVE_RE[directive].search(body))


class Test_RuleInventory(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rules = all_rules()
        assert cls.rules, "no rules found under viralconseq/scripts"

    def test_entry_files_present(self):
        names = {p.name for p in ENTRY_FILES}
        self.assertEqual(
            names,
            {
                "consensus_illumina.smk",
                "consensus_illumina_segmented.smk",
                "consensus_nanopore.smk",
                "consensus_nanopore_segmented.smk",
            },
        )

    def test_entry_files_include_common_first_and_default_to_all(self):
        for entry in ENTRY_FILES:
            text = entry.read_text()
            with self.subTest(workflow=entry.name):
                first_include = INCLUDE_RE.search(text)
                self.assertIsNotNone(first_include, "no include: found")
                self.assertEqual(
                    first_include.group(1),
                    "rules/common.smk",
                    "rules/common.smk must be the first include (global wildcard constraints)",
                )
                all_body = dict(rule_blocks(text))["all"]
                self.assertIn("default_target: True", all_body)

    def test_every_rule_declares_log_and_benchmark_via_the_helpers(self):
        """Paths come from LOG()/BENCH() with the rule's own name, so every log
        lands at logs/<rule>/[<segment>/]<target> and the benchmark collector
        can attribute it."""
        for (entry, name), (path, body) in self.rules.items():
            if name in LOG_EXEMPT:
                continue
            with self.subTest(workflow=entry, rule=name, file=path.name):
                self.assertTrue(has("log", body), "missing log:")
                self.assertTrue(has("benchmark", body), "missing benchmark:")
                self.assertIn(f'LOG("{name}"', body, "log: must use LOG(<own rule name>, ...)")
                self.assertIn(
                    f'BENCH("{name}"', body, "benchmark: must use BENCH(<own rule name>, ...)"
                )

    def test_shell_rules_fail_fast_and_write_their_log(self):
        for (entry, name), (path, body) in self.rules.items():
            if not has("shell", body):
                continue
            with self.subTest(workflow=entry, rule=name, file=path.name):
                self.assertIn("set -euo pipefail", body)
                # Any redirection into the log counts: `exec > {log} 2>&1`,
                # `exec 2> {log}`, `cmd > {log} 2>&1`, `2> {log}`.
                self.assertIn("{log", body, "shell body never redirects into {log}")

    def test_analysis_parameters_have_no_literal_fallback(self):
        for (entry, name), (path, body) in self.rules.items():
            fallbacks = set(FALLBACK_RE.findall(body)) & ANALYSIS_KEYS
            with self.subTest(workflow=entry, rule=name, file=path.name):
                self.assertFalse(
                    fallbacks,
                    f"rule {name} reads {sorted(fallbacks)} with a literal default; "
                    "use config[key] (the CLI owns the default)",
                )

    def test_only_memory_rules_declare_mem_mb(self):
        """mem_mb binds under --max-memory; a figure that was never measured
        would only throttle the run (ResourceDefaults.MEMORY_RULES is the list,
        scoped by each workflow's own resource rule list: Illumina's
        infer_consensus_sequence is samtools consensus, not Clair3)."""
        lists = {
            "consensus_illumina.smk": ResourceDefaults.CONSENSUS_ILLUMINA_RULES,
            "consensus_illumina_segmented.smk": ResourceDefaults.CONSENSUS_ILLUMINA_RULES,
            "consensus_nanopore.smk": ResourceDefaults.CONSENSUS_NANOPORE_RULES,
            "consensus_nanopore_segmented.smk": ResourceDefaults.CONSENSUS_NANOPORE_RULES,
        }
        for (entry, name), (path, body) in self.rules.items():
            declares = "mem_mb" in body
            should = name in ResourceDefaults.MEMORY_RULES and name in lists[entry]
            with self.subTest(workflow=entry, rule=name, file=path.name):
                self.assertEqual(declares, should, "mem_mb declared iff a listed memory rule")
                if declares:
                    self.assertIn(f'ram_mb("{name}")', body)

    def test_rules_read_only_their_own_resource_keys(self):
        for (entry, name), (path, body) in self.rules.items():
            keys = {a or b for a, b in RESOURCE_KEY_RE.findall(body)}
            with self.subTest(workflow=entry, rule=name, file=path.name):
                self.assertTrue(
                    keys <= {name},
                    f"rule {name} reads resource keys of {sorted(keys - {name})}",
                )

    def test_common_threaded_rules_mirror_resource_defaults(self):
        text = (SCRIPTS_DIR / "rules" / "common.smk").read_text()
        for rule in (
            ResourceDefaults.CONSENSUS_ILLUMINA_RULES + ResourceDefaults.CONSENSUS_NANOPORE_RULES
        ):
            self.assertIn(f'"{rule}"', text.split("THREADED_RULES", 1)[1].split("def cpus", 1)[0])

    def test_resource_rule_lists_match_the_workflows(self):
        expectations = {
            "consensus_illumina.smk": ResourceDefaults.CONSENSUS_ILLUMINA_RULES,
            "consensus_illumina_segmented.smk": ResourceDefaults.CONSENSUS_ILLUMINA_RULES,
            "consensus_nanopore.smk": ResourceDefaults.CONSENSUS_NANOPORE_RULES,
            "consensus_nanopore_segmented.smk": ResourceDefaults.CONSENSUS_NANOPORE_RULES,
        }
        for entry in ENTRY_FILES:
            rules = rules_for_entry(entry)
            listed = set(expectations[entry.name])
            reading = {name for name, (_, body) in rules.items() if RESOURCE_KEY_RE.search(body)}
            with self.subTest(workflow=entry.name):
                self.assertTrue(
                    listed <= set(rules),
                    f"ResourceDefaults names rules absent from the DAG: {sorted(listed - set(rules))}",
                )
                self.assertEqual(
                    reading,
                    listed,
                    "rules reading a _cpus/_ram key must be exactly the ResourceDefaults list",
                )


class Test_NanoporeSanitiserCopiesAgree(unittest.TestCase):
    """The contig-name rewrite is written three times: the ``sed`` in both nanopore
    entry workflows, ``integrity._NANOPORE_SANITIZE_RE`` (validation) and
    ``build_report.NANOPORE_UNSAFE`` (the report's BED rewrite). Clair3 passes
    contig names to /bin/sh unquoted, so all three must agree byte for byte."""

    SED_LINE = (
        "LC_ALL=C sed '1s/^\\\\xef\\\\xbb\\\\xbf//; "
        "/^>/{{s/^>[[:space:]]*//; s/[[:space:]]*$//; s/[^A-Za-z0-9._-]/_/g; s/^/>/}}' "
        "{input:q} > {output.fasta:q}"
    )

    def test_both_nanopore_workflows_use_the_exact_sed_line(self):
        for entry in ENTRY_FILES:
            if "nanopore" not in entry.name:
                continue
            with self.subTest(workflow=entry.name):
                self.assertIn(self.SED_LINE, entry.read_text())

    def test_sed_and_python_agree_on_real_headers(self):
        """Run the exact sed over awkward headers and compare with the Python copy."""
        if shutil.which("sed") is None:
            self.skipTest("sed not on PATH")
        headers = [
            "MN908947.3 Severe acute respiratory syndrome coronavirus 2, complete genome",
            "Influenza A virus (A/New York/392/2004(H3N2)) segment 5, complete sequence",
            "hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124",
            "  MN908947.3 x ",
            ">abc\tq",
            "\u00d1:a>b",
            "a\u00a0",
        ]
        # Undo the Snakemake-side escaping: doubled braces and the doubled
        # backslashes of the BOM bytes (Python resolves them before bash).
        program = (
            self.SED_LINE.split("'")[1].replace("{{", "{").replace("}}", "}").replace("\\\\", "\\")
        )
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "ref.fasta"
            body = "\ufeff" + "".join(f">{h}\r\nACGT\r\n" for h in headers)
            fasta.write_text(body, encoding="utf-8", newline="")
            out = subprocess.run(
                ["sed", program, str(fasta)],
                capture_output=True,
                check=True,
                env={"LC_ALL": "C", "PATH": os.environ.get("PATH", "")},
            ).stdout
        sed_names = [line[1:] for line in out.decode("ascii").split("\n") if line.startswith(">")]
        self.assertEqual(sed_names, [integrity.sanitize_nanopore_contig(h) for h in headers])
        self.assertEqual(
            sed_names[0],
            "MN908947.3_Severe_acute_respiratory_syndrome_coronavirus_2__complete_genome",
        )

    def test_python_copies_share_one_pattern(self):
        self.assertEqual(integrity._NANOPORE_SANITIZE_RE.pattern, rb"[^A-Za-z0-9._-]")
        self.assertEqual(
            integrity._NANOPORE_SANITIZE_RE.pattern, build_report.NANOPORE_UNSAFE.pattern
        )


class Test_DryrunConfigsCarryRequiredKeys(unittest.TestCase):
    """Every regression config must carry the keys the rules read as required,
    or the dry-run suite would stop guarding the DAG on that workflow."""

    def test_required_keys_present(self):
        configs = sorted(DRYRUN_CONFIG_DIR.glob("*.yaml"))
        self.assertTrue(configs, "no dry-run configs found")
        for path in configs:
            with open(path) as fh:
                config = yaml.safe_load(fh)
            required = list(REQUIRED_COMMON)
            required += REQUIRED_ILLUMINA if config["data"] == "illumina" else REQUIRED_NANOPORE
            if config.get("run_viralqc", True):
                required.append("run_viralqc_ram")
            missing = [k for k in required if k not in config]
            with self.subTest(config=path.name):
                self.assertFalse(missing, f"missing {missing}")
