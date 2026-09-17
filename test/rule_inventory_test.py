"""Static conventions over the Snakemake rule files.

Pure text checks, no Snakemake import: the four entry workflows under
``viralconseq/scripts/`` and everything they ``include:`` are scanned rule by
rule. The conventions enforced here are the ones that, when broken, hide
failures rather than cause them (a rule without a log, a shell body that never
writes to its log) or move resources the operator never named (a rule reading
another rule's ``_cpus`` key).
"""

import re
import unittest
from pathlib import Path
from typing import Dict, Iterator, List, Tuple

import yaml

import viralconseq
from viralconseq.constants import ResourceDefaults

SCRIPTS_DIR = Path(viralconseq.__file__).resolve().parent / "scripts"
ENTRY_FILES = sorted(SCRIPTS_DIR.glob("consensus_*.smk"))

RULE_RE = re.compile(r"^rule (\w+):", re.MULTILINE)
INCLUDE_RE = re.compile(r'^include:\s*"([^"]+\.smk)"', re.MULTILINE)
# A rule body ends at the next line that starts in column 0.
TOP_LEVEL_RE = re.compile(r"^(?=\S)", re.MULTILINE)
RESOURCE_KEY_RE = re.compile(r'config\.get\("(\w+)_(?:cpus|ram)"')
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
    "minimum_length",
    "af_threshold",
]
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
    "chunk_size",
    "clair3_model",
    "variant_quality",
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

    def test_every_rule_declares_log_and_benchmark(self):
        for (entry, name), (path, body) in self.rules.items():
            if name in LOG_EXEMPT:
                continue
            with self.subTest(workflow=entry, rule=name, file=path.name):
                self.assertTrue(has("log", body), "missing log:")
                self.assertTrue(has("benchmark", body), "missing benchmark:")

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

    def test_rules_read_only_their_own_resource_keys(self):
        for (entry, name), (path, body) in self.rules.items():
            keys = set(RESOURCE_KEY_RE.findall(body))
            with self.subTest(workflow=entry, rule=name, file=path.name):
                self.assertTrue(
                    keys <= {name},
                    f"rule {name} reads resource keys of {sorted(keys - {name})}",
                )

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
            missing = [k for k in required if k not in config]
            with self.subTest(config=path.name):
                self.assertFalse(missing, f"missing {missing}")
