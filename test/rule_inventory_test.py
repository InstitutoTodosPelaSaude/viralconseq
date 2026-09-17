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
