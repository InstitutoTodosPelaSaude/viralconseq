"""Tests for scripts/python/collect_benchmarks.py, driven as the workflow does
(``main(argv)``) over a fake ``logs/`` tree."""

import os
import tempfile
import unittest

from viralconseq.scripts.python import collect_benchmarks

HEADER = "s\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time\n"
ROW = "1.5\t0:00:01\t10.0\t20.0\t5.0\t6.0\t0.0\t0.0\t0.0\t0.5\n"


def _write(root, rel, text):
    path = os.path.join(root, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


def _read(path):
    with open(path) as fh:
        return [line.rstrip("\n").split("\t") for line in fh]


class Test_CollectBenchmarks(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.logs = os.path.join(self._tmp.name, "logs")
        self.out = os.path.join(self._tmp.name, "benchmark.tsv")

    def tearDown(self):
        self._tmp.cleanup()

    def _run(self, *extra):
        rc = collect_benchmarks.main(
            ["--logs-dir", self.logs, "--output", self.out, "--samples", "sample-b", "sample-a"]
            + list(extra)
        )
        self.assertEqual(rc, 0)
        return _read(self.out)

    def test_single_reference_layout(self):
        _write(self.logs, "map_reads/sample-a.benchmark.txt", HEADER + ROW)
        _write(self.logs, "map_reads/sample-b.benchmark.txt", HEADER + ROW)
        _write(self.logs, "run_viralqc/run_viralqc.benchmark.txt", HEADER + ROW)
        rows = self._run("--cpus", "map_reads=4", "run_viralqc=2")
        self.assertEqual(rows[0][:5], ["sample", "rule", "target", "threads", "s"])
        self.assertEqual(rows[0][-1], "cpu_time")
        # configured sample order first, run-level rows last
        self.assertEqual([r[0] for r in rows[1:]], ["sample-b", "sample-a", "All"])
        self.assertEqual(rows[1][1:4], ["map_reads", "sample-b", "4"])
        self.assertEqual(rows[3][1:4], ["run_viralqc", "run_viralqc", "2"])
        self.assertEqual(rows[1][4], "1.5")

    def test_segmented_layout_adds_segment_column(self):
        _write(self.logs, "map_reads/S/sample-a.benchmark.txt", HEADER + ROW)
        _write(self.logs, "map_reads/L/sample-a.benchmark.txt", HEADER + ROW)
        _write(self.logs, "versions/versions.benchmark.txt", HEADER + ROW)
        _write(
            self.logs,
            "align_consensus_to_reference_genome/S/align_consensus_to_reference_genome.benchmark.txt",
            HEADER + ROW,
        )
        rows = self._run("--segmented")
        self.assertEqual(rows[0][:5], ["sample", "segment", "rule", "target", "threads"])
        by_key = {(r[0], r[1], r[2]): r for r in rows[1:]}
        self.assertIn(("sample-a", "L", "map_reads"), by_key)
        self.assertIn(("sample-a", "S", "map_reads"), by_key)
        self.assertEqual(by_key[("All", "-", "versions")][3], "versions")
        self.assertEqual(by_key[("All", "S", "align_consensus_to_reference_genome")][4], "")

    def test_header_only_and_excluded_files_are_skipped(self):
        _write(self.logs, "map_reads/sample-a.benchmark.txt", HEADER + ROW)
        _write(self.logs, "killed_rule/sample-a.benchmark.txt", HEADER)
        own = _write(self.logs, "collect_benchmarks/collect_benchmarks.benchmark.txt", HEADER + ROW)
        rows = self._run("--exclude", own)
        self.assertEqual([r[1] for r in rows[1:]], ["map_reads"])

    def test_unexpected_layout_is_skipped_not_fatal(self):
        _write(self.logs, "deep/er/than/expected.benchmark.txt", HEADER + ROW)
        _write(self.logs, "map_reads/sample-a.benchmark.txt", HEADER + ROW)
        rows = self._run()
        self.assertEqual(len(rows), 2)

    def test_columns_follow_the_files_header(self):
        _write(
            self.logs,
            "map_reads/sample-a.benchmark.txt",
            "s\th:m:s\tmax_rss\tnew_metric\n1\t0:00:01\t2\t3\n",
        )
        rows = self._run()
        self.assertEqual(rows[0][4:], ["s", "h:m:s", "max_rss", "new_metric"])
        self.assertEqual(rows[1][-1], "3")

    def test_bad_cpus_entry_exits(self):
        with self.assertRaises(SystemExit):
            collect_benchmarks.main(
                ["--logs-dir", self.logs, "--output", self.out, "--cpus", "nonsense"]
            )


if __name__ == "__main__":
    unittest.main()
