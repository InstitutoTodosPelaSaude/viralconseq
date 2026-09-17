"""Machine detection: cores, memory, and the arithmetic that turns memory into
a budget. Every file read goes through one seam, ``constants._read_text``, so
these tests never touch /proc or /sys."""

import os
import unittest
from unittest.mock import patch

from viralconseq import constants
from viralconseq.constants import ResourceDefaults, detect_cores, detect_memory_mb, memory_budget_mb

MEMINFO_32G = "MemTotal:       33554432 kB\nMemFree:        1000000 kB\n"
GIB = 1024**3


def _files(mapping):
    return patch.object(constants, "_read_text", side_effect=lambda path: mapping.get(path))


class Test_DetectCores(unittest.TestCase):
    def test_one_core_is_left_free(self):
        with patch.object(os, "sched_getaffinity", lambda _pid: set(range(32)), create=True):
            self.assertEqual(detect_cores(), 31)

    def test_a_single_core_machine_still_gets_one(self):
        with patch.object(os, "sched_getaffinity", lambda _pid: {0}, create=True):
            self.assertEqual(detect_cores(), 1)

    def test_affinity_preferred_over_host_cpu_count(self):
        with (
            patch.object(os, "sched_getaffinity", lambda _pid: set(range(4)), create=True),
            patch.object(os, "cpu_count", lambda: 64),
        ):
            self.assertEqual(detect_cores(), 3)

    def test_blocked_affinity_call_falls_back(self):
        def blocked(_pid):
            raise OSError("blocked by seccomp")

        with (
            patch.object(os, "sched_getaffinity", blocked, create=True),
            patch.object(os, "cpu_count", lambda: 4),
        ):
            self.assertEqual(detect_cores(), 3)

    def test_unknowable_cpu_count_still_gives_one(self):
        def blocked(_pid):
            raise OSError("blocked")

        with (
            patch.object(os, "sched_getaffinity", blocked, create=True),
            patch.object(os, "cpu_count", lambda: None),
        ):
            self.assertEqual(detect_cores(), 1)


class Test_DetectMemory(unittest.TestCase):
    def test_v2_max_means_no_limit(self):
        with _files(
            {
                "/proc/meminfo": MEMINFO_32G,
                "/proc/self/cgroup": "0::/user.slice\n",
                "/sys/fs/cgroup/user.slice/memory.max": "max\n",
            }
        ):
            self.assertEqual(detect_memory_mb(), 32768)

    def test_v2_numeric_limit_below_memtotal_caps(self):
        with _files(
            {
                "/proc/meminfo": MEMINFO_32G,
                "/proc/self/cgroup": "0::/job\n",
                "/sys/fs/cgroup/job/memory.max": f"{8 * GIB}\n",
            }
        ):
            self.assertEqual(detect_memory_mb(), 8192)

    def test_an_ancestor_limit_is_honoured(self):
        with _files(
            {
                "/proc/meminfo": MEMINFO_32G,
                "/proc/self/cgroup": "0::/a/b\n",
                "/sys/fs/cgroup/a/b/memory.max": "max\n",
                "/sys/fs/cgroup/a/memory.max": f"{4 * GIB}\n",
            }
        ):
            self.assertEqual(detect_memory_mb(), 4096)

    def test_v1_sentinel_means_no_limit(self):
        with _files(
            {
                "/proc/meminfo": MEMINFO_32G,
                "/proc/self/cgroup": "12:pids:/docker/x\n3:memory:/docker/x\n",
                "/sys/fs/cgroup/memory/docker/x/memory.limit_in_bytes": "9223372036854771712\n",
            }
        ):
            self.assertEqual(detect_memory_mb(), 32768)

    def test_v1_numeric_limit_caps(self):
        with _files(
            {
                "/proc/meminfo": MEMINFO_32G,
                "/proc/self/cgroup": "3:memory:/docker/x\n",
                "/sys/fs/cgroup/memory/docker/x/memory.limit_in_bytes": f"{2 * GIB}\n",
            }
        ):
            self.assertEqual(detect_memory_mb(), 2048)

    def test_missing_cgroup_files_fall_back_to_meminfo(self):
        with _files({"/proc/meminfo": MEMINFO_32G, "/proc/self/cgroup": "0::/nowhere\n"}):
            self.assertEqual(detect_memory_mb(), 32768)

    def test_non_linux_uses_sysconf(self):
        conf = {"SC_PAGE_SIZE": 4096, "SC_PHYS_PAGES": 4 * 1024 * 1024}
        with _files({}), patch.object(os, "sysconf", lambda name: conf[name]):
            self.assertEqual(detect_memory_mb(), 16384)

    def test_undetectable_memory_returns_none(self):
        def unavailable(_name):
            raise ValueError("no such sysconf name")

        with _files({}), patch.object(os, "sysconf", unavailable):
            self.assertIsNone(detect_memory_mb())


class Test_MemoryBudget(unittest.TestCase):
    def test_budget_takes_ten_percent_headroom(self):
        self.assertEqual(ResourceDefaults.MEMORY_HEADROOM, 0.10)
        self.assertEqual(memory_budget_mb(32768, 2), (29491, False))

    def test_budget_is_clamped_up_to_the_largest_rule(self):
        self.assertEqual(memory_budget_mb(1024, 2), (2048, True))

    def test_memory_rules_are_the_measured_ones(self):
        self.assertEqual(
            set(ResourceDefaults.MEMORY_RULES), {"infer_consensus_sequence", "run_viralqc"}
        )
        self.assertEqual(ResourceDefaults.ram_for("infer_consensus_sequence"), 2)
        self.assertEqual(ResourceDefaults.ram_for("run_viralqc"), 1)
        self.assertEqual(
            ResourceDefaults.memory_rules_for(ResourceDefaults.CONSENSUS_ILLUMINA_RULES),
            ["run_viralqc"],
        )
        self.assertEqual(
            ResourceDefaults.memory_rules_for(ResourceDefaults.CONSENSUS_NANOPORE_RULES),
            ["infer_consensus_sequence", "run_viralqc"],
        )


if __name__ == "__main__":
    unittest.main()
