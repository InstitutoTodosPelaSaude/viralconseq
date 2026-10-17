"""Tests for scripts/python/normalise_depth.py on synthetic SAM text (no samtools)."""

import importlib.util
import os
import random
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

_SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "viralconseq"
    / "scripts"
    / "python"
    / "normalise_depth.py"
)
_spec = importlib.util.spec_from_file_location("normalise_depth", _SCRIPT)
assert _spec is not None and _spec.loader is not None
nd = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nd)

CONTIG = "ref"


def header(length=1000, contig=CONTIG, pg=True):
    lines = ["@HD\tVN:1.6\tSO:unknown\n", f"@SQ\tSN:{contig}\tLN:{length}\n"]
    if pg:
        lines.append("@PG\tID:minimap2\tPN:minimap2\n")
    return lines


def rec(qname, pos, length, flag=0, cigar=None, rnext="*", pnext=0, contig=CONTIG, mapq=60):
    """One SAM line; ``pos`` is 0-based here, 1-based in the line."""
    cigar = cigar or f"{length}M"
    return f"{qname}\t{flag}\t{contig}\t{pos + 1}\t{mapq}\t{cigar}\t{rnext}\t{pnext + 1 if pnext else 0}\t0\t*\t*\n"


def pair(qname, pos1, pos2, length, contig=CONTIG, pnext_shift=0):
    """A proper pair: forward mate at pos1, reverse mate at pos2. ``pnext_shift``
    makes PNEXT lie, as ampliconclip --hard-clip does."""
    a = rec(qname, pos1, length, flag=99, rnext="=", pnext=pos2 + pnext_shift, contig=contig)
    b = rec(qname, pos2, length, flag=147, rnext="=", pnext=pos1 + pnext_shift, contig=contig)
    return [a, b]


def run(lines, target, scheme=None):
    recs, results, values, rows, hdr = nd.normalise(lines, target, scheme)
    return recs, results, values, rows


def depths(result):
    before = [result.before_fwd[p] + result.before_rev[p] for p in range(result.length)]
    after = [result.after_fwd[p] + result.after_rev[p] for p in range(result.length)]
    return before, after


def sort_records(records):
    return sorted(records, key=lambda line: int(line.split("\t")[3]))


class TestCigar(unittest.TestCase):
    def test_match_deletion_cover_insertion_clips_do_not(self):
        self.assertEqual(nd.cigar_intervals(10, "5M2D3M"), [(10, 20)])
        self.assertEqual(nd.cigar_intervals(10, "5S5M2I5M5H"), [(10, 20)])
        self.assertEqual(nd.cigar_intervals(0, "4=1X5M"), [(0, 10)])

    def test_skip_splits_coverage(self):
        self.assertEqual(nd.cigar_intervals(0, "5M10N5M"), [(0, 5), (15, 20)])

    def test_star_gives_nothing(self):
        self.assertEqual(nd.cigar_intervals(0, "*"), [])


class TestSplitTarget(unittest.TestCase):
    def test_proportional_largest_remainder(self):
        self.assertEqual(nd.split_target(200, [900, 100]), [180, 20])
        self.assertEqual(nd.split_target(100, [50, 50, 50]), [34, 33, 33])

    def test_never_above_part_and_passthrough_when_under_target(self):
        self.assertEqual(nd.split_target(100, [30, 20]), [30, 20])
        self.assertEqual(nd.split_target(0, [5, 5]), [0, 0])

    def test_two_strands_round_half_up_to_forward(self):
        self.assertEqual(nd.split_target(5, [10, 10]), [3, 2])


class TestPositionMode(unittest.TestCase):
    def test_staggered_fixed_length_reads_capped_exactly(self):
        # Reads of 50 at every position: original depth 50 inside, cap 20.
        lines = header(400) + [rec(f"r{p}", p, 50) for p in range(0, 350)]
        recs, results, values, _ = run(lines, 20)
        before, after = depths(results[0])
        # Exactly T away from the coverage boundaries; at the tail the reads
        # kept to replace expiring ones reach back over finished positions
        # (the 2T effect), which the table reports honestly.
        interior = range(60, 300)
        self.assertTrue(all(before[p] == 50 for p in interior))
        self.assertTrue(all(after[p] == 20 for p in interior), set(after[p] for p in interior))
        self.assertTrue(all(after[p] >= min(20, before[p]) for p in range(400)))
        self.assertTrue(all(after[p] <= 40 for p in range(400)))
        self.assertEqual(values["mode"], "position")
        self.assertEqual(values["self_check"], "ok")
        self.assertEqual(int(values["records_kept"]), sum(recs.kept))

    def test_low_coverage_keeps_everything(self):
        lines = header(300) + [rec(f"r{i}", 10 * i, 30) for i in range(20)]
        recs, results, values, _ = run(lines, 100)
        self.assertEqual(sum(recs.kept), 20)
        self.assertEqual(values["records_dropped"], "0")
        self.assertEqual(values["positions_capped"], "0")

    def test_random_reads_hold_lower_bound_and_two_t_per_strand(self):
        rng = random.Random(7)
        records = []
        for i in range(3000):
            pos = rng.randrange(0, 900)
            length = rng.randrange(40, 120)
            flag = 16 if rng.random() < 0.5 else 0
            records.append(rec(f"q{i}", pos, length, flag=flag))
        lines = header(1200) + sort_records(records)
        _, results, _, _ = run(lines, 50)
        res = results[0]
        before, after = depths(res)
        for p in range(1200):
            self.assertGreaterEqual(after[p], min(50, before[p]), p)
            self.assertLessEqual(res.after_fwd[p], 100, p)
            self.assertLessEqual(res.after_rev[p], 100, p)

    def test_strand_ratio_preserved(self):
        records = [rec(f"f{i}", 0, 100) for i in range(900)] + [
            rec(f"r{i}", 0, 100, flag=16) for i in range(100)
        ]
        lines = header(200) + records
        _, results, _, _ = run(lines, 200)
        res = results[0]
        self.assertEqual((res.after_fwd[50], res.after_rev[50]), (180, 20))

    def test_secondary_dropped_supplementary_capped(self):
        records = [rec(f"s{i}", 0, 100) for i in range(30)]
        records += [rec(f"sup{i}", 0, 100, flag=2048) for i in range(30)]
        records += [rec("sec", 0, 100, flag=256)]
        recs, results, values, _ = run(header(200) + records, 20)
        self.assertEqual(values["secondary_dropped"], "1")
        self.assertEqual(depths(results[0])[1][50], 20)
        self.assertEqual(recs.records_in, 61)

    def test_unmapped_records_dropped_and_counted(self):
        records = [rec(f"s{i}", 0, 100) for i in range(5)] + ["u\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"]
        recs, _, values, _ = run(header(200) + records, 20)
        self.assertEqual(values["unmapped_dropped"], "1")
        self.assertEqual(values["records_kept"], "5")

    def test_unsorted_input_refused(self):
        lines = header(200) + [rec("a", 50, 10), rec("b", 10, 10)]
        with self.assertRaises(nd.UnsortedInputError):
            run(lines, 10)

    def test_contig_blocks_must_be_contiguous(self):
        lines = ["@HD\tVN:1.6\n", "@SQ\tSN:A\tLN:100\n", "@SQ\tSN:B\tLN:100\n"]
        lines += [
            rec("a", 0, 10, contig="A"),
            rec("b", 0, 10, contig="B"),
            rec("c", 5, 10, contig="A"),
        ]
        with self.assertRaises(nd.UnsortedInputError):
            run(lines, 10)

    def test_multi_contig(self):
        lines = ["@HD\tVN:1.6\n", "@SQ\tSN:A\tLN:100\n", "@SQ\tSN:B\tLN:50\n"]
        lines += [rec(f"a{i}", 0, 100, contig="A") for i in range(30)]
        lines += [rec(f"b{i}", 0, 50, contig="B") for i in range(5)]
        recs, results, values, rows = run(lines, 10)
        self.assertEqual([r.name for r in results], ["A", "B"])
        self.assertEqual(depths(results[0])[1][10], 10)
        self.assertEqual(depths(results[1])[1][10], 5)
        self.assertEqual([row[0] for row in rows], ["A", "B"])

    def test_skip_operation_does_not_cover(self):
        lines = header(100) + [rec("a", 0, 50, cigar="20M10N20M")]
        _, results, _, _ = run(lines, 10)
        before, after = depths(results[0])
        self.assertEqual(before[25], 0)
        self.assertEqual(before[45], 1)
        self.assertEqual(after[45], 1)


class TestPairs(unittest.TestCase):
    def test_pairs_kept_or_dropped_together_matched_by_qname(self):
        records = []
        for i in range(60):
            records += pair(f"p{i}", 0, 150, 100, pnext_shift=35)  # PNEXT lies by 35 bp
        lines = header(400) + sort_records(records)
        recs, results, values, _ = run(lines, 20)
        self.assertTrue(all(recs.mate[i] >= 0 for i in range(recs.n)))
        self.assertTrue(all(recs.kept[i] == recs.kept[recs.mate[i]] for i in range(recs.n)))
        before, after = depths(results[0])
        self.assertTrue(all(after[p] >= min(20, before[p]) for p in range(400)))
        # the far mate's coverage credits its own positions
        self.assertEqual(after[200], after[50])
        self.assertEqual(values["orphan_mates"], "0")

    def test_orphan_mate_is_a_singleton(self):
        records = [pair(f"p{i}", 0, 150, 100)[0] for i in range(30)]  # forward mates only
        recs, results, values, _ = run(header(400) + records, 10)
        self.assertEqual(values["orphan_mates"], "30")
        self.assertEqual(depths(results[0])[1][50], 10)

    def test_first_mate_kept_through_far_mate_credits_past_positions(self):
        # 40 forward singles cover [0,100); 40 pairs have a forward mate at 0
        # ending at 90 and a reverse mate at 300, the only coverage there.
        # Furthest-first keeps the singles at position 0, so at 300 the pairs'
        # first mates have expired and are kept through the far mate: the
        # finished rows must show that retroactive coverage.
        records = [rec(f"s{i}", 0, 100) for i in range(40)]
        for i in range(40):
            a = rec(f"p{i}", 0, 90, flag=99, rnext="=", pnext=300)
            b = rec(f"p{i}", 300, 100, flag=147, rnext="=", pnext=0)
            records += [a, b]
        recs, results, values, _ = run(header(500) + sort_records(records), 20)
        before, after = depths(results[0])
        self.assertEqual(values["self_check"], "ok")
        self.assertEqual(after[350], 20)
        kept_pairs = sum(1 for i in range(recs.n) if recs.kept[i] and recs.mate[i] >= 0) // 2
        self.assertEqual(kept_pairs, 20)
        self.assertEqual(after[95], 20)  # singles only
        self.assertEqual(after[50], 40)  # singles plus the retroactively credited first mates


class TestEdgeShapes(unittest.TestCase):
    def test_leading_skip_operation_is_ordered_by_pos(self):
        # POS 11 with a leading N covers from 16; POS 13 follows it in the file.
        lines = header(100) + [rec("a", 10, 15, cigar="5N10M"), rec("b", 12, 10)]
        _, results, values, _ = run(lines, 10)
        before, after = depths(results[0])
        self.assertEqual(before[12], 1)  # only b
        self.assertEqual(before[16], 2)
        self.assertEqual(values["records_kept"], "2")

    def test_record_past_contig_length_is_an_error(self):
        lines = header(50) + [rec("a", 40, 20)]
        with self.assertRaises(nd.NormaliseError):
            run(lines, 10)

    def test_overlapping_mates_stay_whole(self):
        records = []
        for i in range(50):
            records += pair(f"p{i}", 0, 20, 100)  # mates overlap over 20..100
        recs, results, _, _ = run(header(300) + sort_records(records), 10)
        before, after = depths(results[0])
        self.assertTrue(all(after[p] >= min(10, before[p]) for p in range(300)))
        self.assertTrue(all(recs.kept[i] == recs.kept[recs.mate[i]] for i in range(recs.n)))

    def test_empty_middle_contig(self):
        lines = ["@HD\tVN:1.6\n", "@SQ\tSN:A\tLN:50\n", "@SQ\tSN:B\tLN:50\n", "@SQ\tSN:C\tLN:50\n"]
        lines += [rec("a", 0, 20, contig="A"), rec("c", 0, 20, contig="C")]
        _, results, values, rows = run(lines, 10)
        self.assertEqual([r.name for r in results], ["A", "B", "C"])
        self.assertEqual(sum(depths(results[1])[0]), 0)
        self.assertEqual(values["records_kept"], "2")


class TestSchemeMode(unittest.TestCase):
    BED = (
        "ref\t10\t30\tscheme_1_LEFT\t1\t+\n"
        "ref\t280\t300\tscheme_1_RIGHT\t1\t-\n"
        "ref\t250\t270\tscheme_2_LEFT\t2\t+\n"
        "ref\t520\t540\tscheme_2_RIGHT\t2\t-\n"
        "ref\t255\t275\tscheme_2_LEFT_alt\t2\t+\n"
        "ref\t490\t510\tscheme_3_LEFT\t1\t+\n"
        "ref\t760\t780\tscheme_3_RIGHT\t1\t-\n"
    )

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.bed = os.path.join(self._tmp.name, "scheme.bed")
        with open(self.bed, "w") as fh:
            fh.write(self.BED)
        self.scheme = nd.parse_scheme(self.bed)

    def tearDown(self):
        self._tmp.cleanup()

    def test_parse_merges_alts_and_orders(self):
        amps = self.scheme.amplicons["ref"]
        self.assertEqual([a.name for a in amps], ["scheme_1", "scheme_2", "scheme_3"])
        a2 = amps[1]
        self.assertEqual((a2.outer_start, a2.outer_end), (250, 540))
        self.assertEqual(
            (a2.insert_start, a2.insert_end), (275, 520)
        )  # the alt LEFT pushes the insert start
        self.assertEqual(a2.pool, "2")
        self.assertEqual(self.scheme.notes, [])

    def test_assignment_classes(self):
        assigner = nd._ContigAssigner(self.scheme.amplicons["ref"], self.scheme.primers["ref"])
        self.assertEqual(assigner.assign(30, 280), (0, nd.CLASS_NORMAL))
        # contains amplicon 1's RIGHT primer but fits amplicon 2: normal in a tiling scheme
        self.assertEqual(assigner.assign(270, 500), (1, nd.CLASS_NORMAL))
        # a skip-amplicon product: fits none, contains whole primers
        self.assertEqual(assigner.assign(30, 530), (-1, nd.CLASS_INTERIOR_PRIMER))
        # fits none, contains no primer
        self.assertEqual(assigner.assign(900, 950), (-1, nd.CLASS_UNASSIGNED))
        # overlap region, fits both: nearest primer ends (amplicon 1 inner 30..280)
        self.assertEqual(assigner.assign(35, 275), (0, nd.CLASS_NORMAL))
        self.assertEqual(assigner.assign(276, 519)[0], 1)
        # 10 bp tolerance beyond the outer span
        self.assertEqual(assigner.assign(5, 285)[1], nd.CLASS_NORMAL)

    def test_overlap_shares_target_by_amplicon_and_may_reach_two_t(self):
        records = [rec(f"a{i}", 30, 250) for i in range(300)]  # amplicon 1: 30..280
        records += [rec(f"b{i}", 270, 250, flag=16) for i in range(100)]  # amplicon 2: 270..520
        lines = header(1000) + sort_records(records)
        recs, results, values, rows = run(lines, 40, self.scheme)
        before, after = depths(results[0])
        self.assertEqual(values["mode"], "amplicon")
        self.assertTrue(all(after[p] >= min(40, before[p]) for p in range(1000)))
        # Amplicon 1 is capped at 40 by position 70; amplicon 2's reads kept
        # for its own region (from 280 on) reach back into the overlap, so the
        # overlap 270..280 shows both amplicons at their caps: 2T, as in ARTIC.
        self.assertEqual((results[0].after_fwd[150], results[0].after_rev[150]), (40, 0))
        self.assertEqual((results[0].after_fwd[400], results[0].after_rev[400]), (0, 40))
        self.assertEqual((results[0].after_fwd[275], results[0].after_rev[275]), (40, 40))
        self.assertTrue(all(after[p] <= 80 for p in range(1000)))
        by_name = {row[1]: row for row in rows}
        self.assertEqual(by_name["scheme_1"][7:9], ["300", "40"])
        self.assertEqual(by_name["scheme_2"][7:9], ["100", "40"])
        self.assertIn("unassigned", by_name)
        self.assertIn("interior_primer", by_name)

    def test_other_class_is_never_dropped_and_counted(self):
        records = [rec(f"a{i}", 30, 250) for i in range(100)]
        records += [rec(f"skip{i}", 30, 500) for i in range(5)]  # skip-amplicon products
        records += [rec(f"far{i}", 900, 50) for i in range(3)]  # fit no amplicon, no primer
        lines = header(1000) + sort_records(records)
        recs, results, values, rows = run(lines, 20, self.scheme)
        self.assertEqual(values["interior_primer_before"], "5")
        self.assertEqual(values["interior_primer_after"], "5")
        self.assertEqual(values["unassigned_before"], "3")
        self.assertEqual(values["unassigned_after"], "3")
        before, after = depths(results[0])
        self.assertTrue(all(after[p] >= min(20, before[p]) for p in range(1000)))
        # the always-kept skip reads count against the target: 5 kept + 15 normal
        self.assertEqual(after[100], 20)
        by_name = {row[1]: row for row in rows}
        self.assertEqual(by_name["interior_primer"][7:11], ["5", "5", "5", "5"])
        self.assertEqual(by_name["unassigned"][7:9], ["3", "3"])

    def test_pair_with_one_skip_amplicon_mate_is_kept_whole(self):
        records = [rec(f"a{i}", 30, 250) for i in range(100)]
        for i in range(5):
            a = rec(f"x{i}", 30, 500, flag=99, rnext="=", pnext=560)  # skip-amplicon forward mate
            b = rec(f"x{i}", 560, 100, flag=147, rnext="=", pnext=30)  # normal mate in amplicon 3
            records += [a, b]
        recs, results, values, _ = run(header(1000) + sort_records(records), 20, self.scheme)
        self.assertEqual(values["interior_primer_after"], "5")
        self.assertTrue(
            all(recs.kept[i] == recs.kept[recs.mate[i]] for i in range(recs.n) if recs.mate[i] >= 0)
        )
        before, after = depths(results[0])
        self.assertTrue(all(after[p] >= min(20, before[p]) for p in range(1000)))

    def test_scheme_fallback_by_strand_column(self):
        bed = os.path.join(self._tmp.name, "odd.bed")
        with open(bed, "w") as fh:
            fh.write(
                "ref\t10\t30\tprimerA\t1\t+\nref\t280\t300\tprimerB\t1\t-\nref\t250\t270\tprimerC\t2\t+\nref\t520\t540\tprimerD\t2\t-\n"
            )
        scheme = nd.parse_scheme(bed)
        self.assertEqual(len(scheme.amplicons["ref"]), 2)
        self.assertTrue(scheme.notes and "strand column" in scheme.notes[0])

    def test_scheme_without_names_or_strand_is_refused(self):
        bed = os.path.join(self._tmp.name, "bad.bed")
        with open(bed, "w") as fh:
            fh.write("ref\t10\t30\tprimerA\nref\t280\t300\tprimerB\n")
        with self.assertRaises(nd.NormaliseError):
            nd.parse_scheme(bed)

    def test_scheme_contig_absent_from_bam_keeps_everything_as_unassigned(self):
        lines = ["@HD\tVN:1.6\n", "@SQ\tSN:other\tLN:300\n"] + [
            rec(f"a{i}", 0, 100, contig="other") for i in range(30)
        ]
        recs, _, values, _ = run(lines, 10, self.scheme)
        self.assertEqual(values["unassigned_before"], "30")
        self.assertEqual(sum(recs.kept), 30)


class TestDeterminismAndEmission(unittest.TestCase):
    def _records(self, seed):
        rng = random.Random(seed)
        records = []
        for i in range(800):
            records.append(
                rec(f"q{i}", rng.randrange(0, 300), rng.randrange(50, 90), flag=16 if i % 3 else 0)
            )
        return records

    def test_same_input_same_kept_set_and_same_output_text(self):
        lines = header(500) + sort_records(self._records(3))
        recs1, *_ = run(lines, 30)
        recs2, *_ = run(lines, 30)
        self.assertEqual(bytes(recs1.kept), bytes(recs2.kept))

    def test_permuting_same_start_records_keeps_the_same_depth_profile(self):
        # The tiebreak among units with the same end is input order, so a
        # permutation may swap which of two identical-shaped reads is kept; the
        # depth profile and the kept count must not move.
        records = sort_records(self._records(5))
        _, results1, values1, _ = run(header(500) + records, 30)
        rng = random.Random(11)
        groups = {}
        for line in records:
            groups.setdefault(line.split("\t")[3], []).append(line)
        permuted = []
        for key in sorted(groups, key=int):
            group = groups[key][:]
            rng.shuffle(group)
            permuted += group
        _, results2, values2, _ = run(header(500) + permuted, 30)
        self.assertEqual(depths(results1[0]), depths(results2[0]))
        self.assertEqual(values1["records_kept"], values2["records_kept"])

    def test_emit_writes_header_with_pg_and_kept_lines_in_order(self):
        lines = header(300) + [rec(f"a{i}", 0, 100) for i in range(5)]
        recs, *_ = run(lines, 2)
        import io

        out = io.StringIO()
        written = nd.emit(lines, recs, out, 2, "position")
        text = out.getvalue().splitlines()
        self.assertEqual(written, 2)
        self.assertEqual(text[0], "@HD\tVN:1.6\tSO:coordinate")
        self.assertTrue(text[2].startswith("@PG\tID:minimap2"))
        self.assertEqual(
            text[3],
            "@PG\tID:normalise_depth\tPN:normalise_depth\tPP:minimap2\tCL:normalise_depth.py --target 2 --mode position",
        )
        self.assertEqual(len(text), 4 + 2)
        self.assertNotIn("/", text[3])
        kept_names = [line.split("\t")[0] for line in text[4:]]
        self.assertEqual(kept_names, sorted(kept_names, key=lambda n: int(n[1:])))

    def test_emit_refuses_a_changed_input(self):
        lines = header(300) + [rec(f"a{i}", 0, 100) for i in range(5)]
        recs, *_ = run(lines, 2)
        import io

        with self.assertRaises(nd.NormaliseError):
            nd.emit(lines[:-1], recs, io.StringIO(), 2, "position")


class TestSelfCheck(unittest.TestCase):
    def test_forced_violation_is_loud(self):
        lines = header(300) + [rec(f"a{i}", 0, 100) for i in range(30)]
        recs, results, *_ = run(lines, 10)
        # drop a kept record behind the sweep's back
        first_kept = next(i for i in range(recs.n) if recs.kept[i])
        recs.kept[first_kept] = 0
        with self.assertRaises(nd.SelfCheckError):
            nd.self_check(recs, results, {CONTIG: 0}, 10)


class TestTables(unittest.TestCase):
    def test_depth_table_rows(self):
        import io

        lines = header(3) + [rec("a", 0, 2), rec("b", 1, 2, flag=16)]
        _, results, *_ = run(lines, 10)
        out = io.StringIO()
        nd.write_depth_table(results, out)
        rows = out.getvalue().splitlines()
        self.assertEqual(rows[0].split("\t"), nd.DEPTH_COLUMNS)
        self.assertEqual(rows[1], "ref\t1\t1\t1\t1\t0\t1\t0")
        self.assertEqual(rows[2], "ref\t2\t2\t2\t1\t1\t1\t1")
        self.assertEqual(rows[3], "ref\t3\t1\t1\t0\t1\t0\t1")

    def test_summary_keys_and_placeholders(self):
        with tempfile.TemporaryDirectory() as tmp:
            d, a, s = (os.path.join(tmp, n) for n in ("d.tsv", "a.tsv", "s.txt"))
            nd.write_placeholders(d, a, s)
            self.assertEqual(open(d).read().rstrip("\n").split("\t"), nd.DEPTH_COLUMNS)
            self.assertEqual(open(a).read().rstrip("\n").split("\t"), nd.AMPLICON_COLUMNS)
            values = dict(line.split("\t") for line in open(s).read().splitlines())
            self.assertEqual(list(values), nd.SUMMARY_KEYS)
            self.assertEqual(
                (values["mode"], values["target"], values["self_check"]), ("off", "0", "NA")
            )


class TestMainOff(unittest.TestCase):
    def test_target_zero_copies_without_opening_the_bam(self):
        with tempfile.TemporaryDirectory() as tmp:
            bam = os.path.join(tmp, "in.bam")
            with open(bam, "wb") as fh:
                fh.write(b"not a bam")
            with open(bam + ".bai", "wb") as fh:
                fh.write(b"idx")
            out = os.path.join(tmp, "sub", "out.bam")
            argv = [
                "--bam",
                bam,
                "--output-bam",
                out,
                "--target",
                "0",
                "--depth-tsv",
                os.path.join(tmp, "d.tsv"),
                "--amplicons-tsv",
                os.path.join(tmp, "a.tsv"),
                "--summary",
                os.path.join(tmp, "s.txt"),
            ]
            with mock.patch.object(nd.subprocess, "Popen") as popen:
                self.assertEqual(nd.main(argv), 0)
            popen.assert_not_called()
            self.assertEqual(open(out, "rb").read(), b"not a bam")
            self.assertEqual(open(out + ".bai", "rb").read(), b"idx")
            self.assertIn("mode\toff", open(os.path.join(tmp, "s.txt")).read())

    def test_negative_target_rejected(self):
        argv = ["--bam", "x", "--output-bam", "y", "--target", "-1"]
        argv += ["--depth-tsv", "d", "--amplicons-tsv", "a", "--summary", "s"]
        self.assertEqual(nd.main(argv), 2)


def _main_argv(tmp, target="10", scheme=None):
    argv = ["--bam", os.path.join(tmp, "in.bam"), "--output-bam", os.path.join(tmp, "out", "n.bam")]
    argv += ["--target", target, "--depth-tsv", os.path.join(tmp, "d.tsv")]
    argv += ["--amplicons-tsv", os.path.join(tmp, "a.tsv"), "--summary", os.path.join(tmp, "s.txt")]
    if scheme:
        argv += ["--scheme", scheme]
    return argv


class TestMainWithoutSamtools(unittest.TestCase):
    """main() with the samtools plumbing patched to serve SAM text."""

    def _run(self, lines, argv):
        with mock.patch.object(nd, "_samtools_view", return_value=None):
            with mock.patch.object(nd, "_iter_lines", side_effect=lambda proc: iter(lines)):
                with mock.patch.object(nd.subprocess, "Popen") as popen:
                    rc = nd.main(argv)
        return rc, popen

    def test_unsorted_input_exits_three_and_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            lines = header(200) + [rec("a", 50, 10), rec("b", 10, 10)]
            rc, popen = self._run(lines, _main_argv(tmp))
            self.assertEqual(rc, nd.EXIT_UNSORTED)
            popen.assert_not_called()
            self.assertFalse(os.path.exists(os.path.join(tmp, "d.tsv")))
            self.assertFalse(os.path.exists(os.path.join(tmp, "s.txt")))

    def test_self_check_failure_exits_four_and_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            lines = header(200) + [rec(f"a{i}", 0, 100) for i in range(5)]
            with mock.patch.object(nd, "self_check", side_effect=nd.SelfCheckError("forced")):
                rc, popen = self._run(lines, _main_argv(tmp))
            self.assertEqual(rc, nd.EXIT_SELF_CHECK)
            popen.assert_not_called()
            self.assertFalse(os.path.exists(os.path.join(tmp, "d.tsv")))

    def test_header_without_sq_is_a_clean_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            rc, popen = self._run(["@HD\tVN:1.6\n"], _main_argv(tmp))
            self.assertEqual(rc, 1)
            popen.assert_not_called()


@unittest.skipUnless(shutil.which("samtools"), "samtools not on PATH")
class TestMainWithSamtools(unittest.TestCase):
    """End to end through real samtools: BAM in, indexed BAM and tables out."""

    def test_round_trip(self):
        import subprocess

        with tempfile.TemporaryDirectory() as tmp:
            records = [rec(f"a{i}", 0, 100) for i in range(30)] + [
                rec(f"b{i}", 50, 100, flag=16) for i in range(30)
            ]
            sam = os.path.join(tmp, "in.sam")
            with open(sam, "w") as fh:
                fh.writelines(header(300) + records)
            bam = os.path.join(tmp, "in.bam")
            subprocess.run(["samtools", "view", "-b", "-o", bam, sam], check=True)
            subprocess.run(["samtools", "index", bam], check=True)
            self.assertEqual(nd.main(_main_argv(tmp, target="10")), 0)
            out = os.path.join(tmp, "out", "n.bam")
            subprocess.run(["samtools", "index", out], check=True)
            kept = subprocess.run(
                ["samtools", "view", "-c", out], capture_output=True, text=True, check=True
            )
            values = dict(
                line.split("\t") for line in open(os.path.join(tmp, "s.txt")).read().splitlines()
            )
            self.assertEqual(int(kept.stdout), int(values["records_kept"]))
            recs, *_ = run(header(300) + records, 10)
            self.assertEqual(int(values["records_kept"]), sum(recs.kept))
            depth = subprocess.run(
                ["samtools", "depth", "-a", "-J", out], capture_output=True, text=True, check=True
            )
            table = {
                int(f[1]): int(f[3])
                for f in (
                    row.split("\t")
                    for row in open(os.path.join(tmp, "d.tsv")).read().splitlines()[1:]
                )
            }
            for line in depth.stdout.splitlines():
                _c, pos, d = line.split("\t")
                self.assertEqual(int(d), table[int(pos)], pos)
            hdr = subprocess.run(
                ["samtools", "view", "-H", out], capture_output=True, text=True, check=True
            ).stdout
            self.assertIn("SO:coordinate", hdr)
            self.assertIn("ID:normalise_depth", hdr)


if __name__ == "__main__":
    unittest.main()
