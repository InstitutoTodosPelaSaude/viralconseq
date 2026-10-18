"""Tests for the seqkit-backed FASTQ detector (viralconseq.seqkit)."""

import gzip
import os
import random
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

from viralconseq import seqkit

# The seqkit on PATH runs the integration class (CI installs it); without
# one the class is skipped.
_RED_ERRO = "\x1b[31m[ERRO]\x1b[0m "
_WARN_LINE = "\x1b[33m[WARN]\x1b[0m you may switch on flag -g/--remove-gaps to remove spaces\n"


def _completed(returncode, stdout="", stderr=""):
    return subprocess.CompletedProcess(
        args=["seqkit"], returncode=returncode, stdout=stdout, stderr=stderr
    )


class _TmpBase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def _write(self, name, data):
        path = os.path.join(self.tmp, name)
        with open(path, "wb") as fh:
            fh.write(data)
        return path

    def _write_gz(self, name, data):
        path = os.path.join(self.tmp, name)
        with gzip.open(path, "wb") as fh:
            fh.write(data)
        return path


class TestStripAnsi(unittest.TestCase):
    def test_removes_colour_codes(self):
        self.assertEqual(seqkit.strip_ansi(_RED_ERRO + "boom"), "[ERRO] boom")

    def test_plain_text_unchanged(self):
        self.assertEqual(seqkit.strip_ansi("[ERRO] boom"), "[ERRO] boom")


class TestFindSeqkit(unittest.TestCase):
    def test_missing_binary_raises(self):
        with mock.patch("viralconseq.seqkit.shutil.which", return_value=None):
            with self.assertRaises(seqkit.SeqkitMissingError) as ctx:
                seqkit.find_seqkit()
        msg = str(ctx.exception)
        self.assertIn("seqkit", msg)
        self.assertIn("--skip-input-validation", msg)

    def test_found_binary_is_absolute(self):
        with mock.patch("viralconseq.seqkit.shutil.which", return_value="bin/seqkit"):
            self.assertTrue(os.path.isabs(seqkit.find_seqkit()))

    def test_module_imports_without_seqkit(self):
        # Nothing at import time touches the binary; only find_seqkit/detect do.
        self.assertEqual(seqkit.SEQKIT_ARGS[0], "seq")
        self.assertIn("--max-len", seqkit.SEQKIT_ARGS)


class TestProbeFastq(_TmpBase):
    def test_zero_byte_file_is_empty(self):
        p = self._write("empty.fastq", b"")
        res = seqkit.probe_fastq(p)
        self.assertEqual(res.verdict, "empty")
        self.assertFalse(res.ran_seqkit)

    def test_empty_gzip_is_empty(self):
        p = self._write_gz("empty.fastq.gz", b"")
        self.assertEqual(seqkit.probe_fastq(p).verdict, "empty")

    def test_missing_path_is_unreadable(self):
        res = seqkit.probe_fastq(os.path.join(self.tmp, "nope.fastq"))
        self.assertEqual(res.verdict, "unreadable")
        self.assertIn("nope.fastq", res.message)

    def test_directory_is_unreadable(self):
        res = seqkit.probe_fastq(self.tmp)
        self.assertEqual(res.verdict, "unreadable")
        self.assertTrue(res.message)

    def test_gzip_magic_garbage_is_unreadable(self):
        p = self._write("garbage.fastq.gz", b"\x1f\x8b" + b"not a gzip stream" * 4)
        res = seqkit.probe_fastq(p)
        self.assertEqual(res.verdict, "unreadable")
        self.assertTrue(res.message)

    def test_plain_fastq_defers_to_seqkit(self):
        p = self._write("ok.fastq", b"@r1\nACGT\n+\nIIII\n")
        self.assertIsNone(seqkit.probe_fastq(p))

    def test_gzip_fastq_defers_to_seqkit(self):
        p = self._write_gz("ok.fastq.gz", b"@r1\nACGT\n+\nIIII\n")
        self.assertIsNone(seqkit.probe_fastq(p))

    def test_single_byte_file_defers_to_seqkit(self):
        p = self._write("one.fastq", b"@")
        self.assertIsNone(seqkit.probe_fastq(p))

    def test_fasta_given_as_fastq_is_suspect(self):
        # seqkit reads a FASTA happily; the Python scanner must name it.
        for name, writer in (("a.fastq", self._write), ("a.fastq.gz", self._write_gz)):
            with self.subTest(fixture=name):
                res = seqkit.probe_fastq(writer(name, b">r1\nACGT\n>r2\nACGTN\n"))
                self.assertEqual(res.verdict, "suspect")
                self.assertIn("FASTA", res.message)
                self.assertFalse(res.ran_seqkit)

    def test_wrapped_fastq_is_suspect(self):
        # fastp emits zero reads for a multi-line record; seqkit accepts it.
        wrapped = b"@r1\nACGTACGTAC\nGTACGTACGT\n+\nIIIIIIIIII\nIIIIIIIIII\n"
        for name, writer in (("w.fastq", self._write), ("w.fastq.gz", self._write_gz)):
            with self.subTest(name=name):
                res = seqkit.probe_fastq(writer(name, wrapped))
                self.assertIsNotNone(res)
                self.assertEqual(res.verdict, "suspect")
                self.assertIn("line-wrapped", res.message)

    def test_wrapped_fastq_after_bom_and_blank_lines_is_suspect(self):
        p = self._write("w.fastq", b"\xef\xbb\xbf\n\n@r1\nACGT\nACGT\n+\nIIIIIIII\n")
        self.assertEqual(seqkit.probe_fastq(p).verdict, "suspect")

    def test_crlf_and_indented_separator_defer_to_seqkit(self):
        self.assertIsNone(
            seqkit.probe_fastq(self._write("c.fastq", b"@r1\r\nACGT\r\n+\r\nIIII\r\n"))
        )
        self.assertIsNone(seqkit.probe_fastq(self._write("p.fastq", b"@r1\nACGT\n+r1\nIIII\n")))

    def test_short_or_giant_first_record_defers_to_seqkit(self):
        # Fewer than three complete lines: a truncated file is seqkit's call.
        self.assertIsNone(seqkit.probe_fastq(self._write("t.fastq", b"@r1\nACGT\n+")))
        self.assertIsNone(seqkit.probe_fastq(self._write("t2.fastq", b"@r1\nACGT\n")))
        # A sequence line longer than the probe window is not a wrapped record.
        giant = b"@r1\n" + b"A" * (seqkit._PROBE_LIMIT + 10) + b"\n+\n" + b"I" * 10 + b"\n"
        self.assertIsNone(seqkit.probe_fastq(self._write("g.fastq", giant)))

    def test_other_compression_is_unreadable(self):
        import bz2
        import lzma

        cases = {
            "r.fastq.bz2": (bz2.compress(b"@r1\nACGT\n+\nIIII\n"), "bzip2"),
            "r.fastq.xz": (lzma.compress(b"@r1\nACGT\n+\nIIII\n"), "xz"),
            "r.fastq.zst": (b"\x28\xb5\x2f\xfd" + b"\x00" * 16, "zstd"),
        }
        for name, (blob, label) in cases.items():
            with self.subTest(fixture=name):
                res = seqkit.probe_fastq(self._write(name, blob))
                self.assertEqual(res.verdict, "unreadable")
                self.assertIn(label, res.message)

    def test_whitespace_or_bom_only_is_empty(self):
        cases = {
            "ws.fastq.gz": self._write_gz("ws.fastq.gz", b"\n\n\n"),
            "bom.fastq": self._write("bom.fastq", b"\xef\xbb\xbf"),
            "bom_ws.fastq": self._write("bom_ws.fastq", b"\xef\xbb\xbf\r\n  \n"),
        }
        for name, path in cases.items():
            with self.subTest(fixture=name):
                self.assertEqual(seqkit.probe_fastq(path).verdict, "empty")

    def test_bom_then_records_defers_to_seqkit(self):
        p = self._write("bom_ok.fastq", b"\xef\xbb\xbf@r1\nACGT\n+\nIIII\n")
        self.assertIsNone(seqkit.probe_fastq(p))

    def test_long_whitespace_prefix_defers_to_seqkit(self):
        p = self._write("long_ws.fastq", b"\n" * (2 << 20) + b"@r1\nACGT\n+\nIIII\n")
        self.assertIsNone(seqkit.probe_fastq(p))


class TestDetectMapping(_TmpBase):
    """detect() with subprocess.run patched: exercises the verdict mapping."""

    def setUp(self):
        super().setUp()
        self.path = self._write("in.fastq", b"@r1\nACGT\n+\nIIII\n")

    def _detect(self, completed):
        with mock.patch("viralconseq.seqkit.subprocess.run", return_value=completed) as run:
            res = seqkit.detect(self.path, seqkit="/opt/bin/seqkit")
        return res, run

    def test_clean(self):
        res, run = self._detect(_completed(0, stdout="", stderr=_WARN_LINE))
        self.assertEqual(res.verdict, "clean")
        self.assertEqual(res.returncode, 0)
        self.assertEqual(res.message, "")
        self.assertTrue(res.ran_seqkit)
        self.assertEqual(res.command, ["/opt/bin/seqkit", *seqkit.SEQKIT_ARGS, "--", self.path])
        self.assertEqual(run.call_args.args[0], res.command)
        self.assertIs(run.call_args.kwargs["stdin"], subprocess.DEVNULL)

    def test_bare_newline_on_stdout_is_suspect(self):
        # A file cut right after "@": a zero-length read with an empty id.
        res, _ = self._detect(_completed(0, stdout="\n"))
        self.assertEqual(res.verdict, "suspect")

    def test_suspect_via_zero_length_read_on_stdout(self):
        res, _ = self._detect(_completed(0, stdout="h2\n", stderr=_WARN_LINE))
        self.assertEqual(res.verdict, "suspect")
        self.assertEqual(res.message, "")

    def test_suspect_via_bad_fastq_format(self):
        stderr = _WARN_LINE + _RED_ERRO + "fastx: bad fastq format\n"
        res, _ = self._detect(_completed(255, stderr=stderr))
        self.assertEqual(res.verdict, "suspect")
        self.assertEqual(res.returncode, 255)
        self.assertEqual(res.message, "[ERRO] fastx: bad fastq format")

    def test_suspect_via_unequal_sequence_and_quality(self):
        stderr = _RED_ERRO + "fastx: unequal sequence and quality\n"
        res, _ = self._detect(_completed(255, stderr=stderr))
        self.assertEqual(res.verdict, "suspect")

    def test_unreadable_via_unexpected_eof(self):
        stderr = _WARN_LINE + _RED_ERRO + "unexpected EOF\n"
        res, _ = self._detect(_completed(255, stderr=stderr))
        self.assertEqual(res.verdict, "unreadable")
        self.assertEqual(res.message, "[ERRO] unexpected EOF")

    def test_unreadable_via_invalid_gzip_header(self):
        stderr = _RED_ERRO + "fastx: gzip: invalid header\n"
        res, _ = self._detect(_completed(255, stderr=stderr))
        self.assertEqual(res.verdict, "unreadable")

    def test_nonzero_without_erro_line_is_suspect(self):
        res, _ = self._detect(_completed(1, stderr="something odd\n"))
        self.assertEqual(res.verdict, "suspect")
        self.assertEqual(res.message, "")

    def test_default_binary_comes_from_find_seqkit(self):
        with mock.patch("viralconseq.seqkit.find_seqkit", return_value="/x/seqkit") as find:
            with mock.patch("viralconseq.seqkit.subprocess.run", return_value=_completed(0)):
                res = seqkit.detect(self.path)
        find.assert_called_once_with()
        self.assertEqual(res.command[0], "/x/seqkit")

    def test_explicit_missing_binary_raises(self):
        with mock.patch("viralconseq.seqkit.subprocess.run", side_effect=FileNotFoundError()):
            with self.assertRaises(seqkit.SeqkitMissingError):
                seqkit.detect(self.path, seqkit="/nowhere/seqkit")

    def test_probe_short_circuits_before_seqkit(self):
        empty = self._write("empty.fastq", b"")
        with mock.patch("viralconseq.seqkit.subprocess.run") as run:
            res = seqkit.detect(empty, seqkit="/opt/bin/seqkit")
        run.assert_not_called()
        self.assertEqual(res.verdict, "empty")


def _real_seqkit():
    return shutil.which("seqkit")


_REAL_SEQKIT = _real_seqkit()

_RECORDS = b"@r1\nACGT\n+\nIIII\n@r2\nACGTN\n+\nIIIII\n"


def _truncated_gzip_blob(n_reads=500):
    """A gzip of ``n_reads`` 100-mers cut in half: the first byte decompresses
    (so the probe defers to seqkit) but the stream has no trailer, which
    seqkit reports as ``unexpected EOF``."""
    rng = random.Random(0)
    reads = []
    for i in range(n_reads):
        seq = "".join(rng.choice("ACGT") for _ in range(100))
        reads.append(f"@read{i}\n{seq}\n+\n{'I' * 100}\n")
    blob = gzip.compress("".join(reads).encode())
    return blob[: len(blob) // 2]


@unittest.skipUnless(_REAL_SEQKIT, "seqkit not on PATH (install it from bioconda)")
class TestDetectIntegration(unittest.TestCase):
    """Runs the real seqkit over fixtures built from a small synthetic FASTQ."""

    EXPECTED = {
        "clean": [
            "clean.fastq",
            "clean.fastq.gz",
            "crlf.fastq",
            "bom.fastq",
            "blanklines.fastq",
            "multimember.fastq.gz",
            "plain_named_gz.fastq.gz",
            "gz_named_plain.fastq",
        ],
        "empty": ["empty0.fastq", "emptygz.fastq.gz"],
        "unreadable": ["trunc_gzip.fastq.gz", "garbage.fastq.gz"],
        "suspect": [
            "trunc_after_seq.fastq",
            "trunc_after_seq.fastq.gz",
            "trunc_after_sep.fastq",
            "trunc_after_header.fastq",
            "trunc_bytes.fastq",
            "lenmismatch.fastq",
            "badseq.fastq",
            "badheader.fastq",
            "badsep.fastq",
        ],
    }

    # Documented seqkit blind spots: a space in the quality line, and
    # whitespace inside a sequence line when the lengths still match.
    BLIND_SPOTS = ["badqual_space.fastq", "alpha_space.fastq"]

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.fixtures = cls._tmp.name
        plain = {
            "clean.fastq": _RECORDS,
            "crlf.fastq": _RECORDS.replace(b"\n", b"\r\n"),
            "bom.fastq": b"\xef\xbb\xbf" + _RECORDS,
            "blanklines.fastq": b"@r1\nACGT\n+\nIIII\n\n@r2\nACGT\n+\nIIII\n",
            "multimember.fastq.gz": gzip.compress(_RECORDS) + gzip.compress(_RECORDS),
            "plain_named_gz.fastq.gz": _RECORDS,
            "empty0.fastq": b"",
            "trunc_gzip.fastq.gz": _truncated_gzip_blob(),
            "garbage.fastq.gz": b"\x1f\x8b" + b"not a gzip stream" * 4,
            "trunc_after_seq.fastq": _RECORDS + b"@r3\nACGT\n",
            "trunc_after_sep.fastq": _RECORDS + b"@r3\nACGT\n+\n",
            "trunc_after_header.fastq": _RECORDS + b"@r3\n",
            "trunc_bytes.fastq": _RECORDS + b"@r3\nAC",
            "lenmismatch.fastq": b"@r1\nACGT\n+\nIII\n",
            "badseq.fastq": b"@r1\nACGZ\n+\nIIII\n",
            "badheader.fastq": b">r1\nACGT\n+\nIIII\n",
            "badsep.fastq": b"@r1\nACGT\n*\nIIII\n",
            "badqual_space.fastq": b"@r1\nACGT\n+\nII I\n",
            "alpha_space.fastq": b"@r1\nAC GT\n+\nIIIII\n",
        }
        gzipped = {
            "clean.fastq.gz": _RECORDS,
            "gz_named_plain.fastq": _RECORDS,
            "emptygz.fastq.gz": b"",
            "trunc_after_seq.fastq.gz": _RECORDS + b"@r3\nACGT\n",
        }
        for name, data in plain.items():
            with open(os.path.join(cls.fixtures, name), "wb") as fh:
                fh.write(data)
        for name, data in gzipped.items():
            with gzip.open(os.path.join(cls.fixtures, name), "wb") as fh:
                fh.write(data)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def _detect(self, name):
        return seqkit.detect(os.path.join(self.fixtures, name), seqkit=_REAL_SEQKIT)

    def _check(self, verdict):
        for name in self.EXPECTED[verdict]:
            with self.subTest(fixture=name):
                res = self._detect(name)
                self.assertEqual(res.verdict, verdict, res)

    def test_clean(self):
        self._check("clean")

    def test_empty(self):
        self._check("empty")

    def test_unreadable(self):
        self._check("unreadable")

    def test_suspect(self):
        self._check("suspect")

    def test_clean_result_fields(self):
        res = self._detect("clean.fastq.gz")
        self.assertEqual(res.returncode, 0)
        self.assertEqual(res.message, "")
        self.assertTrue(res.ran_seqkit)
        self.assertEqual(res.command[0], _REAL_SEQKIT)
        self.assertEqual(res.command[1:-2], seqkit.SEQKIT_ARGS)
        self.assertEqual(res.command[-2], "--")

    def test_empty_verdicts_come_from_the_probe(self):
        for name in self.EXPECTED["empty"]:
            with self.subTest(fixture=name):
                self.assertFalse(self._detect(name).ran_seqkit)

    def test_truncated_gzip_is_seqkit_unexpected_eof(self):
        res = self._detect("trunc_gzip.fastq.gz")
        self.assertTrue(res.ran_seqkit)
        self.assertEqual(res.returncode, 255)
        self.assertEqual(res.message, "[ERRO] unexpected EOF")

    def test_zero_length_read_is_suspect_with_exit_zero(self):
        # A cut after a header line yields a zero-length read that --max-len 0
        # prints: exit 0 with a read id on stdout.
        res = self._detect("trunc_after_header.fastq")
        self.assertEqual(res.returncode, 0)
        self.assertEqual(res.message, "")
        self.assertEqual(res.verdict, "suspect")

    def test_error_message_is_ansi_stripped(self):
        res = self._detect("badseq.fastq")
        self.assertEqual(res.returncode, 255)
        self.assertNotIn("\x1b", res.message)
        self.assertTrue(res.message.startswith("[ERRO]"), res.message)
        self.assertIn("invalid DNAredundant letter", res.message)

    def test_blind_spots_are_clean(self):
        # Pinned so a future seqkit that starts catching them is noticed (the
        # Python scanner would then be redundant for those codes).
        for name in self.BLIND_SPOTS:
            with self.subTest(fixture=name):
                self.assertEqual(self._detect(name).verdict, "clean")


if __name__ == "__main__":
    unittest.main()
