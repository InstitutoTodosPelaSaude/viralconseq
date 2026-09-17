"""Tests for scripts/python/fetch_clair3_model.py with the network replaced by
an in-memory opener."""

import io
import os
import tempfile
import unittest
import urllib.error
import zipfile

from viralconseq.scripts.python import fetch_clair3_model as fetch


def _valid_checkpoint_bytes(payload=b"0"):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("model/data.pkl", b"\x80\x02N.")
        archive.writestr("model/data/0", payload)
    return buffer.getvalue()


class _Opener:
    """Maps URL -> bytes (or an exception to raise); records every request."""

    def __init__(self, served):
        self.served = served
        self.requests = []

    def __call__(self, url, timeout=None):
        self.requests.append(url)
        if url not in self.served:
            raise urllib.error.URLError(f"404 {url}")
        value = self.served[url]
        if isinstance(value, Exception):
            raise value
        return io.BytesIO(value)


PRIMARY = "https://primary/models/m1"
BACKUP = "https://backup/m1"


def _serve(base, pileup=None, full=None):
    return {
        f"{base}/pileup.pt": pileup if pileup is not None else _valid_checkpoint_bytes(b"p"),
        f"{base}/full_alignment.pt": full if full is not None else _valid_checkpoint_bytes(b"f"),
    }


class Test_FetchModel(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.model_dir = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def _listing(self):
        return sorted(os.listdir(self.model_dir))

    def test_success_publishes_both_checkpoints_and_source(self):
        opener = _Opener(_serve(PRIMARY))
        source = fetch.fetch_model(self.model_dir, "m1", [PRIMARY, BACKUP], opener=opener)
        self.assertEqual(source, PRIMARY)
        self.assertEqual(self._listing(), ["m1"])  # no staging dir left behind
        for checkpoint in fetch.CHECKPOINTS:
            self.assertTrue(
                fetch.checkpoint_is_valid(os.path.join(self.model_dir, "m1", checkpoint))
            )
        with open(os.path.join(self.model_dir, "m1", "source.txt")) as fh:
            self.assertIn(f"url\t{PRIMARY}", fh.read())
        self.assertFalse(any(u.startswith(BACKUP) for u in opener.requests))

    def test_primary_failure_falls_back_to_backup(self):
        opener = _Opener(_serve(BACKUP))
        source = fetch.fetch_model(
            self.model_dir, "m1", [PRIMARY, BACKUP], opener=opener, attempts=1
        )
        self.assertEqual(source, BACKUP)
        self.assertEqual(self._listing(), ["m1"])

    def test_invalid_content_is_never_published(self):
        served = _serve(PRIMARY, full=b"<html>not a model</html>")
        opener = _Opener(served)
        with self.assertRaises(RuntimeError) as ctx:
            fetch.fetch_model(self.model_dir, "m1", [PRIMARY], opener=opener, attempts=1)
        self.assertIn("not a valid PyTorch checkpoint", str(ctx.exception))
        self.assertEqual(self._listing(), [])  # staging removed, nothing published

    def test_existing_model_replaced_only_after_validation(self):
        old = os.path.join(self.model_dir, "m1", "pileup.pt")
        os.makedirs(os.path.dirname(old))
        with open(old, "wb") as fh:
            fh.write(b"stale tensorflow leftover")
        # a failing fetch keeps the old directory
        with self.assertRaises(RuntimeError):
            fetch.fetch_model(self.model_dir, "m1", [PRIMARY], opener=_Opener({}), attempts=1)
        with open(old, "rb") as fh:
            self.assertEqual(fh.read(), b"stale tensorflow leftover")
        # a succeeding fetch replaces it wholesale
        fetch.fetch_model(self.model_dir, "m1", [PRIMARY], opener=_Opener(_serve(PRIMARY)))
        self.assertTrue(fetch.checkpoint_is_valid(old))
        self.assertEqual(self._listing(), ["m1"])

    def test_retries_then_reports_every_attempt(self):
        opener = _Opener({})
        with self.assertRaises(RuntimeError) as ctx:
            fetch.fetch_model(self.model_dir, "m1", [PRIMARY, BACKUP], opener=opener, attempts=2)
        message = str(ctx.exception)
        self.assertIn("attempt 2", message)
        self.assertIn(PRIMARY, message)
        self.assertIn(BACKUP, message)
        self.assertEqual(len(opener.requests), 4)  # 2 sources x 2 attempts, first file each

    def test_main_exit_codes(self):
        from unittest.mock import patch

        with patch.object(fetch.urllib.request, "urlopen", _Opener(_serve(PRIMARY))):
            rc = fetch.main(["--model-dir", self.model_dir, "--name", "m1", "--url", PRIMARY])
        self.assertEqual(rc, 0)
        with patch.object(fetch.urllib.request, "urlopen", _Opener({})):
            rc = fetch.main(
                ["--model-dir", self.model_dir, "--name", "m2", "--url", PRIMARY, "--attempts", "1"]
            )
        self.assertEqual(rc, 1)


if __name__ == "__main__":
    unittest.main()
