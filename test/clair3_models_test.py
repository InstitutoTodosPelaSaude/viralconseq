"""Tests for viralconseq.clair3_models: basecall-tag parsing, the artic-equivalent
model selection, per-sample resolution over real FASTQ files, and checkpoint
validation."""

import gzip
import os
import tempfile
import unittest

from viralconseq import clair3_models
from viralconseq.constants import Clair3Models
from viralconseq.exceptions import Clair3ModelMixedError, Clair3ModelUnresolvedError

# Worked examples, several checked against artic 1.11.2's choose_model.
EXPECTED = {
    "dna_r10.4.1_e8.2_400bps_hac@v5.0.0": "r1041_e82_400bps_hac_v500",
    "dna_r10.4.1_e8.2_400bps_hac@v4.2.0": "r1041_e82_400bps_hac_v420",
    "dna_r10.4.1_e8.2_400bps_hac@v5.2.0": "r1041_e82_400bps_hac_v520",
    "dna_r10.4.1_e8.2_400bps_hac@v6.0.0": "r1041_e82_400bps_hac_v600",
    "dna_r10.4.1_e8.2_400bps_sup@v4.3.0": "r1041_e82_400bps_sup_v430",
    "dna_r10.4.1_e8.2_400bps_sup@v5.0.0": "r1041_e82_400bps_sup_v500",
    "dna_r10.4.1_e8.2_260bps_hac@v4.1.0": "r1041_e82_260bps_hac_v410",
    "dna_r9.4.1_e8_hac@v3.3": "r941_prom_hac_g360+g422",
    "dna_r9.4.1_e8_sup@v3.3": "r941_prom_sup_g5014",
    "dna_r9.4.1_450bps_hac": "r941_prom_hac_g360+g422",  # 4 tokens, no '@'
}

FIXTURE_HEADER = (
    "a19879db-1909-4b73-8a74-adc0549bcbb4 runid=a90b74d226 read=126 ch=292 "
    "start_time=2022-07-26T18:25:23Z flow_cell_id=FAR87659 barcode=barcode05"
)


def _fastq(path, headers, gz=True):
    body = "".join(f"@{h}\nACGT\n+\nIIII\n" for h in headers)
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(body)
    else:
        with open(path, "w") as fh:
            fh.write(body)
    return path


class Test_BasecallModelId(unittest.TestCase):
    def test_key_value_header(self):
        h = "r1 runid=x basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.0.0 barcode=barcode01"
        self.assertEqual(clair3_models.basecall_model_id(h), "dna_r10.4.1_e8.2_400bps_hac@v5.0.0")

    def test_rg_tag_header(self):
        h = "r1 RG:Z:8f2a1c_dna_r10.4.1_e8.2_400bps_hac@v5.2.0_barcode01 qs:i:20"
        self.assertEqual(clair3_models.basecall_model_id(h), "dna_r10.4.1_e8.2_400bps_hac@v5.2.0")

    def test_untagged_header(self):
        self.assertIsNone(clair3_models.basecall_model_id(FIXTURE_HEADER))


class Test_ModelForBasecallId(unittest.TestCase):
    def test_worked_examples(self):
        for model_id, expected in EXPECTED.items():
            with self.subTest(model_id=model_id):
                self.assertEqual(clair3_models.model_for_basecall_id(model_id), expected)

    def test_every_result_is_a_manifest_name(self):
        for model_id in EXPECTED:
            self.assertIn(clair3_models.model_for_basecall_id(model_id), Clair3Models.MANIFEST)

    def test_guppy_guard(self):
        """sup@v5.1.0 has no versioned model; the only survivors are Guppy-era."""
        with self.assertRaises(Clair3ModelUnresolvedError) as ctx:
            clair3_models.model_for_basecall_id("dna_r10.4.1_e8.2_400bps_sup@v5.1.0")
        self.assertIn("Guppy-era", str(ctx.exception))

    def test_unknown_chemistry_and_shape(self):
        with self.assertRaises(Clair3ModelUnresolvedError):
            clair3_models.model_for_basecall_id("dna_r12.0.0_e9_400bps_hac@v1.0.0")
        with self.assertRaises(Clair3ModelUnresolvedError):
            clair3_models.model_for_basecall_id("weird@v1")

    def test_manifest_names_are_valid_and_have_urls(self):
        for name in Clair3Models.MANIFEST:
            self.assertRegex(name, Clair3Models.NAME_RE)
            primary, backup = Clair3Models.urls_for(name)
            self.assertTrue(primary.startswith(Clair3Models.PRIMARY_BASE + "/"))
            self.assertTrue(primary.endswith("/" + name))
            self.assertTrue(backup.endswith("/" + name))
        for name in Clair3Models.DEFAULT_SETUP_MODELS:
            self.assertIn(name, Clair3Models.MANIFEST)


class Test_ResolveFromFastq(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def _path(self, name="reads.fastq.gz"):
        return os.path.join(self.tmp, name)

    def test_uniform_tags_resolve(self):
        h = "r{i} basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.0.0"
        path = _fastq(self._path(), [h.format(i=i) for i in range(5)])
        self.assertEqual(clair3_models.resolve_model_from_fastq(path), "r1041_e82_400bps_hac_v500")

    def test_plain_fastq_and_rg_tags(self):
        h = "r{i} RG:Z:uuid_dna_r10.4.1_e8.2_400bps_sup@v5.0.0_barcode03"
        path = _fastq(self._path("reads.fastq"), [h.format(i=i) for i in range(3)], gz=False)
        self.assertEqual(clair3_models.resolve_model_from_fastq(path), "r1041_e82_400bps_sup_v500")

    def test_mixed_tags_are_an_error(self):
        headers = [
            "r1 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.0.0",
            "r2 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_sup@v5.0.0",
        ]
        with self.assertRaises(Clair3ModelMixedError):
            clair3_models.resolve_model_from_fastq(_fastq(self._path(), headers))

    def test_untagged_reads_name_the_flag(self):
        path = _fastq(self._path(), [FIXTURE_HEADER] * 3)
        with self.assertRaises(Clair3ModelUnresolvedError) as ctx:
            clair3_models.resolve_model_from_fastq(path)
        self.assertIn("--clair3-model", str(ctx.exception))
        self.assertIn("r941_prom_hac_g360+g422", str(ctx.exception))

    def test_empty_file_is_none(self):
        path = _fastq(self._path(), [])
        self.assertIsNone(clair3_models.resolve_model_from_fastq(path))

    def test_only_first_n_reads_are_read(self):
        h_ok = "r{i} basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.0.0"
        headers = [h_ok.format(i=i) for i in range(20)] + [FIXTURE_HEADER]
        path = _fastq(self._path(), headers)
        # the 21st, untagged read is never looked at
        self.assertEqual(
            clair3_models.resolve_model_from_fastq(path, n=20), "r1041_e82_400bps_hac_v500"
        )

    def test_models_by_sample(self):
        tagged = _fastq(
            self._path("a.fastq.gz"),
            ["r1 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.0.0"],
        )
        empty = _fastq(self._path("b.fastq.gz"), [])
        samples = {"a": [tagged], "b": [empty], "c": [tagged]}
        out = clair3_models.models_by_sample(samples, "auto")
        self.assertEqual(
            out, {"a": "r1041_e82_400bps_hac_v500", "b": None, "c": "r1041_e82_400bps_hac_v500"}
        )
        out = clair3_models.models_by_sample(samples, {"a": "x_model", "c": "auto"})
        self.assertEqual(out["a"], "x_model")
        self.assertIsNone(out["b"])  # unlisted -> auto -> empty file
        self.assertEqual(out["c"], "r1041_e82_400bps_hac_v500")
        out = clair3_models.models_by_sample(samples, "r941_prom_hac_g360+g422")
        self.assertEqual(set(out.values()), {"r941_prom_hac_g360+g422"})


class Test_Checkpoints(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def test_fake_checkpoint_is_valid_and_others_are_not(self):
        good = os.path.join(self.tmp, "m", "pileup.pt")
        clair3_models.write_fake_checkpoint(good)
        self.assertTrue(clair3_models.checkpoint_is_valid(good))
        empty = os.path.join(self.tmp, "m", "empty.pt")
        open(empty, "w").close()
        self.assertFalse(clair3_models.checkpoint_is_valid(empty))
        text = os.path.join(self.tmp, "m", "text.pt")
        with open(text, "w") as fh:
            fh.write("not a zip")
        self.assertFalse(clair3_models.checkpoint_is_valid(text))
        self.assertFalse(clair3_models.checkpoint_is_valid(os.path.join(self.tmp, "absent.pt")))

    def test_missing_checkpoints(self):
        self.assertEqual(clair3_models.missing_checkpoints(self.tmp, "nope"), ["<directory>"])
        os.makedirs(os.path.join(self.tmp, "m"))
        clair3_models.write_fake_checkpoint(os.path.join(self.tmp, "m", "pileup.pt"))
        self.assertEqual(clair3_models.missing_checkpoints(self.tmp, "m"), ["full_alignment.pt"])
        clair3_models.write_fake_checkpoint(os.path.join(self.tmp, "m", "full_alignment.pt"))
        self.assertEqual(clair3_models.missing_checkpoints(self.tmp, "m"), [])


class Test_SampleFastqHeaders(unittest.TestCase):
    def test_reads_headers_and_tolerates_blank_lines(self):
        from viralconseq.integrity import sample_fastq_headers

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "r.fastq")
            with open(path, "w") as fh:
                fh.write("\n@r1 a=b\nACGT\n+\nIIII\n\n@r2 c=d\nAC\n+\nII\n")
            self.assertEqual(sample_fastq_headers(path), ["r1 a=b", "r2 c=d"])
            self.assertEqual(sample_fastq_headers(path, n=1), ["r1 a=b"])


class Test_SingleCheckpointDefinition(unittest.TestCase):
    """``viralconseq setup`` and the nanopore pre-flight check must never
    disagree about what a complete Clair3 model is, so the predicate and the
    checkpoint names have exactly one definition (in the download script, which
    is stdlib-only and cannot import the package)."""

    def test_the_validator_and_the_names_are_shared(self):
        from viralconseq.constants import Clair3Models
        from viralconseq.scripts.python import fetch_clair3_model

        self.assertIs(clair3_models.checkpoint_is_valid, fetch_clair3_model.checkpoint_is_valid)
        self.assertIs(Clair3Models.CHECKPOINTS, fetch_clair3_model.CHECKPOINTS)


if __name__ == "__main__":
    unittest.main()
