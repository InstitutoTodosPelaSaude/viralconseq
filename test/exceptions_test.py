"""Tests for the structured error codes on viralconseq exceptions."""

import unittest

from viralconseq.exceptions import (
    SampleSheetError,
    ValidationError,
    ViralConseqError,
    ViralConseqFileNotFoundError,
    ViralQCDatabaseNotFoundError,
)


class Test_StructuredErrors(unittest.TestCase):
    def test_base_has_code(self):
        self.assertEqual(ViralConseqError("boom").code, "viralconseq_error")

    def test_subclasses_have_distinct_codes(self):
        self.assertEqual(ValidationError("x").code, "validation_error")
        self.assertEqual(ViralConseqFileNotFoundError("x").code, "file_not_found")
        self.assertEqual(SampleSheetError("x").code, "sample_sheet_error")
        self.assertEqual(ViralQCDatabaseNotFoundError("x").code, "viralqc_database_not_found")
        self.assertIsInstance(ViralQCDatabaseNotFoundError("x"), ViralConseqError)

    def test_message_preserved(self):
        err = ValidationError("bad thing")
        self.assertEqual(str(err), "bad thing")
        self.assertEqual(err.message, "bad thing")

    def test_to_dict_is_structured(self):
        err = ViralConseqFileNotFoundError("missing ref")
        self.assertEqual(
            err.to_dict(),
            {
                "error": "ViralConseqFileNotFoundError",
                "code": "file_not_found",
                "message": "missing ref",
            },
        )

    def test_per_instance_code_override(self):
        err = ValidationError("x", code="custom_code")
        self.assertEqual(err.code, "custom_code")


if __name__ == "__main__":
    unittest.main()
