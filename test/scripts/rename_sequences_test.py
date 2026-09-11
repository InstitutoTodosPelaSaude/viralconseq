import os
import tempfile
import unittest
from unittest.mock import mock_open, patch

from viralconseq.scripts.python.rename_sequences import rename_sequences


class TestRenameSequences(unittest.TestCase):
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data=">header\nATCG\nCGTA\nCCGC\nAACG\n",
    )
    def test_rename_sequences(self, mock_open):
        input = "assembly/consensus/final_consensus/1234.consensus.fasta"
        output = "assembly/consensus/final_consensus/1234.consensus.renamed.fasta"

        rename_sequences(input, output)

        mock_open.assert_any_call("assembly/consensus/final_consensus/1234.consensus.fasta")
        mock_open.assert_any_call(
            "assembly/consensus/final_consensus/1234.consensus.renamed.fasta", "w"
        )
        handle = mock_open()
        handle.write.assert_called_once_with(">1234\nATCG\nCGTA\nCCGC\nAACG\n")


class TestRenameSequencesMultiRecord(unittest.TestCase):
    """A multi-contig FASTA must keep one record per contig, not fuse them."""

    def _run(self, fasta_text):
        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "1234.consensus.fasta")
            out = os.path.join(tmp, "1234.consensus.renamed.fasta")
            with open(inp, "w") as fh:
                fh.write(fasta_text)
            rename_sequences(inp, out)
            with open(out) as fh:
                return fh.read()

    def test_multi_record_preserved(self):
        result = self._run(">segA some description\nACGT\n>segB\nTTTT\n")
        # Two records survive, each headed by sample + original contig id.
        self.assertEqual(result, ">1234_segA\nACGT\n>1234_segB\nTTTT\n")
        self.assertEqual(result.count(">"), 2)

    def test_single_record_uses_bare_sample_name(self):
        result = self._run(">whatever\nACGTACGT\n")
        self.assertEqual(result, ">1234\nACGTACGT\n")


if __name__ == "__main__":
    unittest.main()
