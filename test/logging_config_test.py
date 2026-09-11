"""Tests for viralconseq.logging_config."""

import json
import logging
import unittest

from viralconseq.logging_config import configure_logging, new_run_id


class Test_ConfigureLogging(unittest.TestCase):
    def tearDown(self):
        # Leave the root logger clean for other tests.
        root = logging.getLogger()
        for handler in list(root.handlers):
            root.removeHandler(handler)

    def test_installs_single_handler_and_is_idempotent(self):
        configure_logging(level="INFO")
        configure_logging(level="INFO")
        root = logging.getLogger()
        self.assertEqual(len(root.handlers), 1)

    def test_sets_level(self):
        configure_logging(level="DEBUG")
        self.assertEqual(logging.getLogger().level, logging.DEBUG)

    def test_returns_and_stamps_run_id(self):
        run_id = configure_logging(level="INFO", run_id="abc123")
        self.assertEqual(run_id, "abc123")
        handler = logging.getLogger().handlers[0]
        record = logging.LogRecord("x", logging.INFO, __file__, 1, "hi", None, None)
        for f in handler.filters:
            f.filter(record)
        self.assertEqual(record.run_id, "abc123")

    def test_json_logs_are_valid_json(self):
        configure_logging(level="INFO", json_logs=True, run_id="rid42")
        handler = logging.getLogger().handlers[0]
        record = logging.LogRecord("mylogger", logging.INFO, __file__, 1, "hello", None, None)
        for f in handler.filters:
            f.filter(record)
        payload = json.loads(handler.formatter.format(record))
        self.assertEqual(payload["message"], "hello")
        self.assertEqual(payload["run_id"], "rid42")
        self.assertEqual(payload["level"], "INFO")

    def test_new_run_id_is_unique_ish(self):
        self.assertNotEqual(new_run_id(), new_run_id())


if __name__ == "__main__":
    unittest.main()
