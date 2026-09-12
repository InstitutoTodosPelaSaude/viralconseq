"""Central logging configuration for viralconseq.

Historically the pipeline created module loggers (``logging.getLogger(__name__)``)
and emitted ``logger.info``/``logger.error`` everywhere, but never configured a
handler. Python's last-resort handler only prints ``WARNING`` and above,
unformatted, to stderr, so a failed run could print nothing useful. This module
installs a single stderr handler, attaches a per-run correlation id to every
record, and supports human-readable or JSON output — the minimum an embedding
service needs to trace a job's logs.
"""

from __future__ import annotations

import json
import logging
import sys
import uuid
from typing import Optional

DEFAULT_LEVEL = "INFO"


def new_run_id() -> str:
    """Return a short, unique-enough correlation id for a single run."""
    return uuid.uuid4().hex[:12]


class _RunIdFilter(logging.Filter):
    """Stamp every record with the current run id so formatters can show it."""

    def __init__(self, run_id: str):
        super().__init__()
        self.run_id = run_id

    def filter(self, record: logging.LogRecord) -> bool:
        record.run_id = self.run_id
        return True


class _JsonFormatter(logging.Formatter):
    """Emit one JSON object per log record (for ingestion by a service)."""

    def format(self, record: logging.LogRecord) -> str:
        payload = {
            "ts": self.formatTime(record, "%Y-%m-%dT%H:%M:%S"),
            "level": record.levelname,
            "logger": record.name,
            "run_id": getattr(record, "run_id", None),
            "message": record.getMessage(),
        }
        if record.exc_info:
            payload["exception"] = self.formatException(record.exc_info)
        return json.dumps(payload)


def _coerce_level(level) -> int:
    if isinstance(level, int):
        return level
    return getattr(logging, str(level).upper(), logging.INFO)


def configure_logging(
    level=DEFAULT_LEVEL,
    json_logs: bool = False,
    run_id: Optional[str] = None,
) -> str:
    """Install a single stderr log handler on the root logger.

    Idempotent: existing handlers are removed first, so calling this more than
    once (e.g. per job in a long-lived service process) does not duplicate log
    lines.

    Args:
        level: Log level name or int (default ``"INFO"``).
        json_logs: If ``True``, emit structured JSON instead of text.
        run_id: Correlation id to stamp on every record. Generated if omitted.

    Returns:
        The run id in use (generated if none was passed).
    """
    run_id = run_id or new_run_id()

    root = logging.getLogger()
    for handler in list(root.handlers):
        root.removeHandler(handler)

    handler = logging.StreamHandler(stream=sys.stderr)
    handler.addFilter(_RunIdFilter(run_id))
    if json_logs:
        handler.setFormatter(_JsonFormatter())
    else:
        handler.setFormatter(
            logging.Formatter(
                "%(asctime)s %(levelname)s [%(name)s] [run:%(run_id)s] %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )
        )

    root.addHandler(handler)
    root.setLevel(_coerce_level(level))
    return run_id
