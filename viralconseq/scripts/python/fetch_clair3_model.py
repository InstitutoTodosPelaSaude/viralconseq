"""Download one Clair3 2.x model (``pileup.pt`` + ``full_alignment.pt``).

Used by ``viralconseq setup --clair3-models`` (imported, ``fetch_model``) and
runnable on its own::

    python fetch_clair3_model.py --model-dir ~/.cache/viralconseq/clair3-models \\
        --name r941_prom_hac_g360+g422 \\
        --url https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/r941_prom_hac_g360+g422 \\
        --backup-url https://artic-example-datasets.s3.climb.ac.uk/clair3-models/r941_prom_hac_g360+g422

Standard library only (no ``viralconseq`` import), so it also works from a
bare checkout. Each checkpoint is streamed into a staging directory next to
the destination, validated as a PyTorch zip archive (a truncated download or
an HTML error page is not a model), and the staged directory is renamed into
place only when both files pass, so a model directory is either complete or
absent. A ``source.txt`` records where and when the files came from.
"""

import argparse
import datetime
import os
import shutil
import sys
import tempfile
import urllib.error
import urllib.request
import zipfile
from typing import Callable, List, Optional, Sequence

CHECKPOINTS = ("pileup.pt", "full_alignment.pt")


def checkpoint_is_valid(path: str) -> bool:
    """A loadable Clair3 2.x checkpoint: a zip with ``*/data.pkl`` and a ``data/`` payload."""
    if not os.path.isfile(path):
        return False
    try:
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            return (
                any(n.endswith("/data.pkl") for n in names)
                and any("/data/" in n for n in names)
                and archive.testzip() is None
            )
    except (zipfile.BadZipFile, OSError):
        return False


def _download(url: str, destination: str, timeout: float, opener: Callable) -> None:
    with opener(url, timeout=timeout) as response, open(destination, "wb") as handle:
        shutil.copyfileobj(response, handle, 1 << 20)


def fetch_model(
    model_dir: str,
    name: str,
    urls: Sequence[str],
    *,
    timeout: float = 60.0,
    attempts: int = 2,
    opener: Optional[Callable] = None,
    log: Callable[[str], None] = lambda message: None,
) -> str:
    """Fetch ``name`` into ``<model_dir>/<name>/`` from the first ``urls`` base
    that serves both valid checkpoints.

    Returns:
        The base URL the files came from.

    Raises:
        RuntimeError: When every source failed; the message lists each attempt.
    """
    # Resolved at call time so a test can patch urllib.request.urlopen.
    opener = opener or urllib.request.urlopen
    os.makedirs(model_dir, exist_ok=True)
    final = os.path.join(model_dir, name)
    errors: List[str] = []
    for base in urls:
        for attempt in range(1, attempts + 1):
            staging = tempfile.mkdtemp(prefix=f".{name}.fetch-", dir=model_dir)
            staged = os.path.join(staging, name)
            os.makedirs(staged)
            try:
                for checkpoint in CHECKPOINTS:
                    url = f"{base.rstrip('/')}/{checkpoint}"
                    log(f"downloading {url}")
                    _download(url, os.path.join(staged, checkpoint), timeout, opener)
                    if not checkpoint_is_valid(os.path.join(staged, checkpoint)):
                        raise RuntimeError(f"{url} is not a valid PyTorch checkpoint")
                stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
                with open(os.path.join(staged, "source.txt"), "w") as handle:
                    handle.write(f"url\t{base}\nfetched_utc\t{stamp}\n")
                # Publish: move any previous directory aside into the staging
                # dir (removed with it) and rename the validated one into place.
                if os.path.isdir(final):
                    os.rename(final, os.path.join(staging, "previous"))
                os.rename(staged, final)
                shutil.rmtree(staging, ignore_errors=True)
                return base
            except (OSError, urllib.error.URLError, RuntimeError, ValueError) as exc:
                errors.append(f"{base} (attempt {attempt}): {exc}")
                log(f"failed: {base} (attempt {attempt}): {exc}")
                shutil.rmtree(staging, ignore_errors=True)
    raise RuntimeError(
        f"could not fetch Clair3 model {name} from any source:\n  " + "\n  ".join(errors)
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument(
        "--url", required=True, help="base URL holding pileup.pt and full_alignment.pt"
    )
    parser.add_argument("--backup-url", default=None, help="second base URL tried on failure")
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--attempts", type=int, default=2)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    urls = [args.url] + ([args.backup_url] if args.backup_url else [])
    try:
        source = fetch_model(
            args.model_dir,
            args.name,
            urls,
            timeout=args.timeout,
            attempts=args.attempts,
            log=lambda message: print(message, file=sys.stderr),
        )
    except RuntimeError as exc:
        print(f"fetch_clair3_model.py: {exc}", file=sys.stderr)
        return 1
    print(f"{args.name}: fetched from {source} into {os.path.join(args.model_dir, args.name)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
