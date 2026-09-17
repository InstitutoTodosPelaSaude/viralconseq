"""Record the versions of the tools installed in the current conda environment.

Run by the ``versions_<env>`` rules in ``rules/provenance.smk`` inside each
per-rule conda environment (where the ``viralconseq`` package is NOT
installed), so this script is standard-library only and takes every fact as a
command-line argument::

    python tool_versions.py --tools minimap2 samtools --python-dists pandas \\
        --output other/versions/alignment.tsv

Each requested tool is probed with a known command and the first
``major.minor[...]`` token of its stdout+stderr is taken as the version. A tool
that is not on PATH is an environment error and exits 1; output that holds no
version token is recorded verbatim (``unparsed: <first line>``) so the file is
still written and the gap is visible.
"""

import argparse
import importlib.metadata
import re
import shutil
import subprocess
import sys
from typing import Dict, List, Optional, Sequence

# tool name -> command that prints its version. Tools that print usage (and exit
# non-zero) when called bare, like GSAlign, still yield a version string, so the
# exit status is deliberately ignored.
PROBES: Dict[str, List[str]] = {
    "minimap2": ["minimap2", "--version"],
    "samtools": ["samtools", "--version"],
    "bedtools": ["bedtools", "--version"],
    "GSAlign": ["GSAlign"],
    "gofasta": ["gofasta", "--version"],
    "fastp": ["fastp", "--version"],
    "multiqc": ["multiqc", "--version"],
    "lofreq": ["lofreq", "version"],
    "bcftools": ["bcftools", "--version"],
    "clair3": ["run_clair3.sh", "--version"],
    "python": ["python", "--version"],
    "nextclade": ["nextclade", "--version"],
    "blastn": ["blastn", "-version"],
}

# First "N.N..." token not glued to a word or dot before it. Accepts
# "2.30-r1287", "v2.31.1", "2.16.0+", "1.23.1".
VERSION_RE = re.compile(r"(?<![\w.])v?(\d+\.\d+[\w.+-]*)")


def parse_version(text: str) -> Optional[str]:
    match = VERSION_RE.search(text)
    return match.group(1) if match else None


def probe_tool(tool: str) -> str:
    """Return the version string of ``tool`` or ``unparsed: <first line>``.

    Raises:
        FileNotFoundError: If the probe command's executable is not on PATH.
    """
    try:
        command = PROBES[tool]
    except KeyError:
        raise SystemExit(f"tool_versions.py: no probe defined for {tool!r}")
    if shutil.which(command[0]) is None:
        raise FileNotFoundError(command[0])
    result = subprocess.run(command, capture_output=True, text=True, check=False, timeout=120)
    output = (result.stdout or "") + "\n" + (result.stderr or "")
    version = parse_version(output)
    if version is not None:
        return version
    first_line = next((line for line in output.splitlines() if line.strip()), "")
    return f"unparsed: {first_line.strip()}"


def probe_dist(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tools", nargs="*", default=[], help="tool names from PROBES")
    parser.add_argument(
        "--python-dists", nargs="*", default=[], help="Python distributions to look up"
    )
    parser.add_argument("--output", required=True, help="TSV to write (tool<TAB>version)")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    rows: List[List[str]] = [["tool", "version"]]
    missing: List[str] = []
    for tool in args.tools:
        try:
            rows.append([tool, probe_tool(tool)])
        except FileNotFoundError as exc:
            missing.append(f"{tool} ({exc})")
    for dist in args.python_dists:
        rows.append([dist, probe_dist(dist)])
    if missing:
        print(
            "tool_versions.py: not on PATH in this environment: " + ", ".join(missing),
            file=sys.stderr,
        )
        return 1
    with open(args.output, "w") as handle:
        for row in rows:
            handle.write("\t".join(row) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
