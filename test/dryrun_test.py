import os
import re
import subprocess

import pytest

# Resolve paths relative to this file, not the process cwd, so the suite runs
# from anywhere (not only the repository root).
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(REPO_ROOT, "test", "dryrun_configs")
SCRIPTS_DIR = os.path.join(REPO_ROOT, "viralconseq", "scripts")
PLACEHOLDER_SCRIPT = os.path.join(CONFIG_DIR, "create_dryrun_placeholders.sh")


@pytest.fixture(scope="session", autouse=True)
def setup_dryrun_placeholders():
    """Runs the placeholder script once per test session."""
    if os.path.exists(PLACEHOLDER_SCRIPT):
        # Run from the repo root because the script uses relative paths like 'data/reads'.
        result = subprocess.run(
            ["bash", PLACEHOLDER_SCRIPT],
            check=True,
            capture_output=True,
            text=True,
            cwd=REPO_ROOT,
        )
        print(result.stdout)
    else:
        pytest.fail(f"Placeholder script not found at {PLACEHOLDER_SCRIPT}")


def get_dryrun_configs():
    """Discovers all .yaml files in the dryrun_configs directory."""
    if not os.path.exists(CONFIG_DIR):
        return []
    return [f for f in os.listdir(CONFIG_DIR) if f.endswith(".yaml")]


def get_workflow_file(config_name):
    """Maps a configuration file name to its corresponding Snakemake workflow file.

    A ``__variant`` suffix lets several configs target the same workflow, e.g.
    both ``consensus_illumina.yaml`` and ``consensus_illumina__isnv_primers.yaml``
    resolve to ``consensus_illumina.smk``.
    """
    # Example: consensus_illumina.yaml -> viralconseq/scripts/consensus_illumina.smk
    base_name = os.path.splitext(config_name)[0].split("__")[0]
    workflow_file = os.path.join(SCRIPTS_DIR, f"{base_name}.smk")
    return workflow_file


@pytest.mark.parametrize("config_filename", get_dryrun_configs())
def test_snakemake_dryrun(config_filename):
    """Executes a Snakemake dry-run for a given configuration file."""
    config_path = os.path.join(CONFIG_DIR, config_filename)
    workflow_path = get_workflow_file(config_filename)

    # Check if the workflow file exists
    if not os.path.exists(workflow_path):
        pytest.fail(f"Workflow file not found: {workflow_path} for config {config_filename}")

    # Run Snakemake dry-run (-n) printing commands (-p). No explicit target:
    # every entry workflow marks ``rule all`` with ``default_target: True``, so
    # a bare ``snakemake -s <workflow>`` must plan the whole DAG even where
    # another rule (nanopore ``sanitize_reference``) is defined before it.
    cmd = [
        "snakemake",
        "-s",
        workflow_path,
        "--configfile",
        config_path,
        "-n",
        "-p",
        "--cores",
        "1",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO_ROOT)

    # Assert success (return code 0)
    assert (
        result.returncode == 0
    ), f"Snakemake dry-run failed for {config_filename}\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}"

    # The report is planned exactly when run_report is on.
    planned = re.search(r"^rule report:", result.stdout, re.M) is not None
    assert planned == (
        "__no_report" not in config_filename
    ), f"rule report {'planned' if planned else 'absent'} for {config_filename}\n{result.stdout}"


def test_missing_required_key_is_reported_at_parse_time(tmp_path):
    """A hand-edited YAML lacking an analysis key fails with the key named,
    before any rule is planned (rules/common.smk guard)."""
    import yaml

    with open(os.path.join(CONFIG_DIR, "consensus_illumina.yaml")) as fh:
        config = yaml.safe_load(fh)
    del config["minimum_depth"]
    broken = tmp_path / "broken.yaml"
    broken.write_text(yaml.safe_dump(config))

    cmd = [
        "snakemake",
        "-s",
        get_workflow_file("consensus_illumina.yaml"),
        "--configfile",
        str(broken),
        "-n",
        "--cores",
        "1",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO_ROOT)
    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert "missing required key(s): minimum_depth" in combined, combined


def test_unsafe_flag_string_is_rejected_at_parse_time(tmp_path):
    """A flag string with a shell metacharacter never reaches a shell."""
    import yaml

    with open(os.path.join(CONFIG_DIR, "consensus_illumina.yaml")) as fh:
        config = yaml.safe_load(fh)
    config["viralqc_extra_flags"] = "--x; rm -rf /"
    broken = tmp_path / "broken.yaml"
    broken.write_text(yaml.safe_dump(config))

    cmd = [
        "snakemake",
        "-s",
        get_workflow_file("consensus_illumina.yaml"),
        "--configfile",
        str(broken),
        "-n",
        "--cores",
        "1",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO_ROOT)
    assert result.returncode != 0
    assert "viralqc_extra_flags contains a shell metacharacter" in result.stdout + result.stderr
