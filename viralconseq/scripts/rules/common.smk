# Shared preamble for the four consensus workflows. Included FIRST by every
# entry-point .smk, before ``rule all`` and before any other include, because
# global wildcard constraints only apply to rules parsed after them.
#
# Provides:
# - exact-sample wildcard constraints (a sample wildcard matches only the ids
#   in config["samples"], so ids that are prefixes of one another can never
#   produce an ambiguous match)
# - the run log ``<output>/logs/run.log`` written by onstart/onsuccess/onerror
#
# Expects ``config["samples"]`` and ``config["output"]`` (with its trailing
# slash, as ConfigGenerator writes it).

import datetime
import os
import re

wildcard_constraints:
    sample="|".join(re.escape(s) for s in config["samples"]),
    segment=r"[^/]+",
    ref_key=r"[^/]+",


RUN_LOG = config["output"] + "logs/run.log"


def append_run_log(message):
    """Append one timestamped line to the run log. Never raises: this is
    provenance, and a full disk or read-only directory must not turn a
    finished analysis into a failure."""
    try:
        os.makedirs(os.path.dirname(RUN_LOG), exist_ok=True)
        stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        with open(RUN_LOG, "a") as handle:
            handle.write(f"{stamp}\t{message}\n")
    except OSError:
        pass


onstart:
    append_run_log(
        "start\tviralconseq "
        + str(config.get("viralconseq_version", "unknown"))
        + "\tdata=" + str(config.get("data", "unknown"))
    )
    append_run_log(f"samples\t{len(config['samples'])}")
    append_run_log(f"output\t{config['output']}")
    append_run_log(f"snakemake_log\t{log}")


onsuccess:
    append_run_log("end\tsuccess")
    print("viralconseq: workflow completed successfully.")


onerror:
    append_run_log("end\tfailed")
    print(f"viralconseq: workflow FAILED - see the Snakemake log: {log} and {RUN_LOG}")
