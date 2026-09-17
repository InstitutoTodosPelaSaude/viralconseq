#!/usr/bin/env bash
# Set up viralconseq: create its conda environment, install the package into it, and
# run `viralconseq setup` (per-rule conda envs, viralQC databases, Clair3 models).
#
# Usage:
#   bash setup.sh                                  # env + install + viralconseq setup
#   bash setup.sh --no-setup                       # env + install only
#   bash setup.sh --skip-clair3-models             # no Clair3 models (Illumina-only sites)
#   bash setup.sh --skip-viralqc-db                # no viralQC databases (about 1 GB)
#   bash setup.sh --conda-prefix ~/envs/viralconseq --viralqc-db /data/viralqc-db
#   bash setup.sh --clair3-models r1041_e82_400bps_hac_v500 r1041_e82_400bps_sup_v500
#
# Options:
#   --env-name NAME          conda env to create/update  [$VIRALCONSEQ_ENV_NAME, viralconseq]
#   --cores N                threads for `viralconseq setup`  [$VIRALCONSEQ_CORES, 4]
#   --no-setup               stop after `pip install`; do not run `viralconseq setup`
#   --conda-prefix DIR       cache for the per-rule envs (viralconseq setup --conda-prefix)
#   --viralqc-db DIR         viralQC database directory (viralconseq setup --viralqc-db)
#   --skip-viralqc-db        do not download the viralQC databases
#   --clair3-model-dir DIR   Clair3 model directory (viralconseq setup --clair3-model-dir)
#   --clair3-models NAME...  Clair3 models to download; repeatable, space- or comma-
#                            separated; 'all' fetches the whole manifest
#   --skip-clair3-models     do not download Clair3 models
#   -h, --help               print this text
#
# conda, mamba or micromamba is found on PATH, or in the usual install roots. Set
# VIRALCONSEQ_CONDA_BIN to point at a particular one.
#
# Everything is idempotent: re-running updates the environment and the install, and
# `viralconseq setup` skips envs, databases and models that are already present.
#
# On WSL, keep the per-rule env cache (--conda-prefix) on the Linux filesystem, not
# under /mnt/c: conda environments and Snakemake metadata are several times slower
# there. The default, ~/.cache/viralconseq/conda-envs, already is.

set -euo pipefail

ENV_NAME="${VIRALCONSEQ_ENV_NAME:-viralconseq}"
CORES="${VIRALCONSEQ_CORES:-4}"
RUN_SETUP=1
CONDA_PREFIX_ARG=""
VIRALQC_DB=""
CLAIR3_MODEL_DIR=""
CLAIR3_MODELS=()
SKIP_VIRALQC_DB=0
SKIP_CLAIR3_MODELS=0
ENV_PREFIX=""

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

need_value() {
    if [[ $# -lt 2 || "$2" == --* ]]; then
        echo "ERROR: $1 requires a value" >&2
        exit 2
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-setup) RUN_SETUP=0; shift ;;
        --env-name) need_value "$@"; ENV_NAME="$2"; shift 2 ;;
        --cores) need_value "$@"; CORES="$2"; shift 2 ;;
        --conda-prefix) need_value "$@"; CONDA_PREFIX_ARG="$2"; shift 2 ;;
        --viralqc-db) need_value "$@"; VIRALQC_DB="$2"; shift 2 ;;
        --skip-viralqc-db) SKIP_VIRALQC_DB=1; shift ;;
        --clair3-model-dir) need_value "$@"; CLAIR3_MODEL_DIR="$2"; shift 2 ;;
        --clair3-models)
            need_value "$@"; shift
            # Swallow every following word that is not an option: both
            # `--clair3-models a b` and `--clair3-models a,b` (and repeats) work.
            while [[ $# -gt 0 && "$1" != --* && "$1" != -h ]]; do
                CLAIR3_MODELS+=("$1"); shift
            done ;;
        --skip-clair3-models) SKIP_CLAIR3_MODELS=1; shift ;;
        -h|--help) awk 'NR>1 && /^#/ {sub(/^# ?/, ""); print; next} NR>1 {exit}' \
                       "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "unknown option: $1 (see --help)" >&2; exit 2 ;;
    esac
done

if ! [[ "${CORES}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --cores must be a positive integer, got '${CORES}'" >&2
    exit 2
fi

# --- the conda front end --------------------------------------------------------------
#
# Anything conda-flavoured will do. mamba and micromamba share one command line; conda
# spells two things differently, so the flavour is recorded and branched on below.
#
# A non-interactive shell often has none of them on PATH even when the interactive one
# does -- `conda init` only edits the interactive part of the shell profile -- so fall
# back to the pointers a conda install exports, and then to the usual install roots.

# VIRALCONSEQ_CONDA_BIN, when set, is taken as given; otherwise search.
CONDA_BIN="${VIRALCONSEQ_CONDA_BIN:-}"
if [[ -n "${CONDA_BIN}" && ! -x "${CONDA_BIN}" ]]; then
    echo "ERROR: VIRALCONSEQ_CONDA_BIN=${CONDA_BIN} is not an executable" >&2
    exit 1
fi

if [[ -z "${CONDA_BIN}" ]]; then
    for candidate in \
            "$(command -v mamba 2>/dev/null || true)" \
            "$(command -v conda 2>/dev/null || true)" \
            "$(command -v micromamba 2>/dev/null || true)" \
            "${MAMBA_EXE:-}" \
            "${CONDA_EXE:-}" \
            "${HOME}"/{miniforge3,mambaforge,miniconda3,anaconda3}/bin/{mamba,conda} \
            /opt/{miniforge3,mambaforge,conda,miniconda3}/bin/{mamba,conda}; do
        if [[ -n "${candidate}" && -x "${candidate}" ]]; then
            CONDA_BIN="${candidate}"
            break
        fi
    done
fi

if [[ -z "${CONDA_BIN}" ]]; then
    echo "ERROR: found no conda, mamba or micromamba. Install miniforge first:" >&2
    echo "  https://github.com/conda-forge/miniforge" >&2
    echo "If one is installed, either put it on PATH or point VIRALCONSEQ_CONDA_BIN at it:" >&2
    echo "  VIRALCONSEQ_CONDA_BIN=~/miniforge3/bin/mamba bash setup.sh" >&2
    exit 1
fi

case "$(basename "${CONDA_BIN}")" in
    mamba|micromamba) FLAVOUR="mamba" ;;
    *)                FLAVOUR="conda" ;;
esac
echo "==> using ${CONDA_BIN}"

# Where the environment lives, empty if it does not exist yet. Environments are
# addressed by path rather than by name from here on: a machine with two installs (a
# miniforge and a micromamba root, say) lists environments belonging to both, and a name
# resolved against the wrong root fails with EnvironmentLocationNotFound.
env_prefix() {
    "${CONDA_BIN}" env list | awk -v n="${ENV_NAME}" '$1 == n {print $NF; exit}'
}

# Run a command inside the environment. This is `conda run`, not `conda activate`:
# activation needs a shell hook sourced into this process, which is exactly the thing
# that is missing whenever conda is not on PATH.
in_env() {
    local target=(--name "${ENV_NAME}")
    [[ -n "${ENV_PREFIX}" ]] && target=(--prefix "${ENV_PREFIX}")
    if [[ "${FLAVOUR}" == "conda" ]]; then
        "${CONDA_BIN}" run --no-capture-output "${target[@]}" "$@"
    else
        "${CONDA_BIN}" run "${target[@]}" "$@"
    fi
}

if [[ "${FLAVOUR}" == "conda" ]]; then
    priority="$("${CONDA_BIN}" config --show channel_priority 2>/dev/null | awk '{print $2}')"
else
    priority="$("${CONDA_BIN}" config list channel_priority 2>/dev/null | awk '{print $2}')"
fi
if [[ "${priority}" != "strict" ]]; then
    echo "NOTE: channel_priority is not 'strict'. With conda-forge and bioconda both"
    echo "      enabled that can produce inconsistent solves. Consider:"
    echo "        $(basename "${CONDA_BIN}") config --set channel_priority strict"
fi

# --- environment --------------------------------------------------------------------

yes_flag=()
[[ "${FLAVOUR}" == "mamba" ]] && yes_flag=(--yes)

ENV_PREFIX="$(env_prefix)"

if [[ -n "${ENV_PREFIX}" && -d "${ENV_PREFIX}" ]]; then
    echo "==> updating environment '${ENV_NAME}' (${ENV_PREFIX})"
    "${CONDA_BIN}" env update --prefix "${ENV_PREFIX}" --file "${here}/environment.yml" \
        --prune "${yes_flag[@]}"
else
    echo "==> creating environment '${ENV_NAME}'"
    "${CONDA_BIN}" env create --name "${ENV_NAME}" --file "${here}/environment.yml" \
        "${yes_flag[@]}"
    ENV_PREFIX="$(env_prefix)"
fi

echo "==> installing viralconseq (editable, with the dev extras) from ${here}"
in_env python -m pip install -e "${here}[dev]"

echo "==> $(in_env viralconseq --version)"

# --- viralconseq setup ------------------------------------------------------------------

setup_args=(--pipelines all --threads "${CORES}")
[[ -n "${CONDA_PREFIX_ARG}" ]] && setup_args+=(--conda-prefix "${CONDA_PREFIX_ARG}")
[[ -n "${VIRALQC_DB}" ]] && setup_args+=(--viralqc-db "${VIRALQC_DB}")
[[ "${SKIP_VIRALQC_DB}" -eq 1 ]] && setup_args+=(--skip-viralqc-db)
[[ -n "${CLAIR3_MODEL_DIR}" ]] && setup_args+=(--clair3-model-dir "${CLAIR3_MODEL_DIR}")
for model in "${CLAIR3_MODELS[@]+"${CLAIR3_MODELS[@]}"}"; do
    setup_args+=(--clair3-models "${model}")
done
[[ "${SKIP_CLAIR3_MODELS}" -eq 1 ]] && setup_args+=(--skip-clair3-models)

if [[ "${RUN_SETUP}" -eq 0 ]]; then
    echo
    echo "Skipped 'viralconseq setup'. Build the per-rule envs and fetch the databases later with:"
    echo "  viralconseq setup ${setup_args[*]}"
else
    echo "==> viralconseq setup ${setup_args[*]}"
    echo "    (per-rule conda envs; viralQC databases are about 1 GB and take 15-60 min)"
    in_env viralconseq setup "${setup_args[@]}"
fi

front_end="$(basename "${CONDA_BIN}")"
activate="${front_end} activate ${ENV_NAME}"
init_hint="${front_end} init bash"
[[ "${front_end}" == "micromamba" ]] && init_hint="micromamba shell init -s bash"

cat <<MSG

Done.

  ${activate}
  viralconseq create-samplesheet --input <run_dir> --output samples.csv
  viralconseq consensus illumina --sample-sheet samples.csv --reference ref.fasta \\
      --primer-scheme primers.bed --run-name run1 --config-file run1.yml --output results/

If activation reports that your shell is not set up for it, initialise it once with
\`${init_hint}\` and open a new terminal.

On WSL, keep the working directory and the per-rule env cache off /mnt/c if you can;
conda environments and Snakemake metadata are much slower there. The default
--conda-prefix (~/.cache/viralconseq/conda-envs) already lives under \$HOME for that
reason; pass a Linux-filesystem path if you override it.
MSG
