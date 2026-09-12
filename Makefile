.PHONY: test test-dryrun test-empirical lint format typecheck run-consensus install install-dev build build-docker run-docker

# Conda env used for the end-to-end suite. viralconseq pins python <3.12, so the
# empirical target runs inside this env rather than whatever interpreter invoked
# make (a python 3.12 base env would fail `pip install -e`). Override as needed:
#   make test-empirical VIRALCONSEQ_ENV=my-env
VIRALCONSEQ_ENV ?= viralconseq

test: install-dev
	python3 -m unittest discover ./test -p *test.py

test-dryrun: install-dev
	pytest test/dryrun_test.py -v

# Opt-in end-to-end suite: downloads real data, runs full pipelines. Runs inside
# the $(VIRALCONSEQ_ENV) conda env (see the variable above) so it works regardless
# of the interpreter that invoked make. Prerequisite: `viralconseq setup
# --pipelines all` once to build the per-rule envs. Override the data cache with
# VIRALCONSEQ_TEST_CACHE=/some/path.
test-empirical:
	@command -v conda >/dev/null 2>&1 || { \
		echo "ERROR: conda not found on PATH. Install conda/miniforge (see README)."; exit 1; }
	@conda run -n $(VIRALCONSEQ_ENV) python --version >/dev/null 2>&1 || { \
		echo "ERROR: conda env '$(VIRALCONSEQ_ENV)' not found. Create it with:"; \
		echo "    conda env create -n $(VIRALCONSEQ_ENV) -f environment.yml"; \
		echo "or override the name: make test-empirical VIRALCONSEQ_ENV=<name>"; exit 1; }
	conda run --no-capture-output -n $(VIRALCONSEQ_ENV) python -m pip install -e ".[dev]"
	conda run --no-capture-output -n $(VIRALCONSEQ_ENV) pytest test/empirical_test.py -v -m empirical

lint: install-dev
	black --check viralconseq/ test/
	ruff check viralconseq/ test/

# Type check (gating in CI).
typecheck: install-dev
	mypy viralconseq/

format: install-dev
	black viralconseq/ test/
	ruff check --fix viralconseq/ test/

install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"

# --no-run-viralqc keeps this developer target runnable without the viralQC
# databases; drop the flag once `viralconseq setup` has populated them.
run-consensus:
	viralconseq consensus illumina \
		--no-run-viralqc \
		--sample-sheet input/samplesheet.csv \
		--config-file output/config_consensus.yml \
		--run-name test-consensus \
		--reference input/references/SARS-CoV-2_RefSeq.fasta \
		--primer-scheme input/references/scheme.bed \
		--threads 1 \
		--threads-total 1 \
		--output output/test-consensus-example

# Build the sdist + wheel for PyPI (see RELEASING.md). Requires `build` + `twine`.
build:
	python -m build
	twine check dist/*

build-docker:
	docker build -t filiperomero2/viralconseq:latest .

run-docker:
	docker run --rm -i -t filiperomero2/viralconseq:latest
