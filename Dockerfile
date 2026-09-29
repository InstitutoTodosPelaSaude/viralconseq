FROM continuumio/miniconda3:25.1.1-2
LABEL version="0.2.1" \
      description="Reference-guided consensus sequence inference for viral high-throughput sequencing data."

WORKDIR /app

COPY environment.yml /app/environment.yml
# Install the base environment (viralconseq + core runtime deps).
RUN conda env create --quiet -f environment.yml && conda clean -a -y

ENV PATH=/opt/conda/envs/viralconseq/bin:$PATH

COPY . /app/viralconseq

WORKDIR /app/viralconseq
RUN pip install . && rm -rf /root/.cache/pip
RUN viralconseq --version > /app/viralconseq-version.txt

# Run as an unprivileged user rather than root.
RUN useradd --create-home --uid 1000 viralconseq \
    && chown -R viralconseq:viralconseq /app
USER viralconseq

WORKDIR /tmp/
ENTRYPOINT ["viralconseq"]
CMD ["--help"]

# NOTE: per-rule conda envs are still created on first run and the viralQC
# databases (`viralconseq setup --viralqc-db`) are not baked in either, so a
# `docker run ... consensus ...` needs either a mounted database
# (`-v /host/viralqc-db:/db -e VIRALCONSEQ_VIRALQC_DB=/db`) or `--no-run-viralqc`. To
# pre-build the envs into the image, run `viralconseq setup --pipelines all
# --conda-prefix <fixed path> --skip-viralqc-db` here and make the pipeline reuse that same --conda-prefix at runtime.
# This needs a fixed, world-readable prefix (not $HOME-derived) so the non-root
# runtime user can reuse the root-built envs; it requires build-time network
# access, so it is intentionally left as a separate change.
