"""Choose and check Clair3 models for nanopore runs.

Two concerns, both pure functions over :class:`viralconseq.constants.Clair3Models`:

* **Selection** — read the basecall model tag Dorado/MinKNOW write into every
  FASTQ header (``basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.0.0``
  or the newer ``RG:Z:<uuid>_dna_r10.4.1_e8.2_400bps_hac@v5.0.0_barcode01``)
  and map it to a model name. The mapping is a stdlib port of
  ``artic.utils.choose_model`` (artic 1.11.2) so viralconseq picks the same
  model artic would, without importing artic or Biopython. Several reads are
  inspected, not one, so a barcode merged from two basecalling runs is caught.

* **Completeness** — a model directory is usable only if both PyTorch
  checkpoints are valid zip archives holding a ``data.pkl`` member. Clair3
  itself only checks that the directory exists and then dies inside
  ``torch.load`` with a bare exit 2, so this is checked before any work starts.
"""

import os
import re
import zipfile
from typing import Dict, List, Optional

from viralconseq import integrity
from viralconseq.constants import Clair3Models
from viralconseq.exceptions import Clair3ModelMixedError, Clair3ModelUnresolvedError

_RG_MODEL_RE = re.compile(r"(?:dna|rna)\w*_[^@\s]+@v[\d.]+")

HINT = "pass --clair3-model NAME (see 'viralconseq setup --help' for the model list)"


def basecall_model_id(header: str) -> Optional[str]:
    """The basecall model id in a FASTQ header line, or ``None``.

    ``header`` is the header without its leading ``@`` (read id first, then
    the description tokens).
    """
    tokens = header.split()[1:]
    kv = dict(t.split("=", 1) for t in tokens if "=" in t)
    if "basecall_model_version_id" in kv:
        return kv["basecall_model_version_id"]
    for token in tokens:
        if token.startswith("RG:Z:"):
            match = _RG_MODEL_RE.search(token[5:])
            if match:
                return match.group(0)
    return None


def model_for_basecall_id(model_id: str) -> str:
    """Map a Dorado basecall model id to a Clair3 model name.

    Port of artic 1.11.2 ``choose_model``: 4-token ids
    (``dna_r9.4.1_e8_hac@v3.3``) filter the manifest by pore type, 5-token ids
    (``dna_r10.4.1_e8.2_400bps_hac@v5.0.0``) by ``<pore>_<kit>_<speed>``; then
    by preset (``hac``/``sup``/``fast``), then the first manifest entry whose
    name holds the basecaller version (``v5.0.0`` -> ``v500``). Guppy-era
    candidates (``_g###``) are never a silent fallback for a Dorado tag.

    Raises:
        Clair3ModelUnresolvedError: No, or no unambiguous, model.
    """
    names = list(Clair3Models.MANIFEST)
    tokens = model_id.split("_")
    if len(tokens) == 4:
        _, pore_type, _kit, model = tokens
        search = pore_type.replace(".", "")
    elif len(tokens) == 5:
        _, pore_type, kit_id, speed, model = tokens
        search = f"{pore_type}_{kit_id}_{speed}".replace(".", "")
    else:
        raise Clair3ModelUnresolvedError(
            f"basecall model id {model_id!r} has an unexpected shape; {HINT}."
        )
    if "@" in model:
        preset, version = model.split("@", 1)
        version = version.replace(".", "")
    else:
        preset, version = model, ""

    def none_found() -> Clair3ModelUnresolvedError:
        return Clair3ModelUnresolvedError(
            f"no Clair3 model in the manifest ({Clair3Models.MANIFEST_SOURCE}) matches "
            f"basecall model id {model_id!r}; {HINT}."
        )

    candidates = [name for name in names if search in name]
    if not candidates:
        raise none_found()
    if len(candidates) == 1:
        return candidates[0]
    candidates = [name for name in candidates if preset in name]
    if not candidates:
        raise none_found()
    if len(candidates) == 1:
        return candidates[0]
    if version:
        for name in candidates:
            if version in name:
                return name
    guppy = [name for name in candidates if name.split("_")[-1].startswith("g")]
    if guppy:
        raise Clair3ModelUnresolvedError(
            f"no versioned Clair3 model matches Dorado basecall model {model_id!r}; the only "
            f"candidates ({', '.join(guppy)}) are Guppy-era models and are incompatible with "
            f"Dorado data. Re-basecall with hac or sup, or {HINT}."
        )
    raise Clair3ModelUnresolvedError(
        f"basecall model id {model_id!r} matches several Clair3 models "
        f"({', '.join(candidates)}); {HINT}."
    )


def resolve_model_from_fastq(path: str, n: int = 20) -> Optional[str]:
    """The Clair3 model for the reads in ``path``, from their basecall tags.

    Inspects the first ``n`` reads. Returns ``None`` for a file with no reads
    (the caller decides what an empty sample gets).

    Raises:
        Clair3ModelUnresolvedError: No read carries a tag, or the tag maps to
            no model.
        Clair3ModelMixedError: The reads resolve to different models.
    """
    headers = integrity.sample_fastq_headers(path, n)
    if not headers:
        return None
    ids = {basecall_model_id(h) for h in headers} - {None}
    if not ids:
        raise Clair3ModelUnresolvedError(
            f"{path}: none of the first {len(headers)} reads carries a basecall_model_version_id "
            f"or RG:Z tag, so the Clair3 model cannot be chosen from the data (Guppy-era or "
            f"re-headered reads); {HINT}, e.g. r941_prom_hac_g360+g422 for R9.4.1 Guppy hac reads."
        )
    resolved = {model_for_basecall_id(i) for i in ids if i}
    if len(resolved) > 1:
        raise Clair3ModelMixedError(
            f"{path}: reads resolve to different Clair3 models {sorted(resolved)} (basecall ids "
            f"{sorted(i for i in ids if i)}); split the sample, or {HINT} to force one."
        )
    return resolved.pop()


def checkpoint_is_valid(path: str) -> bool:
    """True if ``path`` is a PyTorch checkpoint Clair3 2.x can load: a zip
    archive with a ``*/data.pkl`` member, a ``*/data/`` payload and no CRC
    errors. Detects truncated downloads and TensorFlow-era leftovers without
    unpickling anything."""
    if not os.path.isfile(path):
        return False
    try:
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            if not any(name.endswith("/data.pkl") for name in names):
                return False
            if not any("/data/" in name for name in names):
                return False
            return archive.testzip() is None
    except (zipfile.BadZipFile, OSError):
        return False


def missing_checkpoints(model_dir: str, name: str) -> List[str]:
    """Checkpoints of ``<model_dir>/<name>`` that are absent or invalid
    (``["<directory>"]`` when the model directory itself is missing)."""
    directory = os.path.join(model_dir, name)
    if not os.path.isdir(directory):
        return ["<directory>"]
    return [
        checkpoint
        for checkpoint in Clair3Models.CHECKPOINTS
        if not checkpoint_is_valid(os.path.join(directory, checkpoint))
    ]


def write_fake_checkpoint(path: str) -> None:
    """Write a minimal archive that passes :func:`checkpoint_is_valid`.

    For tests and for the placeholder model directory ``viralconseq setup``
    builds so Snakemake can walk the DAG before any real model exists.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("model/data.pkl", b"\x80\x02N.")
        archive.writestr("model/data/0", b"0")


def models_by_sample(
    samples: Dict[str, List[str]], requested: object, n: int = 20
) -> Dict[str, Optional[str]]:
    """Resolve the model for every sample.

    ``requested`` is the ``--clair3-model`` value: ``"auto"``, a model name, or
    a mapping ``{sample id: name or "auto"}``. Detection reads the sample's
    first FASTQ. A sample whose FASTQ holds no reads maps to ``None``.
    """
    if isinstance(requested, dict):
        per_sample = {s: str(requested.get(s, Clair3Models.AUTO)) for s in samples}
    else:
        per_sample = {s: str(requested) for s in samples}
    out: Dict[str, Optional[str]] = {}
    for sample, name in per_sample.items():
        if name == Clair3Models.AUTO:
            out[sample] = resolve_model_from_fastq(samples[sample][0], n)
        else:
            out[sample] = name
    return out
