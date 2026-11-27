"""Translate released CPST names to and from their pinned historical paths."""

from functools import lru_cache
import json
from pathlib import Path
import re

_PATH_COMPONENTS = {
    "stage2_CarPS-VAE": "cpst",
    "07_cpsv_alpha": "07_cpst_alpha",
    "08_cpsv_chem-emb": "08_cpst_chem-emb",
    "09_cpsv_axis": "09_cpst_axis",
    "10_cpsv_dose": "10_cpst_dose",
    "stage2_util": "cpst_util",
    "carps_vae_alpha.py": "cpst_alpha.py",
    "alpha_model.py": "cpst_model.py",
    "carps-vae-receipts.sha256": "cpst-receipts.sha256",
    "carps-vae-verification.md": "cpst-verification.md",
}
_ORIGINAL_COMPONENTS = {current: original for original, current in _PATH_COMPONENTS.items()}


def _translate_path(path, components):
    return "/".join(components.get(part, part) for part in path.split("/"))


def current_path(path):
    """Translate complete path components in an original CPST identity."""
    return _translate_path(path, _PATH_COMPONENTS)


def historical_path(path):
    """Return the path component spellings used by the original Git revisions."""
    return _translate_path(path, _ORIGINAL_COMPONENTS)


_METADATA_NAMES = re.compile(
    rb"(?<![A-Za-z0-9_])(?:"
    + b"|".join(re.escape(name.encode()) for name in _PATH_COMPONENTS)
    + rb")(?![A-Za-z0-9_])"
)


def migrate_names(contents):
    """Translate current paths and model identity labels in textual release metadata."""
    contents = _METADATA_NAMES.sub(lambda match: current_path(match[0].decode()).encode(), contents)
    return contents.replace(b"carps-vae", b"cpst").replace(b"CarPS-VAE", b"CarPS-Transcriptome")


@lru_cache(maxsize=1)
def _registration_hashes():
    evidence = Path(__file__).with_name("registration-hashes.json")
    lookup = json.loads(evidence.read_bytes())
    if lookup["schema"] != "cpst-registration-hashes-v1":
        raise ValueError("unknown CPST registration hash schema")
    return lookup


def original_alpha_release_sha256():
    """The exact alpha file digest embedded in the accepted pre-rename gates."""
    return _registration_hashes()["alpha_original_release_sha256"]


def original_axis_manifest_sha256(manifest_sha256):
    """Resolve the path-only Denoise-DEG manifest reissue for sealed input checks."""
    return _registration_hashes()["denoise_manifest_sha256"].get(manifest_sha256, manifest_sha256)


def checkpoint_registration_sha256(registration_sha256):
    """Resolve a reissued registration to the original digest stored in its checkpoint."""
    return _registration_hashes()["registrations_sha256"].get(registration_sha256, registration_sha256)
