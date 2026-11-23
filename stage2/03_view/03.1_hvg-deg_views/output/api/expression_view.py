"""Bounded reads of the published raw-first HVG-DEG expression views."""
from __future__ import annotations

import json
import operator
from hashlib import file_digest, sha256 as hash_bytes
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

SCHEMA = "stage2-hvg-deg-views-raw-first-v1"
DEFAULT = Path(__file__).resolve().parents[1] / "raw-first-v1"
REPO = next(parent for parent in DEFAULT.parents if (parent / ".git").exists())
CORPUS_MANIFEST = REPO / "stage2/01_corpus/output/native/manifest.json"
CANONICAL_AXIS = REPO / "stage2/02_gene_axis/02.1_hvg-deg/output/release/gene_axis_5000.csv"
NORMALIZATION = "sum selected-axis integer counts, then separate CP10K and log1p"
CONTROL = "union matched vehicle rows by source observation identity"
DELTA = "X minus control, or the corresponding merge delta for pooled labels"
SOLO_SEED = 20261122
RECIPE = {"default_family": "full", "normalization": NORMALIZATION, "control": CONTROL,
          "delta": DELTA, "solo_seed": SOLO_SEED}
MAX_ROWS = 256
ARMS = {
    "tenv": ("full", "merge", "solo", "solo-deepest", "full-pooled-label"),
    "tahoe_bulk": ("full", "merge", "solo", "full-pooled-label"),
}
LAYERS = ("X", "control", "delta")


def sha256(path):
    with Path(path).open("rb") as handle:
        return file_digest(handle, "sha256").hexdigest()


def layer_path(layer):
    return layer if layer == "X" else f"layers/{layer}"


class View:
    """One source and view family. The full family is the default expression release."""

    def __init__(self, root, source, family, corpus_sha256, axis_sha256):
        if source not in ARMS or family not in ARMS[source]:
            raise ValueError(f"unknown expression view arm: {source}/{family}")
        base = root / "manifest.json"
        payload = base.read_bytes()
        manifest = json.loads(payload)
        if manifest.get("schema") != SCHEMA or (
            corpus_sha256 is not None and manifest["corpus_sha256"] != corpus_sha256
        ) or (axis_sha256 is not None and manifest["axis_sha256"] != axis_sha256):
            raise ValueError("expression view generation does not match the requested corpus or axis")
        if any(manifest.get(key) != value for key, value in RECIPE.items()):
            raise ValueError("expression view recipe differs from the accepted raw-first recipe")
        self.manifest = manifest
        self.manifest_sha256 = hash_bytes(payload).hexdigest()
        self.source = source
        self.family = family
        self.root = root
        self._handles = {}
        try:
            key = f"{family}/{source}.h5ad"
            handle = self._open(key, family)
            self.obs = read_elem(handle["obs"])
            self.var = read_elem(handle["var"])
            self.shape = (len(self.obs), len(self.var))
            if (self.obs.index.has_duplicates or self.var.index.has_duplicates
                    or self.shape[1] != manifest["axis_size"]):
                raise ValueError("expression view identities or axis size differ")
            self.exclusions = pd.read_csv(self._checked(f"{source}.exclusions.csv"),
                                          keep_default_na=False)
            self._arrays = {
                layer: handle[layer_path(layer)] for layer in LAYERS
                if layer_path(layer) in handle
            }
            if family == "full-pooled-label":
                full = self._open(f"full/{source}.h5ad", "full")
                if not self.obs.equals(read_elem(full["obs"])) or not self.var.equals(read_elem(full["var"])):
                    raise ValueError("pooled-label view does not bind identical full-view rows and genes")
                self._arrays.update({layer: full[layer_path(layer)] for layer in ("X", "control")})
            if set(self._arrays) != set(LAYERS) or any(
                block.shape != self.shape or block.dtype != np.float32 for block in self._arrays.values()
            ):
                raise ValueError("expression view layers have unsupported shapes or dtypes")
        except BaseException:
            self.close()
            raise

    def _checked(self, relative):
        path = self.root / relative
        if sha256(path) != self.manifest["outputs"][relative]:
            raise ValueError(f"expression view payload digest differs: {relative}")
        return path

    def _open(self, relative, family):
        handle = h5py.File(self._checked(relative), "r")
        self._handles[relative] = handle
        for key in ("schema", "corpus_sha256", "axis_sha256"):
            if handle.attrs[key] != self.manifest[key]:
                raise ValueError(f"mixed expression view generation: {relative}")
        if handle.attrs["source"] != self.source or handle.attrs["family"] != family:
            raise ValueError(f"expression view arm differs: {relative}")
        return handle

    def read(self, start, stop, layer="X"):
        """Read at most 256 consecutive rows, using X, control, or delta."""
        start, stop = operator.index(start), operator.index(stop)
        if not 0 <= start <= stop <= self.shape[0] or stop - start > MAX_ROWS:
            raise ValueError(f"expression window must contain at most {MAX_ROWS} valid rows")
        if layer not in LAYERS:
            raise ValueError(f"unknown expression view layer: {layer}")
        return self._arrays[layer][start:stop]

    def batches(self, layer="X", batch_size=128):
        """Yield row position and a bounded block for the selected layer."""
        batch_size = operator.index(batch_size)
        if not 1 <= batch_size <= MAX_ROWS:
            raise ValueError(f"batch size must be within 1 to {MAX_ROWS}")
        for start in range(0, self.shape[0], batch_size):
            yield start, self.read(start, min(start + batch_size, self.shape[0]), layer)

    def close(self):
        for handle in self._handles.values():
            handle.close()
        self._handles.clear()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def open_view(path=None, *, source, family="full", corpus_sha256=None, axis_sha256=None):
    root = DEFAULT if path is None else Path(path)
    if root.resolve() == DEFAULT.resolve():
        corpus = CORPUS_MANIFEST.read_bytes()
        if json.loads(corpus).get("schema") != "stage2-native-corpus-v1":
            raise ValueError("canonical native corpus provenance is missing")
        axis = sha256(CANONICAL_AXIS)
        axis_manifest = json.loads(CANONICAL_AXIS.with_name("manifest.json").read_bytes())
        if axis_manifest["axes"]["5000"] != axis:
            raise ValueError("canonical prediction axis payload digest differs")
        canonical = (hash_bytes(corpus).hexdigest(), axis)
        specified = (corpus_sha256, axis_sha256)
        if any(pin is not None and pin != current for pin, current in zip(specified, canonical, strict=True)):
            raise ValueError("requested corpus or axis does not match the canonical release inputs")
        corpus_sha256, axis_sha256 = canonical
    return View(root, source, family, corpus_sha256, axis_sha256)
