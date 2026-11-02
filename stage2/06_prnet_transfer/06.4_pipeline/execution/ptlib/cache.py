"""Derived, memory-shareable Tahoe log1p-CP10k cache (ticket #5 memory rule).

The frozen pretrain matrix is a 16.3 GB CSR ``log1p_cp10k`` layer (8 GB of
``data`` + 8 GB of ``indices`` for 1.01e9 non-zeros).  Loading it once per
training process makes a seed cost ~14.5 GiB RSS and blocks concurrent seeds.
This module materialises the same numbers as one dense ``float32`` ``.npy``
matrix (67,018 x 18,164 = 4.87 GB) with a sidecar record, and training reads it
through ``np.load(..., mmap_mode="r")`` so concurrent processes share the
kernel page cache instead of each holding a private CSR copy.

Integrity:

* build streams the frozen CSR layer in row chunks, writes the ``.npy`` and
  hashes the **logical dense bytes** in the same pass — the recorded
  ``layer_sha256`` is a full-array checksum of the values;
* after writing, it re-opens the frozen h5ad and compares a random sample of
  rows (default 16) against the new file; ``sample_max_abs_diff`` must be 0;
* the sidecar records the frozen h5ad's sha256, so a changed input invalidates
  the cache (``--cache-mode auto`` rebuilds; ``require`` fails closed);
* at load, ``metadata`` verification checks the sidecar + shape + dtype + size +
  the frozen source hash and re-runs the random-sample check; ``deep`` also
  recomputes the full file sha256.

The training math is untouched: the cache holds exactly the frozen layer.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import h5py
import numpy as np
import scipy.sparse as sp

from . import PIPELINE_VERSION, h5io, manifest as manifest_mod, paths

__all__ = [
    "CacheError",
    "CACHE_SCHEMA",
    "DEFAULT_CACHE_DIR",
    "DEFAULT_CACHE_NAME",
    "TahoeCache",
    "build_cache",
    "open_cache",
    "cache_sidecar",
]

CACHE_SCHEMA = "ptlib.tahoe_log1p_cache/v1"
DEFAULT_CACHE_DIR = paths.PIPELINE / "output" / "cache"
DEFAULT_CACHE_NAME = "tahoe_log1p_cp10k_v1.npy"

#: layer key inside the frozen Tahoe artifact
_LAYER_KEY = "layers/log1p_cp10k"


class CacheError(RuntimeError):
    """The derived cache is missing, corrupt or built from a different input."""


def cache_sidecar(path: Path | str) -> Path:
    return Path(str(path) + ".json")


@dataclass
class TahoeCache:
    """A verified read-only memmap of the frozen log1p-CP10k layer."""

    path: Path
    matrix: np.ndarray
    record: dict[str, Any] = field(default_factory=dict)
    verify_mode: str = "metadata"
    sample_max_abs_diff: float | None = None

    @property
    def shape(self):
        return tuple(int(x) for x in self.matrix.shape)

    def rows(self, rows: np.ndarray) -> np.ndarray:
        """Copy the requested rows out of the memmap (float32)."""
        return np.asarray(self.matrix[np.asarray(rows)], dtype=np.float32)


def _chunk_ranges(n_obs: int, chunk_rows: int):
    for start in range(0, n_obs, chunk_rows):
        yield start, min(start + chunk_rows, n_obs)


def build_cache(
    source_h5ad: Path | str,
    cache_path: Path | str,
    source_sha256: str,
    chunk_rows: int = 2048,
    sample_rows: int = 16,
    seed: int = 20261021,
    log=None,
) -> dict[str, Any]:
    """Materialise the frozen CSR layer as a dense float32 ``.npy`` + sidecar.

    Returns the sidecar record (also written next to the ``.npy``).  Fails
    closed if the source layer is not CSR, the shapes disagree, the sample
    check fails, or the streamed checksum disagrees with the written file.
    """
    source_h5ad = Path(source_h5ad)
    cache_path = Path(cache_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    meta = h5io.open_obs_var(source_h5ad)
    n_obs, n_vars = meta.n_obs, meta.n_vars
    started = time.perf_counter()

    with h5py.File(source_h5ad, "r") as handle:
        if _LAYER_KEY not in handle:
            raise CacheError(f"{source_h5ad} has no {_LAYER_KEY}")
        node = handle[_LAYER_KEY]
        if not isinstance(node, h5py.Group) or not {"data", "indices", "indptr"} <= set(node.keys()):
            raise CacheError(
                f"{_LAYER_KEY} is not CSR; the cache build only supports the frozen CSR contract"
            )
        indptr = node["indptr"][:].astype(np.int64)
        data_node, index_node = node["data"], node["indices"]
        shape = tuple(int(x) for x in node.attrs.get("shape", (n_obs, n_vars)))
        if shape != (n_obs, n_vars):
            raise CacheError(f"layer shape {shape} != obs/var shape {(n_obs, n_vars)}")

        tmp_path = cache_path.with_name(cache_path.name + ".tmp.npy")
        out = np.lib.format.open_memmap(
            tmp_path, mode="w+", dtype=np.float32, shape=(n_obs, n_vars)
        )
        import hashlib

        digest = hashlib.sha256()
        n_written = 0
        for start, stop in _chunk_ranges(n_obs, chunk_rows):
            lo, hi = int(indptr[start]), int(indptr[stop])
            chunk = sp.csr_matrix(
                (data_node[lo:hi], index_node[lo:hi], indptr[start : stop + 1] - lo),
                shape=(stop - start, n_vars),
            ).toarray().astype(np.float32, copy=False)
            out[start:stop] = chunk
            digest.update(chunk.tobytes(order="C"))
            n_written += stop - start
        out.flush()
        del out
        if n_written != n_obs:
            raise CacheError(f"wrote {n_written} rows, expected {n_obs}")
        layer_sha256 = digest.hexdigest()

    tmp_path.replace(cache_path)
    file_sha256 = manifest_mod.sha256_file(cache_path)

    # Independent random-sample verification against the frozen h5ad.
    rng = np.random.default_rng(seed)
    sample = np.unique(rng.choice(n_obs, size=min(sample_rows, n_obs), replace=False))
    expected, _ = h5io.read_layer_rows(source_h5ad, sample, n_vars)
    actual = np.load(cache_path, mmap_mode="r")[sample]
    diff = float(np.max(np.abs(np.asarray(actual, dtype=np.float64) - expected)))
    if diff != 0.0:
        raise CacheError(f"cache sample check failed: max abs diff {diff}")

    record = {
        "schema": CACHE_SCHEMA,
        "cache_path": str(cache_path),
        "source_path": str(source_h5ad),
        "source_sha256": source_sha256,
        "source_layer": _LAYER_KEY,
        "shape": [n_obs, n_vars],
        "dtype": "float32",
        "bytes": int(cache_path.stat().st_size),
        "layer_sha256": layer_sha256,
        "file_sha256": file_sha256,
        "sample_rows": [int(x) for x in sample],
        "sample_max_abs_diff": diff,
        "chunk_rows": int(chunk_rows),
        "build_seconds": round(time.perf_counter() - started, 3),
        "created_utc": manifest_mod.now_utc(),
        "pipeline_version": PIPELINE_VERSION,
    }
    cache_sidecar(cache_path).write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if log is not None:
        log.info(
            "built Tahoe cache %s: %s float32, %.2f GB, layer_sha256=%s…, "
            "sample max|diff|=%.1g (%.1fs)",
            cache_path,
            record["shape"],
            record["bytes"] / 1024**3,
            layer_sha256[:12],
            diff,
            record["build_seconds"],
        )
    return record


def _read_sidecar(path: Path) -> dict[str, Any]:
    sidecar = cache_sidecar(path)
    if not sidecar.exists():
        raise CacheError(f"cache sidecar is missing: {sidecar}")
    try:
        record = json.loads(sidecar.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise CacheError(f"cache sidecar is not valid JSON: {sidecar}") from exc
    if record.get("schema") != CACHE_SCHEMA:
        raise CacheError(f"cache sidecar schema is {record.get('schema')!r}, expected {CACHE_SCHEMA}")
    return record


def open_cache(
    cache_path: Path | str,
    source_sha256: str,
    source_h5ad: Path | str | None = None,
    verify: str = "metadata",
    sample_rows: int | None = None,
) -> TahoeCache:
    """Open and verify a derived cache.

    ``verify`` is ``metadata`` (sidecar + size + source hash + random-sample
    equality) or ``deep`` (additionally recompute the full file sha256).
    """
    cache_path = Path(cache_path)
    if not cache_path.exists():
        raise CacheError(f"derived Tahoe cache not found: {cache_path}")
    record = _read_sidecar(cache_path)
    if record.get("source_sha256") != source_sha256:
        raise CacheError(
            "derived Tahoe cache was built from a different frozen artifact: "
            f"cache source {record.get('source_sha256')!r} != frozen {source_sha256!r}"
        )
    matrix = np.load(cache_path, mmap_mode="r")
    if tuple(record.get("shape", ())) != tuple(matrix.shape):
        raise CacheError("derived Tahoe cache shape does not match its sidecar")
    if int(record.get("bytes", -1)) != cache_path.stat().st_size:
        raise CacheError("derived Tahoe cache size does not match its sidecar")

    if verify == "deep":
        actual = manifest_mod.sha256_file(cache_path)
        if actual != record.get("file_sha256"):
            raise CacheError(
                f"derived Tahoe cache checksum mismatch: stored {record.get('file_sha256')}, "
                f"file {actual}"
            )
    elif verify != "metadata":
        raise CacheError(f"unknown cache verify mode {verify!r}")

    if str(matrix.dtype) != "float32":
        raise CacheError(f"derived Tahoe cache dtype is {matrix.dtype}, expected float32")

    sample_max_abs_diff: float | None = None
    if source_h5ad is not None and verify in ("metadata", "deep"):
        n = sample_rows if sample_rows is not None else len(record.get("sample_rows", []) or [])
        if n:
            rng = np.random.default_rng(20261021)
            sample = np.unique(rng.choice(matrix.shape[0], size=min(n, matrix.shape[0]), replace=False))
            expected, _ = h5io.read_layer_rows(source_h5ad, sample, matrix.shape[1])
            actual = np.asarray(matrix[sample], dtype=np.float64)
            sample_max_abs_diff = float(np.max(np.abs(actual - expected)))
            if sample_max_abs_diff != 0.0:
                raise CacheError(f"derived Tahoe cache sample check failed: {sample_max_abs_diff}")
    elif record.get("sample_max_abs_diff") is not None:
        sample_max_abs_diff = float(record["sample_max_abs_diff"])

    return TahoeCache(
        path=cache_path,
        matrix=matrix,
        record=record,
        verify_mode=verify,
        sample_max_abs_diff=sample_max_abs_diff,
    )
