"""Minimal HDF5 access for the frozen artifacts.

The frozen CM matrix is a dense ``int32`` ``X`` plus a dense ``float32``
``layers['log1p_cp10k']``; the Tahoe matrix is CSR.  Nothing here densifies a
whole matrix by accident:

* ``open_obs_var`` decodes only the ``obs``/``var`` element groups with h5py —
  never through ``anndata.read_h5ad(..., backed="r")``, whose "backed" mode
  still materialises every ``layers`` entry (8.1 GB on the frozen Tahoe view,
  12.9 GB on the 25.8 GB corpus) before the metadata copy;
* ``read_dense_rows`` reads a requested row subset in contiguous runs (HDF5
  fancy indexing needs sorted, unique indices — the run split is what makes it
  fast);
* ``read_csr_rows`` slices the CSR ``data``/``indices`` arrays by row ranges,
  so a 16 GB pretrain matrix can be sampled without materialising it;
* ``read_csr_full`` is the explicit "load the whole layer" path the formal
  pretrain uses (the same 14.5 GiB RSS pattern as ``05/train_prnet.py``).

All readers fail closed on a missing key, a shape mismatch or an unsupported
encoding; the caller never silently falls back to a different matrix.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import h5py
import numpy as np
import scipy.sparse as sp

__all__ = [
    "H5Error",
    "BackedMeta",
    "open_obs_var",
    "close_backed",
    "detect_encoding",
    "read_dense_rows",
    "read_csr_rows",
    "read_csr_full",
    "read_layer_rows",
    "read_layer_all",
]


class H5Error(RuntimeError):
    """A frozen artifact does not match the contracted encoding."""


@dataclass(frozen=True)
class BackedMeta:
    n_obs: int
    n_vars: int
    obs: Any  # pandas.DataFrame
    var: Any  # pandas.DataFrame
    uns_keys: tuple[str, ...]


def _read_elem(payload: Any) -> Any:
    """``anndata.io.read_elem`` on 0.11+, ``anndata.experimental.read_elem`` before."""
    try:  # anndata >= 0.11
        from anndata.io import read_elem
    except ImportError:  # anndata 0.10.x
        from anndata.experimental import read_elem
    return read_elem(payload)


def _matrix_shape(handle: h5py.File) -> tuple[int, int] | None:
    """``(n_obs, n_vars)`` from the ``X`` node header, without reading values."""
    if "X" not in handle:
        return None
    node = handle["X"]
    if isinstance(node, h5py.Group):
        raw = node.attrs.get("shape")
        if raw is None:  # pragma: no cover - defensive
            return None
        return (int(raw[0]), int(raw[1]))
    if isinstance(node, h5py.Dataset) and len(node.shape) == 2:
        return (int(node.shape[0]), int(node.shape[1]))
    return None  # pragma: no cover - defensive


def open_obs_var(path: Path | str) -> BackedMeta:
    """Decode ``obs`` / ``var`` / ``uns`` keys from an h5ad header.

    Reads the element groups with h5py instead of ``anndata.read_h5ad(backed="r")``:
    backed mode keeps ``X`` lazy but **materialises every layer**, which is 8.1 GB
    for the frozen Tahoe view and 12.9 GB for the CM corpus — more than the WSL2
    budget before training starts (see ``docs/agents/execution-environment.md``).
    Here the peak stays proportional to the metadata (a few MB for 67k rows), and
    ``X`` / ``layers`` are never touched.  The expression values have explicit
    readers (:func:`read_layer_rows`, the ``cache`` memmap) or nothing.
    """
    path = Path(path)
    if not path.exists():
        raise H5Error(f"frozen artifact not found: {path}")
    try:
        with h5py.File(path, "r") as handle:
            obs = _read_elem(handle["obs"]) if "obs" in handle else None
            var = _read_elem(handle["var"]) if "var" in handle else None
            uns_keys = (
                tuple(sorted(str(k) for k in handle["uns"].keys())) if "uns" in handle else ()
            )
            shape = _matrix_shape(handle)
    except (OSError, KeyError) as exc:
        raise H5Error(f"cannot read obs/var from {path}: {exc}") from exc
    if obs is None or var is None:
        raise H5Error(f"{path}: h5ad has no obs/var group")
    n_obs, n_vars = int(len(obs)), int(len(var))
    if shape is not None:
        if shape != (n_obs, n_vars):
            raise H5Error(
                f"{path}: obs/var lengths {(n_obs, n_vars)} disagree with X shape {shape}"
            )
        n_obs, n_vars = shape
    return BackedMeta(n_obs, n_vars, obs, var, uns_keys)


def close_backed(path: Path | str | None = None) -> None:
    """Deprecated no-op kept for API compatibility.

    :func:`open_obs_var` no longer holds a backed handle (a backed open is what
    materialised the layers), so there is nothing to close.
    """
    del path


def _resolve_key(handle: h5py.File, candidates: Sequence[str]) -> str:
    for candidate in candidates:
        if candidate in handle:
            return candidate
    raise H5Error(
        f"none of the expected HDF5 keys is present: {list(candidates)} "
        f"(file has {sorted(handle.keys())})"
    )


def detect_encoding(path: Path | str, key: str) -> str:
    """Return ``"csr"`` or ``"dense"`` for an HDF5 matrix node."""
    with h5py.File(path, "r") as handle:
        node = handle[_resolve_key(handle, [key])]
        if isinstance(node, h5py.Group):
            if {"data", "indices", "indptr"} <= set(node.keys()):
                return "csr"
            raise H5Error(f"HDF5 group {key!r} is not a CSR matrix (keys: {sorted(node.keys())})")
        if isinstance(node, h5py.Dataset):
            return "dense"
        raise H5Error(f"unsupported HDF5 node for {key!r}: {type(node).__name__}")


def _contiguous_runs(rows: np.ndarray) -> list[tuple[int, int, int]]:
    """Split sorted unique *rows* into ``(start, stop, row_pos)`` contiguous runs."""
    rows = np.asarray(rows, dtype=np.int64)
    runs: list[tuple[int, int, int]] = []
    start_pos = 0
    for pos in range(1, len(rows) + 1):
        if pos == len(rows) or rows[pos] != rows[pos - 1] + 1:
            runs.append((int(rows[start_pos]), int(rows[pos - 1]) + 1, start_pos))
            start_pos = pos
    return runs


def _check_rows(rows: np.ndarray, n_obs: int) -> np.ndarray:
    rows = np.asarray(rows)
    if rows.ndim != 1:
        raise H5Error(f"row indices must be 1-D, got shape {rows.shape}")
    if rows.size == 0:
        return rows.astype(np.int64)
    if rows.min() < 0 or rows.max() >= n_obs:
        raise H5Error(f"row index outside 0..{n_obs - 1}: {rows.min()}..{rows.max()}")
    unique = np.unique(rows.astype(np.int64))
    if unique.size != rows.size:
        raise H5Error("row indices contain duplicates")
    return unique


def read_dense_rows(
    path: Path | str, key: str, rows: np.ndarray, n_vars: int | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Read a dense matrix row subset.

    Returns ``(values, sorted_rows)`` where ``values[i]`` is the row
    ``sorted_rows[i]`` in ascending order.
    """
    rows = _check_rows(np.asarray(rows), n_obs=2**62)  # bounds checked below
    with h5py.File(path, "r") as handle:
        node = handle[_resolve_key(handle, [key])]
        if isinstance(node, h5py.Group):
            raise H5Error(f"{key!r} is not a dense dataset")
        shape = tuple(int(x) for x in node.shape)
        if len(shape) != 2:
            raise H5Error(f"expected a 2-D matrix at {key!r}, got shape {shape}")
        if n_vars is not None and shape[1] != int(n_vars):
            raise H5Error(f"{key!r} has {shape[1]} columns, expected {n_vars}")
        if rows.size and rows.max() >= shape[0]:
            raise H5Error(f"row index {rows.max()} outside {key!r} with {shape[0]} rows")
        pieces = [node[r0:r1] for r0, r1, _ in _contiguous_runs(rows)]
    if not pieces:
        return np.empty((0, shape[1]), dtype=np.float32), rows
    return np.vstack(pieces).astype(np.float32, copy=False), rows


def read_csr_rows(
    path: Path | str, key: str, rows: np.ndarray, n_vars: int | None = None
) -> tuple[sp.csr_matrix, np.ndarray]:
    """Read a CSR matrix row subset without touching the rest of the matrix.

    Returns ``(matrix, sorted_rows)``; ``matrix[i]`` is row ``sorted_rows[i]``.
    """
    rows = _check_rows(np.asarray(rows), n_obs=2**62)
    with h5py.File(path, "r") as handle:
        node = handle[_resolve_key(handle, [key])]
        if not isinstance(node, h5py.Group) or not {"data", "indices", "indptr"} <= set(node.keys()):
            raise H5Error(f"{key!r} is not a CSR matrix")
        indptr = node["indptr"][:].astype(np.int64)
        shape = tuple(int(x) for x in node.attrs.get("shape", (indptr.size - 1, n_vars)))
        if len(shape) != 2 or shape[0] != indptr.size - 1:
            raise H5Error(f"inconsistent CSR shape {shape} vs indptr {indptr.size - 1}")
        if n_vars is not None and shape[1] != int(n_vars):
            raise H5Error(f"{key!r} has {shape[1]} columns, expected {n_vars}")
        if rows.size and rows.max() >= shape[0]:
            raise H5Error(f"row index {rows.max()} outside {key!r} with {shape[0]} rows")
        data_node, index_node = node["data"], node["indices"]
        data_chunks: list[np.ndarray] = []
        index_chunks: list[np.ndarray] = []
        local_indptr = [0]
        total = 0
        for r0, r1, _ in _contiguous_runs(rows):
            start, stop = int(indptr[r0]), int(indptr[r1])
            data_chunks.append(data_node[start:stop])
            index_chunks.append(index_node[start:stop])
            local_indptr.extend((indptr[r0 + 1 : r1 + 1] - start + total).tolist())
            total += stop - start
    if not data_chunks:
        return sp.csr_matrix((0, shape[1]), dtype=np.float32), rows
    data = np.concatenate(data_chunks)
    indices = np.concatenate(index_chunks)
    matrix = sp.csr_matrix(
        (data.astype(np.float32, copy=False), indices, np.asarray(local_indptr, dtype=np.int64)),
        shape=(rows.size, shape[1]),
    )
    return matrix, rows


def read_csr_full(path: Path | str, key: str, n_obs: int, n_vars: int) -> sp.csr_matrix:
    """Load a whole CSR layer (the formal-pretrain path; ~8 GB for 67k x 18,164)."""
    with h5py.File(path, "r") as handle:
        node = handle[_resolve_key(handle, [key])]
        if not isinstance(node, h5py.Group) or not {"data", "indices", "indptr"} <= set(node.keys()):
            raise H5Error(f"{key!r} is not a CSR matrix")
        data = node["data"][:]
        indices = node["indices"][:]
        indptr = node["indptr"][:]
        shape = tuple(int(x) for x in node.attrs.get("shape", (n_obs, n_vars)))
    if shape != (n_obs, n_vars):
        raise H5Error(f"{key!r} shape {shape} != expected {(n_obs, n_vars)}")
    matrix = sp.csr_matrix((data.astype(np.float32, copy=False), indices, indptr), shape=shape)
    if int(matrix.indptr[0]) != 0 or int(matrix.indptr[-1]) != matrix.data.size:
        raise H5Error(f"{key!r} indptr is inconsistent")
    return matrix


#: layer names the frozen input may use for the log1p-CP10k expression matrix.
LOG1P_KEYS = (
    "layers/log1p_cp10k",
    "layers/log1p_cp10k_layer",
    "layers/log1p_cpm10k",
)


def read_layer_rows(
    path: Path | str, rows: np.ndarray, n_vars: int | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Read the log1p-CP10k row subset, regardless of dense/CSR storage."""
    try:
        encoding = detect_encoding(path, LOG1P_KEYS[0])
    except H5Error:
        # Try the alternative keys before failing.
        for candidate in LOG1P_KEYS[1:]:
            try:
                encoding = detect_encoding(path, candidate)
            except H5Error:
                continue
            else:
                return _read_layer(path, candidate, encoding, rows, n_vars)
        raise H5Error(
            f"{path} has no recognised log1p-CP10k layer (looked for {list(LOG1P_KEYS)})"
        )
    return _read_layer(path, LOG1P_KEYS[0], encoding, rows, n_vars)


def _read_layer(
    path: Path | str, key: str, encoding: str, rows: np.ndarray, n_vars: int | None
) -> tuple[np.ndarray, np.ndarray]:
    if encoding == "csr":
        matrix, sorted_rows = read_csr_rows(path, key, rows, n_vars)
        return matrix.toarray().astype(np.float32, copy=False), sorted_rows
    return read_dense_rows(path, key, rows, n_vars)


def read_layer_all(path: Path | str, n_obs: int, n_vars: int):
    """Load the whole log1p-CP10k layer (CSR stays sparse; dense stays dense)."""
    for candidate in LOG1P_KEYS:
        try:
            encoding = detect_encoding(path, candidate)
        except H5Error:
            continue
        if encoding == "csr":
            return read_csr_full(path, candidate, n_obs, n_vars)
        with h5py.File(path, "r") as handle:
            node = handle[candidate]
            shape = tuple(int(x) for x in node.shape)
            if shape != (n_obs, n_vars):
                raise H5Error(f"{candidate} shape {shape} != {(n_obs, n_vars)}")
            return node[:]
    raise H5Error(f"{path} has no recognised log1p-CP10k layer (looked for {list(LOG1P_KEYS)})")
