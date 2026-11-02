"""Table reading for the frozen input, with a derived CSV cache for the pinned env.

The pinned GPU env (`arc-stage2-bench-prnet`) has **no parquet engine** (no
pyarrow / fastparquet / polars) and no compression codecs, while ticket #3 ships
its fold / consensus-membership / control-group tables as parquet.  Rather than
mutate that pinned environment, this module reads parquet when an engine is
available and otherwise falls back to a **derived CSV cache** built once from the
same files by ``tools/convert_frozen_tables.py`` under the ``arc`` env.

The cache is a derived artifact, not frozen evidence:

* it lives in the pipeline's gitignored ``output/input_cache/``;
* each CSV records the source parquet's sha256 (the frozen manifest hash) and its
  own sha256; the loader refuses a cache entry whose source hash differs;
* the conversion is byte-preserving over values (parquet -> DataFrame -> CSV).

If neither an engine nor a valid cache exists, reading fails closed with the
exact command that fixes it.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from . import PIPELINE_VERSION, paths

__all__ = [
    "TableCacheError",
    "CACHE_SCHEMA",
    "DEFAULT_TABLE_CACHE",
    "read_table",
    "convert_parquet",
    "parquet_engine_available",
    "cache_manifest_path",
]

CACHE_SCHEMA = "ptlib.table_cache/v1"
DEFAULT_TABLE_CACHE = paths.PIPELINE / "output" / "input_cache"


class TableCacheError(RuntimeError):
    """A parquet table cannot be read and no valid derived cache exists."""


def parquet_engine_available() -> tuple[bool, str]:
    try:
        import pyarrow  # noqa: F401

        return True, "pyarrow"
    except Exception:
        pass
    try:
        import fastparquet  # noqa: F401

        return True, "fastparquet"
    except Exception:
        return False, "none"


def cache_manifest_path(cache_dir: Path | str) -> Path:
    return Path(cache_dir) / "manifest.json"


def _cache_record(cache_dir: Path, parquet_path: Path) -> dict[str, Any] | None:
    path = cache_manifest_path(cache_dir)
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    if payload.get("schema") != CACHE_SCHEMA:
        return None
    return (payload.get("tables") or {}).get(parquet_path.name)


def _read_cached_csv(cache_dir: Path, parquet_path: Path, source_sha256: str) -> pd.DataFrame:
    record = _cache_record(cache_dir, parquet_path)
    if record is None:
        raise TableCacheError(
            f"cannot read {parquet_path.name}: the pinned environment has no parquet engine "
            f"({parquet_engine_available()[1]}) and no derived CSV cache entry exists under "
            f"{cache_dir}. Build the cache once with the `arc` env:\n"
            f"    /home/chang/miniconda3/envs/arc/bin/python "
            f"stage2/06_prnet_transfer/06.4_pipeline/execution/tools/convert_frozen_tables.py "
            f"--data-dir <frozen input>"
        )
    if record.get("source_sha256") != source_sha256:
        raise TableCacheError(
            f"derived CSV cache for {parquet_path.name} was built from a different source parquet "
            f"(cache {record.get('source_sha256')!r} != frozen {source_sha256!r}); rebuild the cache"
        )
    csv_path = Path(cache_dir) / record["csv"]
    if not csv_path.exists():
        raise TableCacheError(f"derived CSV is missing: {csv_path}")
    frame = pd.read_csv(csv_path)
    return frame


def read_table(
    parquet_path: Path | str,
    expected_sha256: str | None = None,
    cache_dir: Path | str | None = DEFAULT_TABLE_CACHE,
) -> pd.DataFrame:
    """Read a frozen table from parquet, or from the derived CSV cache.

    ``expected_sha256`` should be the frozen manifest's hash for the parquet file;
    when given, it is used to validate the cache entry.
    """
    parquet_path = Path(parquet_path)
    if not parquet_path.exists():
        raise TableCacheError(f"frozen table is missing: {parquet_path}")
    available, engine = parquet_engine_available()
    if available:
        try:
            return pd.read_parquet(parquet_path)
        except (ImportError, ValueError) as exc:  # pragma: no cover - defensive
            if cache_dir is None:
                raise TableCacheError(
                    f"parquet engine {engine} is present but reading {parquet_path} failed: {exc}"
                ) from exc
    if cache_dir is None:
        raise TableCacheError(
            f"cannot read {parquet_path.name}: no parquet engine and no table cache directory"
        )
    return _read_cached_csv(Path(cache_dir), parquet_path, expected_sha256 or "")


def convert_parquet(parquet_path: Path | str, cache_dir: Path | str, source_sha256: str) -> dict[str, Any]:
    """Convert one parquet table to CSV in the derived cache; return its entry."""
    parquet_path = Path(parquet_path)
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_parquet(parquet_path)
    csv_name = f"{parquet_path.stem}.csv"
    csv_path = cache_dir / csv_name
    frame.to_csv(csv_path, index=False)
    import hashlib

    digest = hashlib.sha256()
    with open(csv_path, "rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    entry = {
        "source": parquet_path.name,
        "source_sha256": source_sha256,
        "csv": csv_name,
        "csv_sha256": digest.hexdigest(),
        "rows": int(len(frame)),
        "columns": [str(c) for c in frame.columns],
    }
    _update_cache_manifest(cache_dir, {parquet_path.name: entry})
    return entry


def _update_cache_manifest(cache_dir: Path, entries: dict[str, dict[str, Any]]) -> None:
    path = cache_manifest_path(cache_dir)
    payload: dict[str, Any] = {"schema": CACHE_SCHEMA, "pipeline_version": PIPELINE_VERSION, "tables": {}}
    if path.exists():
        try:
            existing = json.loads(path.read_text(encoding="utf-8"))
            if existing.get("schema") == CACHE_SCHEMA:
                payload["tables"] = existing.get("tables") or {}
        except json.JSONDecodeError:
            pass
    payload["tables"].update(entries)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
