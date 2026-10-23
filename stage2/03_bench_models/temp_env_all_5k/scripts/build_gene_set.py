#!/usr/bin/env python
"""Select the 5,000-gene axis for the `temp_env_all_5k` sensitivity benchmark.

The selection is the sum of per-source HVG ranks:

1. For each of the three CM views (`dtoxs_cm`, `screen_ctox`,
   `temp_env_all-smiles`), run Scanpy's Seurat HVG statistic on
   ``layers["log1p_cp10k"]`` over **all** wells, with no ``batch_key``.
2. Per source, genes whose mean expression lies in ``[0.0125, 3]`` are ranked
   ``1..M`` by descending ``dispersions_norm`` (ties by corpus order).  Every
   other gene -- including zero-mean genes with a NaN dispersion -- receives
   the worst rank ``n_genes + 1 = 18,165``.
3. Sum the three ranks and take the 5,000 smallest; ties are broken by corpus
   order.  The selected genes are written in corpus order.

The rule deliberately keeps low-expression genes out of the ranking rather
than letting NaN/zero dispersions masquerade as highly variable.  It is fixed
before any benchmark result is seen.

Usage
-----
    python stage2/03_bench_models/temp_env_all_5k/scripts/build_gene_set.py
    python stage2/03_bench_models/temp_env_all_5k/scripts/build_gene_set.py --check
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DATA = HERE.parent / "data"
CORPUS_GENES = ROOT / "stage2/02_cm_data/corpus/v1/genes.csv"
DATASET_DIR = ROOT / "stage2/02_cm_data/datasets-1"

SOURCES = ("dtoxs_cm", "screen_ctox", "temp_env_all-smiles")
TARGET_SIZE = 5_000
N_GENES = 18_164
WORST_RANK = N_GENES + 1
MIN_MEAN = 0.0125
MAX_MEAN = 3.0
FLAVOR = "seurat"
LAYER = "log1p_cp10k"

OUT_CSV = DATA / "gene_set_5k.csv"
OUT_RANKS = DATA / "gene_set_5k_ranks.csv"
OUT_MANIFEST = DATA / "gene_set_5k.manifest.json"


class GeneSetError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def corpus_symbols() -> np.ndarray:
    genes = pd.read_csv(CORPUS_GENES, usecols=["symbol"])["symbol"].astype(str)
    if len(genes) != N_GENES or genes.duplicated().any():
        raise GeneSetError("corpus/v1/genes.csv is not the frozen 18,164-gene table")
    return genes.to_numpy()


def _read_view(name: str, symbols: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
    """Return `(means, dispersions_norm, provenance)` for one source view."""

    path = DATASET_DIR / f"{name}.h5ad"
    if not path.exists():
        raise GeneSetError(f"missing source view: {path}")
    adata = ad.read_h5ad(path)
    if tuple(adata.var_names.astype(str)) != tuple(symbols):
        raise GeneSetError(f"{name}: var order does not match corpus/v1/genes.csv")
    if LAYER not in adata.layers:
        raise GeneSetError(f"{name}: missing layers[{LAYER!r}]")
    if tuple(adata.shape[1:]) != (N_GENES,):
        raise GeneSetError(f"{name}: unexpected shape {adata.shape}")

    # Scanpy writes the statistic into a fresh frame; the AnnData is discarded
    # immediately so only one view is resident at a time.
    import scanpy as sc

    stats = sc.pp.highly_variable_genes(
        adata,
        layer=LAYER,
        flavor=FLAVOR,
        n_top_genes=None,
        inplace=False,
    )
    means = stats["means"].to_numpy(dtype=np.float64)
    dispersion = stats["dispersions_norm"].to_numpy(dtype=np.float64)
    if means.shape != (N_GENES,) or dispersion.shape != (N_GENES,):
        raise GeneSetError(f"{name}: unexpected HVG frame {stats.shape}")

    info = {
        "name": name,
        "path": str(path.relative_to(ROOT)),
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "n_cells_used": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
    }
    del adata
    return means, dispersion, info


def rank_view(
    means: np.ndarray, dispersion: np.ndarray, corpus_index: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Rank one source under the fixed rule; non-evaluable genes tie last."""

    eligible = (
        np.isfinite(means)
        & (means >= MIN_MEAN)
        & (means <= MAX_MEAN)
        & np.isfinite(dispersion)
    )
    rank = np.full(N_GENES, WORST_RANK, dtype=np.int64)
    order = np.flatnonzero(eligible)
    # Descending dispersion; corpus order breaks ties (stable, explicit key).
    order = order[np.lexsort((corpus_index[order], -dispersion[order]))]
    rank[order] = np.arange(1, len(order) + 1, dtype=np.int64)
    return rank, eligible


def select(
    stats: dict[str, tuple[np.ndarray, np.ndarray]], symbols: np.ndarray
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return `(full ranked table, selected rows in corpus order)`."""

    corpus_index = np.arange(N_GENES, dtype=np.int64)
    frame = pd.DataFrame({"symbol": symbols, "corpus_index": corpus_index})
    rank_sum = np.zeros(N_GENES, dtype=np.int64)

    for name in SOURCES:
        means, dispersion = stats[name]
        rank, eligible = rank_view(means, dispersion, corpus_index)
        rank_sum += rank
        frame[f"{name}__mean"] = means
        frame[f"{name}__dispersions_norm"] = dispersion
        frame[f"{name}__eligible"] = eligible
        frame[f"{name}__rank"] = rank

    frame["rank_sum"] = rank_sum

    # Pick the 5,000 smallest rank sums; corpus order breaks ties.  Then store
    # the result in corpus order so the benchmark gene axis is a stable subset.
    order = np.lexsort((corpus_index, rank_sum))
    selected = np.sort(order[:TARGET_SIZE])
    out = frame.iloc[selected].reset_index(drop=True)

    # Fail closed if the persisted rows are not exactly that boundary.
    if not np.array_equal(selected, out["corpus_index"].to_numpy()):
        raise GeneSetError("selected rows are not in corpus order")
    if not np.array_equal(
        np.sort(order[:TARGET_SIZE]), out["corpus_index"].to_numpy()
    ):
        raise GeneSetError("selected rows are not the 5,000 lowest rank sums")
    return frame, out


def build() -> dict:
    symbols = corpus_symbols()
    stats = {}
    inputs = []
    for name in SOURCES:
        means, dispersion, info = _read_view(name, symbols)
        stats[name] = (means, dispersion)
        inputs.append(info)

    full, frame = select(stats, symbols)
    if len(frame) != TARGET_SIZE or frame["symbol"].duplicated().any():
        raise GeneSetError(f"selection has {len(frame)} rows / duplicates")

    DATA.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT_CSV, index=False, lineterminator="\n")
    full.to_csv(OUT_RANKS, index=False, lineterminator="\n")
    csv_hash = sha256_file(OUT_CSV)
    ranks_hash = sha256_file(OUT_RANKS)

    from importlib.metadata import version as _pkg_version

    manifest = {
        "tool": "scripts/build_gene_set.py",
        "inputs": inputs,
        "corpus": {
            "path": str(CORPUS_GENES.relative_to(ROOT)),
            "sha256": sha256_file(CORPUS_GENES),
            "n_genes": N_GENES,
        },
        "method": {
            "scanpy_version": _pkg_version("scanpy"),
            "anndata_version": _pkg_version("anndata"),
            "numpy_version": np.__version__,
            "pandas_version": pd.__version__,
            "function": "scanpy.pp.highly_variable_genes",
            "layer": LAYER,
            "flavor": FLAVOR,
            "n_top_genes": None,
            "batch_key": None,
            "cells": "all wells of each view (treated + controls)",
            "statistic": "dispersions_norm",
            "min_mean": MIN_MEAN,
            "max_mean": MAX_MEAN,
        },
        "rank_rule": (
            "per source: genes with finite mean in [0.0125, 3] ranked 1..M by "
            "descending dispersions_norm, ties by corpus order; all other "
            "genes (including zero-mean/NaN) rank 18165; sum the three ranks; "
            "take the 5,000 smallest, ties by corpus order; output in corpus "
            "order"
        ),
        "target_size": TARGET_SIZE,
        "n_genes": N_GENES,
        "worst_rank": WORST_RANK,
        "selection": {
            "path": str(OUT_CSV.relative_to(ROOT)),
            "sha256": csv_hash,
            "bytes": OUT_CSV.stat().st_size,
            "n_genes": int(len(frame)),
            "rank_sum_max": int(frame["rank_sum"].max()),
            "rank_sum_min": int(frame["rank_sum"].min()),
            "all_ranks": {
                "path": str(OUT_RANKS.relative_to(ROOT)),
                "sha256": ranks_hash,
                "bytes": OUT_RANKS.stat().st_size,
                "rows": int(len(full)),
            },
            "n_eligible_genes": {
                name: int(full[f"{name}__eligible"].sum()) for name in SOURCES
            },
            "n_low_expression": {
                name: int((~frame[f"{name}__eligible"]).sum()) for name in SOURCES
            },
        },
    }
    OUT_MANIFEST.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def check() -> None:
    if not OUT_CSV.exists() or not OUT_MANIFEST.exists():
        raise GeneSetError("run build_gene_set.py before --check")
    manifest = json.loads(OUT_MANIFEST.read_text())
    if sha256_file(OUT_CSV) != manifest["selection"]["sha256"]:
        raise GeneSetError("gene_set_5k.csv hash does not match its manifest")
    ranks = manifest["selection"]["all_ranks"]
    if sha256_file(ROOT / ranks["path"]) != ranks["sha256"]:
        raise GeneSetError("gene_set_5k_ranks.csv hash does not match its manifest")
    frame = pd.read_csv(OUT_CSV)
    if len(frame) != TARGET_SIZE or frame["symbol"].duplicated().any():
        raise GeneSetError("gene set is not 5,000 unique symbols")
    if not frame["corpus_index"].is_monotonic_increasing:
        raise GeneSetError("gene set is not in corpus order")
    if manifest["method"]["batch_key"] is not None:
        raise GeneSetError("manifest records a batch_key")
    print(f"OK: {len(frame)} genes, rank_sum {frame['rank_sum'].min()}..{frame['rank_sum'].max()}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    if args.check:
        check()
        return 0
    manifest = build()
    print(json.dumps(manifest["selection"], indent=1))
    print(f"wrote {OUT_CSV.relative_to(ROOT)} and {OUT_MANIFEST.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
