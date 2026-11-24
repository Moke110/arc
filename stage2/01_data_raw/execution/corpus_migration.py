"""Migrate a sealed flat corpus and independently check unchanged datasets."""
from __future__ import annotations

import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from cardiac_sources import merge_primary_tables
from native_corpus import CorpusError, LEGACY_SCHEMA, MAX_ROWS, open_corpus

UNCHANGED = {
    "tenv": "tenv", "tahoe_bulk": "tahoe_bulk", "sctox": "screen_ctox",
    "dtoxs_ipsc_conv": "dtoxs_ipsc_conv", "dtoxs_ipsc_dge": "dtoxs_ipsc_dge",
    "gse247541": "gse247541",
}


def predecessor_record(corpus):
    if corpus.manifest["schema"] != LEGACY_SCHEMA:
        raise CorpusError("the corpus predecessor must be a sealed v1 release")
    return {"manifest_sha256": corpus.manifest_sha256, "unchanged_sources": UNCHANGED}


def stage_predecessor(corpus, candidate, axis):
    """Stage cached counts only after checking the accepted predecessor outputs."""
    predecessor_record(corpus)
    columns = {}
    for name, previous in UNCHANGED.items():
        record = corpus.manifest["sources"][previous]
        shutil.copyfile(corpus.root / f"{previous}.h5ad", candidate / f"{name}.h5ad")
        if "observation_sha256" in record:
            shutil.copyfile(corpus.root / f"{previous}.obs.h5", candidate / f"{name}.obs.h5")
        columns[name] = record["common_columns"]
    for extra in ("tenv_manifest.json", "tenv_probe_mapping.csv"):
        shutil.copyfile(corpus.root / extra, candidate / extra)

    prefix = "LINCS-LDS-1593/source/"
    paths = sorted(path.removeprefix(prefix)
                   for path in corpus.manifest["inputs"]["additional_count_files"]
                   if path.startswith(prefix) and "Read-Counts" in Path(path).name)
    tables = {path: corpus.root / f"dtoxs_primary_{index:02}.h5ad"
              for index, path in enumerate(paths, 1)}
    if len(tables) != 8:
        raise CorpusError("the predecessor must contain eight DToxS-primary tables")
    columns["dtoxs_primary"] = merge_primary_tables(candidate, tables, axis)
    return columns


def _check_dense(before, after, name):
    if before.shape != after.shape or before.dtype != after.dtype:
        raise CorpusError(f"{name}: unchanged count matrix layout differs")
    batch = min(MAX_ROWS, max(1, 8_000_000 // (before.dtype.itemsize * before.shape[1])))
    for start in range(0, before.shape[0], batch):
        if not np.array_equal(before[start:start + batch], after[start:start + batch]):
            raise CorpusError(f"{name}: unchanged count values differ at row {start}")


def verify_equivalence(predecessor, candidate):
    """Verify native counts, order, control and source metadata in bounded windows."""
    with open_corpus(candidate) as successor:
        if successor.manifest.get("predecessor") != predecessor_record(predecessor):
            raise CorpusError("candidate predecessor provenance differs")
        if successor.manifest["common_axis"] != predecessor.manifest["common_axis"]:
            raise CorpusError("candidate common axis differs from predecessor")
        for name, previous in UNCHANGED.items():
            before = predecessor.source(previous)
            after = successor.source(name)
            pd.testing.assert_frame_equal(before.obs, after.obs)
            pd.testing.assert_frame_equal(before.var, after.var)
            if (predecessor.manifest["sources"][previous]["common_columns"]
                    != successor.manifest["sources"][name]["common_columns"]):
                raise CorpusError(f"{name}: candidate gene map differs")
            with (h5py.File(predecessor.root / f"{previous}.h5ad") as old,
                  h5py.File(successor.root / name / f"{name}.h5ad") as new):
                left, right = old["X"], new["X"]
                if isinstance(left, h5py.Dataset) and isinstance(right, h5py.Dataset):
                    _check_dense(left, right, name)
                else:
                    for (start, rows), (new_start, new_rows) in zip(before.batches(), after.batches(), strict=True):
                        if start != new_start or not np.array_equal(rows, new_rows):
                            raise CorpusError(f"{name}: unchanged count values differ at row {start}")
