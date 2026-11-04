"""Native-axis corpus matrices: build, digest and write the four CM h5ads.

ADR 0033 decision 2/3: a corpus is the information-maximal substrate between a
raw download and a model's input.  One row per raw well, ``X`` = raw counts
summed per resolved HGNC id on the source's **native** gene axis, one
``log1p_cp10k`` convenience layer, the source's own row order, no pharmacology
filter, no Δ, no integer encodings.

The four CM corpora share one ``obs`` contract (``paths.CORPUS_OBS``): the 24
ADR-0026 columns in their v1 order plus ``cell_line``, ``n_raw_rows`` and
``is_blank``.  ``dataset_id`` uses the v2 short ids (``sctox``/``tenv``); the
``well`` string is kept **verbatim from the loader** so it still carries the
historical id prefix the v1 evidence used.  Downstream view data stores remap
``dataset_id``/``control_group`` (``sctox``→``screen_ctox``,
``tenv``→``temp_env``) during the G4 reconstruction (ticket #17); the corpus
does not rewrite provenance.

Equivalence is content-level, never file-byte (ADR 0033 decision 7): the digest
is the v1 recipe minus ``uns`` — ``X`` bytes, layer bytes, ``obs.to_csv`` and
``var.to_csv``.  The writer passes ``track_times=False`` and renames a temp file
atomically, so a rebuild is file-byte reproducible (gate G5).
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import write_h5ad as _write_h5ad

from . import paths
from .annotation import HgncIndex
from .drugs import CompoundTable, annotate_source_meta
from .sources import SourceData


class CorpusError(RuntimeError):
    """A corpus violated the native-axis contract."""


def _replicate_index(meta: pd.DataFrame) -> np.ndarray:
    """Within-group replicate index, deterministic by well id.

    Computed on the full source row set *before* the corpus row filter, exactly
    as the v1 view writer did (ADR 0026 decision 3), so the corpus stores the
    same values a filtered view would have derived.
    """
    group = np.where(
        meta["is_control"].to_numpy(),
        meta["control_group"].astype(str) + "|control|" + meta["vehicle"].astype(str),
        meta["control_group"].astype(str) + "|" + meta["drug_key"].astype(str)
        + "|" + meta["dose_um"].astype(str))
    order = np.lexsort((meta["well"].astype(str).to_numpy(), group))
    counts = np.zeros(len(meta), dtype=np.int32)
    last_group, running = None, 0
    for position in order:
        key = group[position]
        if key != last_group:
            running = 0
            last_group = key
        counts[position] = running
        running += 1
    return counts


def native_gene_axis(source: SourceData, index: HgncIndex
                     ) -> tuple[np.ndarray, pd.DataFrame]:
    """Sum raw counts per resolved HGNC id, in first-appearance order.

    The axis is the unique resolved HGNC ids of *this source* (first appearance
    over its raw gene ids), each represented by its canonical HGNC symbol.  Two
    raw ids (duplicate probes, alias symbols, Ensembl versions) that resolve to
    the same HGNC id are summed, which is exactly the integer the v1 writer
    produced before it projected onto the 18,164 CM axis.
    """
    mapping = source.gene_map
    position: dict[str, int] = {}
    hgnc_ids: list[str] = []
    match_types: list[str] = []
    for hgnc_id, match_type in zip(mapping.hgnc_ids, mapping.match_types):
        if hgnc_id is None:
            continue
        if hgnc_id not in position:
            position[hgnc_id] = len(hgnc_ids)
            hgnc_ids.append(hgnc_id)
            match_types.append(str(match_type))
    if not hgnc_ids:
        raise CorpusError(f"{source.dataset_id}: no resolved gene on the native axis")

    column = np.array([position.get(h, -1) for h in mapping.hgnc_ids],
                      dtype=np.int64)
    keep = column >= 0
    n_genes = len(hgnc_ids)
    aggregation = sp.csr_matrix(
        (np.ones(int(keep.sum()), dtype=np.int64),
         (np.arange(int(keep.sum())), column[keep])),
        shape=(int(keep.sum()), n_genes))

    # Chunk the row projection: the source matrix is int32 and small, but the
    # int64 matmul accumulators are the largest transient of the whole build.
    counts = np.zeros((source.counts.shape[0], n_genes), dtype=np.int32)
    limit = np.iinfo(np.int32).max
    step = 512
    for start in range(0, counts.shape[0], step):
        stop = min(start + step, counts.shape[0])
        block = np.asarray(source.counts[start:stop, keep], dtype=np.int64)
        summed = np.asarray(block @ aggregation)
        if summed.size and int(summed.max()) > limit:
            raise CorpusError(
                f"{source.dataset_id}: summed counts overflow int32 in rows "
                f"{start}..{stop}")
        counts[start:stop] = summed.astype(np.int32)

    symbols = [index.ref(hgnc_id).symbol for hgnc_id in hgnc_ids]
    if len(set(symbols)) != len(symbols):
        raise CorpusError(f"{source.dataset_id}: native axis has duplicate symbols")
    var = pd.DataFrame({
        "hgnc_id": hgnc_ids,
        "entrez_id": [index.ref(h).entrez_id for h in hgnc_ids],
        "ensembl_id": [index.ref(h).ensembl_id for h in hgnc_ids],
        "match_type": match_types,
    }, index=pd.Index(symbols, name="symbol"))
    return counts, var


def build_corpus(corpus_id: str, source: SourceData, table: CompoundTable,
                 index: HgncIndex) -> ad.AnnData:
    """Build one CM corpus in memory (nothing is written here)."""
    expected_source = paths.CORPUS_TO_SOURCE.get(corpus_id)
    if expected_source is None:
        raise CorpusError(f"unknown corpus {corpus_id!r}")
    if source.dataset_id != expected_source:
        raise CorpusError(f"{corpus_id}: source is {source.dataset_id!r}, "
                          f"expected {expected_source!r}")

    native_counts, var = native_gene_axis(source, index)

    # Row metadata exactly as the v1 writer derived it, then the corpus row
    # filter: controls + technical blanks + every SMILES-resolvable compound.
    meta = annotate_source_meta(source.meta.copy(), table)
    meta["control_group"] = (corpus_id + "::" + meta["plate"].astype(str))
    meta["replicate"] = _replicate_index(meta)
    keep = (meta["is_control"].to_numpy() | meta["is_blank"].to_numpy()
            | (meta["entity_type"] == "small_molecule").to_numpy())
    kept = meta.loc[keep]
    counts = np.ascontiguousarray(native_counts[keep])
    del native_counts

    lib_size = counts.sum(axis=1, dtype=np.int64)
    if (lib_size <= 0).any():
        bad = kept.loc[lib_size <= 0, "well"].tolist()[:5]
        raise CorpusError(f"{corpus_id}: zero library size for wells {bad}")
    n_detected = (counts > 0).sum(axis=1).astype(np.int32)

    obs = pd.DataFrame({
        "well": kept["well"].astype(str).values,
        "dataset_id": corpus_id,
        "compound_name": kept["compound_name"].astype(str).values,
        "drug_key": kept["drug_key"].fillna("").astype(str).values,
        "inchikey": kept["inchikey"].fillna("").astype(str).values,
        "smiles": kept["smiles"].fillna("").astype(str).values,
        "dose_um": kept["dose_um"].astype(float).values,
        "dose_source": kept["dose_source"].fillna("").astype(str).values,
        "time_h": kept["time_h"].astype(float).values,
        "is_control": kept["is_control"].astype(bool).values,
        "vehicle": kept["vehicle"].astype(str).values,
        "plate": kept["plate"].astype(str).values,
        "control_group": kept["control_group"].astype(str).values,
        "cell_type": kept["cell_type"].astype(str).values,
        "donor": kept["donor"].astype(str).values,
        "replicate": kept["replicate"].astype(np.int32).values,
        "atc": kept["atc"].fillna("").astype(str).values,
        "moa": kept["moa"].fillna("").astype(str).values,
        "mesh_classes": kept["mesh_classes"].fillna("").astype(str).values,
        "pharm_atc_moa": kept["pharm_atc_moa"].astype(bool).values,
        "pharm_mesh": kept["pharm_mesh"].astype(bool).values,
        "annotation_source": kept["annotation_source"].fillna("").astype(str).values,
        "lib_size": lib_size,
        "n_genes_detected": n_detected,
        "cell_line": kept["cell_line"].astype(str).values,
        "n_raw_rows": kept["n_raw_rows"].astype(np.int32).values,
        "is_blank": kept["is_blank"].astype(bool).values,
    }, index=pd.Index(kept["well"].astype(str).values, name="well"))
    if list(obs.columns) != paths.CORPUS_OBS:
        raise CorpusError(f"{corpus_id}: obs contract drifted: {list(obs.columns)}")

    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.layers["log1p_cp10k"] = np.log1p(
        counts / lib_size[:, None] * 1e4).astype(np.float32)
    return adata


def content_digest(adata: ad.AnnData) -> str:
    """v1 content recipe minus ``uns``: X bytes | layer bytes | obs | var."""
    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(adata.X).tobytes())
    digest.update(np.ascontiguousarray(adata.layers["log1p_cp10k"]).tobytes())
    digest.update(adata.obs.to_csv(index=False, lineterminator="\n").encode())
    digest.update(adata.var.to_csv(lineterminator="\n").encode())
    return digest.hexdigest()


def file_sha256(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def write_corpus(adata: ad.AnnData, path: Path) -> dict:
    """Write one corpus atomically (``.tmp`` then rename), byte-reproducibly."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    if temp.exists():
        temp.unlink()
    try:
        # `track_times=False` is what makes a rebuild file-byte reproducible
        # (gate G5); the module-level writer forwards it to every dataset.
        _write_h5ad(temp, adata, compression="gzip", track_times=False)
        digest = file_sha256(temp)
        size = temp.stat().st_size
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()
    return {
        "file": str(path.relative_to(paths.REPO_ROOT)),
        "file_sha256": digest,
        "bytes": size,
    }


def corpus_stats(adata: ad.AnnData, corpus_id: str) -> dict:
    """The counts a consumer contract cares about (G1/G2)."""
    obs = adata.obs
    control = obs["is_control"].to_numpy()
    blank = obs["is_blank"].to_numpy()
    perturbed = ~control & ~blank
    chemicals = int(obs.loc[perturbed, "drug_key"].astype(str).nunique())
    return {
        "rows": int(adata.n_obs),
        "perturbed": int(perturbed.sum()),
        "control": int(control.sum()),
        "blank": int(blank.sum()),
        "chemicals": chemicals,
        "axis": int(adata.n_vars),
        "dataset_id": corpus_id,
        "distinct_wells": int(obs["well"].astype(str).nunique()),
    }
