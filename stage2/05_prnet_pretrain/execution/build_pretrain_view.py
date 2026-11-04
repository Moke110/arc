#!/usr/bin/env python
"""Build the model-facing pretrain view of the Tahoe pseudobulk corpus.

Contract: `.scratch/stage2-05-prnet-pretrain/spec.md` §Implementation Decisions
(agreed 2026-11-01).
This is a *separate* entry point from the corpus builder in
`stage2/02_dataset_corpus/`: it never rebuilds the corpus. It projects the frozen
62,710-gene Tahoe corpus (`stage2/02_dataset_corpus/output/corpus/v2/tahoe.h5ad`)
onto the frozen 18,164-gene CM axis (`.../corpus/v2/tables/genes.csv`) and writes
one h5ad under this step's own `output/pretrain_view/`, whose schema is identical
to the frozen CM views (`stage2/_frozen/cm-data-v1/datasets-1/`), plus the side
files a downstream trainer needs.

Layout differs from CM on purpose (deviation D0): CM views are dense int32; this
view is CSR.  Schema (obs, var, layers, uns) is identical; only the container
differs.  The plan's original size argument ("dense = 4.87 TB") was a units
error: a dense 67,018 x 18,164 int32 X is 4.87 GB (9.74 GB with the float32
layer), and at 83% density (1,010,892,356 of 1,217,314,952 cells) the CSR
artifact is in fact the larger container (16.32 GB).  CSR is kept because it
mirrors the corpus's per-plate CSR parts and lets the builder assemble with
slice writes — a recorded deviation, not a size win.

Usage (run from the repository root; any cwd works, the root is found by
walking up to `.git`):

    python stage2/05_prnet_pretrain/execution/build_pretrain_view.py --check
    python stage2/05_prnet_pretrain/execution/build_pretrain_view.py            # build + QA
    python stage2/05_prnet_pretrain/execution/build_pretrain_view.py --force    # rebuild _view_parts/
    python stage2/05_prnet_pretrain/execution/build_pretrain_view.py --qa       # QA only, non-zero on failure
    python stage2/05_prnet_pretrain/execution/build_pretrain_view.py --verify   # replay the content hash

Output tree (`output/pretrain_view/`):

    tahoe_cm_schema.h5ad   local-only (gitignored by the repo-wide `output/` rule)
    manifest.csv           tracked
    compounds.csv          tracked
    qa.json                tracked
    qa_report.md           tracked
    _view_parts/           local-only intermediate cache (gitignored)

Two passes keep peak RSS small (~2 GiB) on a 31 GiB host: pass 1 gathers the
retained columns of each `_parts/plate<N>.npz` into `_view_parts/`; pass 2 fills
the preallocated h5ad CSR arrays with slice writes, one part in memory at a time.

Contract deviations, recorded in the step README and ADR 0031:
    D0  layout: CSR instead of CM's dense int32 (container, not size — see above).
    D1  control keying: `control_group = tahoe::<plate>::<cell_line_id>` (700
        groups) instead of the corpus's `tahoe::<plate>` (14), because a plate
        pools ~50 cell lines and a plate-level control mean is the mixture
        ADR 0029's "Rejected (a)" rejects.
    D2  gene axis: 18,164 slots; `FAM153B` / ENSG00000182230 has no Tahoe row,
        kept as an all-zero row so the axis stays byte-identical to CM.
    D3  provenance gaps: `time_h` NaN (the atlas has no time column),
        `cell_type`/`donor` NaN, `atc`/`moa`/`mesh_classes`/
        `annotation_source` "", `pharm_*` False — recorded here instead of
        invented constants.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


# ---------------------------------------------------------------------------
# layout: resolve the repo root by walking up to .git (never count levels)
# ---------------------------------------------------------------------------
def find_repo_root(start: Path) -> Path:
    """Walk up from *start* to the repository root (the directory holding `.git`)."""
    for candidate in (start, *start.parents):
        if (candidate / ".git").exists():
            return candidate
    raise SystemExit(f"cannot find the repository root: no .git at or above {start}")


ROOT = find_repo_root(Path(__file__).resolve().parent)
assert (ROOT / "stage2").is_dir(), f"resolved repo root has no stage2/: {ROOT}"

CORPUS_H5AD = ROOT / "stage2" / "02_dataset_corpus" / "output" / "corpus" / "v2" / "tahoe.h5ad"
CORPUS_PARTS = ROOT / "stage2" / "02_dataset_corpus" / "execution" / "_parts"
OUT_DIR = ROOT / "stage2" / "05_prnet_pretrain" / "output" / "pretrain_view"
VIEW_PARTS = OUT_DIR / "_view_parts"
VIEW_H5AD = OUT_DIR / "tahoe_cm_schema.h5ad"
GENES_CSV = (ROOT / "stage2" / "02_dataset_corpus" / "output" / "corpus" / "v2"
             / "tables" / "genes.csv")
CM_COMPOUNDS = (ROOT / "stage2" / "02_dataset_corpus" / "output" / "corpus" / "v2"
                / "tables" / "compounds.csv")
TAHOE_META = ROOT / "pert_data" / "Tahoe-100M" / "source" / "metadata"

VIEW_NAME = "tahoe_cm_schema"
BUILT_FROM = "stage2/02_dataset_corpus/output/corpus/v2/tahoe.h5ad"

#: `content_sha256()` of the view as recorded by the move (ticket #14) and by the
#: 06.3 input contract.  `--verify` replays the projection and compares to this.
RECORDED_CONTENT_SHA256 = "476812c09d06c29e56069a3a1e670834384abb524a9468dc596229b5cbcbf407"
PSEUDOBULK_UNIT = "(plate, well, cell_line)"
CORPUS_VERSION = "v1"

N_VIEW_GENES = 18_164
N_VIEW_OBS = 67_018
N_CONTROLS = 1_400
N_DRUGS = 379
N_CONTROL_GROUPS = 700
CONTROL_DRUG = "DMSO_TF"

# CM's 24 obs columns, in the CM contract's order (the `well` column duplicates the index)
CM_OBS_COLUMNS = [
    "well", "dataset_id", "compound_name", "drug_key", "inchikey", "smiles",
    "dose_um", "dose_source", "time_h", "is_control", "vehicle", "plate",
    "control_group", "cell_type", "donor", "replicate", "atc", "moa",
    "mesh_classes", "pharm_atc_moa", "pharm_mesh", "annotation_source",
    "lib_size", "n_genes_detected",
]
TAHOE_OBS_COLUMNS = [
    "cell_line_id", "cell_name", "n_cells", "moa_broad", "moa_fine",
    "pubchem_cid", "has_structure", "has_matched_control",
]
SHIM_COLUMNS = ["SMILES", "dose", "paired_control_index"]
FORBIDDEN_OBS_COLUMNS = {"pert_idx", "cell_idx", "dose_idx", "pert_dose_idx"}

# recorded overlap routes (ticket 02 / the ADR 0029 amendment), asserted so a drift is
# loud, not silent.  The spec/ticket state only 31 desalted / 24 raw / 36 CID / 34 name;
# the union 42 is a *derived* number (the three identity routes agree on 42 distinct CM
# keys) and is labeled as derived wherever it is written.
EXPECTED_DESALTED_KEY_OVERLAP = 31
EXPECTED_RAW_KEY_OVERLAP = 24
EXPECTED_CID_KEY_OVERLAP = 36
EXPECTED_NORMALISED_NAME_KEY_OVERLAP = 34
EXPECTED_UNION_KEY_OVERLAP = 42  # derived: distinct CM keys in parent | CID | name
CM_N_DRUG_KEYS = 525


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def log(msg: str) -> None:
    print(msg, flush=True)


def short_commit() -> str:
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short=7", "HEAD"],
                              check=True, capture_output=True, text=True).stdout.strip()
    except Exception:
        return "unknown"


# ---------------------------------------------------------------------------
# frozen inputs: CM gene axis, Tahoe gene map, drug metadata, corpus obs
# ---------------------------------------------------------------------------
def load_cm_var() -> pd.DataFrame:
    """The 18,164-row CM `var`, copied verbatim from `corpus/v2/tables/genes.csv`."""
    genes = pd.read_csv(GENES_CSV, dtype=str)
    assert len(genes) == N_VIEW_GENES, len(genes)
    assert list(genes.columns) == ["symbol", "hgnc_id", "entrez_id", "ensembl_id", "match_type"]
    assert not genes["symbol"].duplicated().any()
    var = genes[["entrez_id", "ensembl_id", "match_type"]].copy()
    var.index = pd.Index(genes["symbol"].to_numpy(), name="symbol")
    # CM writes match_type as a categorical; keep the same category order
    var["match_type"] = pd.Categorical(var["match_type"], categories=["alias", "direct", "entrez"])
    return var


def build_colmap(var: pd.DataFrame) -> np.ndarray:
    """Tahoe column -> CM var position, joined on ensembl_id only (-1 = dropped).

    Joining on symbol would work for 17,935 genes; ensembl_id keeps 18,163.
    `FAM153B` (ENSG00000182230) has no Tahoe row and stays -1.
    """
    genes = pq.read_table(TAHOE_META / "gene_metadata.parquet").to_pandas()
    assert len(genes) == 62_710
    assert not genes["ensembl_id"].duplicated().any()
    position = {e: i for i, e in enumerate(var["ensembl_id"])}
    colmap = np.full(len(genes), -1, dtype=np.int32)
    for j, ensembl in enumerate(genes["ensembl_id"]):
        target = position.get(ensembl)
        if target is not None:
            colmap[j] = target
    return colmap


def load_drugs() -> pd.DataFrame:
    """One row per Tahoe drug (379), with raw/desalted InChIKey and annotations.

    Two well names carry a trailing space in `sample_metadata` ('Erdafitinib ',
    'Selinexor '); the corpus kept them verbatim, so the join to `drug_metadata`
    is done on the stripped name.  Two drugs have no structure
    (Sacubitril/Valsartan, Verteporfin).
    """
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    meta = pd.read_parquet(TAHOE_META / "drug_metadata.parquet")
    chooser = rdMolStandardize.LargestFragmentChooser()
    rows = []
    for row in meta.itertuples(index=False):
        smiles = row.canonical_smiles
        smiles = smiles if isinstance(smiles, str) and smiles else ""
        raw_inchikey = parent_inchikey = ""
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise SystemExit(f"drug_metadata: cannot parse SMILES for {row.drug!r}")
            raw_inchikey = Chem.MolToInchiKey(mol)
            parent_inchikey = Chem.MolToInchiKey(chooser.choose(mol))
            for key in (raw_inchikey, parent_inchikey):
                if len(key) != 27 or key.count("-") != 2:
                    raise SystemExit(f"drug_metadata: bad InChIKey {key!r} for {row.drug!r}")
        rows.append({
            "drug": row.drug,
            "smiles": smiles,
            "raw_inchikey": raw_inchikey,
            "parent_inchikey": parent_inchikey,
            "pubchem_cid": "" if pd.isna(row.pubchem_cid) else str(int(row.pubchem_cid)),
        })
    drugs = pd.DataFrame(rows)
    # `itertuples` renames hyphenated columns; re-read the two string columns by name
    drugs["moa_broad"] = meta["moa-broad"].fillna("").astype(str).to_numpy()
    drugs["moa_fine"] = meta["moa-fine"].fillna("").astype(str).to_numpy()
    drugs["has_structure"] = drugs["smiles"].ne("")
    drugs["join_name"] = drugs["drug"].str.strip()
    assert drugs["join_name"].is_unique, "trailing-space strip made drug names ambiguous"
    return drugs


def _decode(value):
    """HDF5 stores strings as bytes; normalise to str."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def read_obs_h5py(path: Path) -> pd.DataFrame:
    """Read only the `obs` group of an h5ad, without materialising X or layers.

    `anndata.read_h5ad(backed="r")` eagerly loads every layer as a real sparse
    matrix; on the 25.8 GB corpus that costs 12.8 GB of RSS for a 67k-row table.
    Obs/var are small, so they are decoded directly here instead.
    """
    import h5py

    with h5py.File(path, "r") as handle:
        group = handle["obs"]
        index_name = _decode(group.attrs["_index"])
        index = pd.Index([_decode(v) for v in group[index_name][:]], name=index_name)  # dataset named after the index
        columns = [_decode(c) for c in group.attrs["column-order"]]
        data = {}
        for column in columns:
            node = group[column]
            if isinstance(node, h5py.Group):  # categorical
                categories = np.array([_decode(v) for v in node["categories"][:]])
                data[column] = pd.Categorical.from_codes(node["codes"][:], categories=categories)
            else:
                values = node[:]
                if values.dtype == object:
                    values = np.array([_decode(v) for v in values], dtype=object)
                data[column] = values
    return pd.DataFrame(data, index=index)[columns]  # keep the stored column order


def read_var_h5py(path: Path) -> pd.DataFrame:
    """Read only the `var` group of an h5ad, with CM's dtypes."""
    import h5py

    with h5py.File(path, "r") as handle:
        group = handle["var"]
        index_name = _decode(group.attrs["_index"])
        index = pd.Index([_decode(v) for v in group[index_name][:]], name=index_name)  # dataset named after the index
        data = {}
        for column in [_decode(c) for c in group.attrs["column-order"]]:
            node = group[column]
            if isinstance(node, h5py.Group):
                categories = np.array([_decode(v) for v in node["categories"][:]])
                data[column] = pd.Categorical.from_codes(node["codes"][:], categories=categories)
            else:
                values = node[:]
                if values.dtype == object:
                    values = np.array([_decode(v) for v in values], dtype=object)
                data[column] = values
    return pd.DataFrame(data, index=index)[
        ["entrez_id", "ensembl_id", "match_type"]]


def read_control_mean_h5py(path: Path) -> pd.DataFrame:
    """Decode `uns['control_mean']` (anndata dataframe encoding, column-wise)."""
    import h5py

    with h5py.File(path, "r") as handle:
        group = handle["uns/control_mean"]
        index_name = _decode(group.attrs["_index"])
        columns = [_decode(c) for c in group.attrs["column-order"]]
        index = [_decode(v) for v in group[index_name][:]]
        values = np.empty((len(index), len(columns)), dtype=np.float64)
        for j, column in enumerate(columns):
            values[:, j] = group[column][:]
    return pd.DataFrame(values, index=pd.Index(index, name=index_name), columns=columns)


def anndata_roundtrip_check(path: Path) -> dict:
    """Prove `anndata.read_h5ad` recovers `uns['control_mean']` as a DataFrame.

    Runs in a subprocess so the QA process does not pay for the eager layer load.
    """
    import subprocess

    code = (
        "import json, anndata as ad\n"
        f"a = ad.read_h5ad({str(path)!r}, backed='r')\n"
        "m = a.uns['control_mean']\n"
        "print(json.dumps({'anndata': ad.__version__, 'type': type(m).__name__,\n"
        "                  'shape': [int(m.shape[0]), int(m.shape[1])],\n"
        "                  'index_name': m.index.name,\n"
        "                  'n_columns': int(len(m.columns)),\n"
        "                  'first_columns': [str(c) for c in list(m.columns[:2])],\n"
        "                  'last_column': str(list(m.columns)[-1])}))\n"
    )
    try:
        result = subprocess.run([sys.executable, "-c", code], check=True,
                                capture_output=True, text=True, timeout=900)
        payload = json.loads(result.stdout.strip().splitlines()[-1])
        payload["ok"] = True
        return payload
    except Exception as error:  # noqa: BLE001 - recorded as data, not raised
        return {"ok": False, "error": repr(error)[:500]}


def load_corpus_obs() -> pd.DataFrame:
    """Corpus obs only (67,018 x 9); never reads X or the 12.8 GB layer."""
    obs = read_obs_h5py(CORPUS_H5AD)
    assert len(obs) == N_VIEW_OBS, len(obs)
    return obs


def view_control_group(obs: pd.DataFrame) -> pd.Series:
    """The view's D1 control key: `tahoe::<plate>::<cell_line_id>` (700 groups).

    The corpus' own `tahoe::<plate>` (14 groups) is not reused: a plate pools
    ~50 cell lines, so it would mix cell lines inside one control mean.
    """
    return ("tahoe::" + obs["plate"].astype(str) + "::"
            + obs["cell_line_id"].astype(str))


def input_parts() -> list[Path]:
    parts = sorted(CORPUS_PARTS.glob("plate*.npz"),
                   key=lambda p: int(p.stem[len("plate"):]))
    if len(parts) != 14:
        raise SystemExit(f"expected 14 corpus parts, found {len(parts)}")
    return parts


def verify_row_order(obs: pd.DataFrame, parts: list[Path]) -> None:
    """The view row order is plate order then part row order; assert the corpus agrees."""
    plates, wells, lines, ncells = [], [], [], []
    for part in parts:
        data = np.load(part, allow_pickle=False)
        plates.append(np.full(len(data["row_ncell"]), str(data["plate"])))
        wells.append(data["wells"][data["row_well_ix"]])
        lines.append(data["line_ids"][data["row_line_ix"]])
        ncells.append(data["row_ncell"])
    plates = np.concatenate(plates)
    wells = np.concatenate(wells)
    lines = np.concatenate(lines)
    ncells = np.concatenate(ncells)
    checks = {
        "plate": (plates, obs["plate"].astype(str).to_numpy()),
        "well": (wells, obs["well"].astype(str).to_numpy()),
        "cell_line_id": (lines, obs["cell_line_id"].astype(str).to_numpy()),
        "n_cells": (ncells, obs["n_cells"].to_numpy()),
    }
    for name, (got, want) in checks.items():
        if len(got) != len(want) or not np.array_equal(got, want):
            raise SystemExit(f"corpus row order disagrees with obs on {name!r}; "
                             "refusing to build a misaligned view")


# ---------------------------------------------------------------------------
# pass 1: gather retained columns part by part (the only full-corpus read)
# ---------------------------------------------------------------------------
def _view_part_path(part: Path) -> Path:
    return VIEW_PARTS / f"{part.stem}.npz"


def _ctrl_part_path(part: Path) -> Path:
    return VIEW_PARTS / f"{part.stem}.ctrl.npz"


def build_view_parts(colmap: np.ndarray, obs: pd.DataFrame, force: bool) -> list[Path]:
    import scipy.sparse as sp

    VIEW_PARTS.mkdir(parents=True, exist_ok=True)
    parts = input_parts()
    view_groups = view_control_group(obs)
    groups = sorted(view_groups[obs["is_control"]].unique())
    if len(groups) != N_CONTROL_GROUPS:
        raise SystemExit(f"expected {N_CONTROL_GROUPS} control groups, got {len(groups)}")
    group_codes = pd.Categorical(view_groups, categories=groups).codes.astype(np.int32)
    is_control = obs["is_control"].to_numpy()
    row0 = 0
    t0 = time.time()
    for part in parts:
        dst, ctrl_dst = _view_part_path(part), _ctrl_part_path(part)
        n_part = len(np.load(part, allow_pickle=False)["indptr"]) - 1
        if dst.exists() and ctrl_dst.exists() and not force:
            log(f"  pass1 {part.stem}: cached ({n_part:,} rows)")
            row0 += n_part
            continue
        src = np.load(part, allow_pickle=False)
        indptr = src["indptr"].astype(np.int64)
        values = src["data"]
        indices = src["indices"]
        if not np.array_equal(np.diff(indptr) > 0, np.ones(n_part, dtype=bool)):
            raise SystemExit(f"{part.stem}: empty source row; the corpus contract assumes none")
        target = colmap[indices]
        keep = (target >= 0) & (values != 0)
        kept_target = target[keep].astype(np.int32)
        kept_values = values[keep].astype(np.int32)
        # per-row retained nnz; np.add.reduceat is safe because no row is empty
        counts = np.add.reduceat(keep, indptr[:-1], dtype=np.int32).astype(np.int64)
        newptr = np.zeros(n_part + 1, dtype=np.int64)
        np.cumsum(counts, out=newptr[1:])
        matrix = sp.csr_matrix((kept_values, kept_target, newptr),
                               shape=(n_part, N_VIEW_GENES))
        matrix.sort_indices()
        del target, keep, kept_target, kept_values, counts, newptr, src, indices, values, indptr
        row_sum = np.asarray(matrix.sum(axis=1)).ravel().astype(np.int64)
        row_detected = matrix.getnnz(axis=1).astype(np.int32)
        if (row_sum <= 0).any():
            raise SystemExit(f"{part.stem}: {int((row_sum <= 0).sum())} rows lost every gene")
        src = np.load(part, allow_pickle=False)  # small row maps, reloaded after freeing buffers
        np.savez(
            dst,
            data=matrix.data.astype(np.int32, copy=False),
            indices=matrix.indices.astype(np.int32, copy=False),
            indptr=matrix.indptr.astype(np.int64),
            nnz=np.int64(matrix.nnz), n_rows=np.int64(n_part),
            row_sum=row_sum, row_ndetected=row_detected,
            row_well_ix=src["row_well_ix"], row_line_ix=src["row_line_ix"],
            row_ncell=src["row_ncell"], wells=src["wells"],
            line_ids=src["line_ids"], plate=np.array(str(src["plate"])),
        )
        # control rows only, for uns["control_mean"] (2 vehicle wells per group)
        ctrl_local = np.flatnonzero(is_control[row0:row0 + n_part])
        sub = matrix[ctrl_local]
        np.savez(
            ctrl_dst,
            groups=group_codes[row0 + ctrl_local],
            data=sub.data.astype(np.int32, copy=False),
            indices=sub.indices.astype(np.int32, copy=False),
            indptr=sub.indptr.astype(np.int64),
            n_rows=np.int64(len(ctrl_local)),
        )
        log(f"  pass1 {part.stem}: {n_part:,} rows, {matrix.nnz:,} retained nnz, "
            f"{len(ctrl_local):,} control rows ({time.time() - t0:.0f}s)")
        del matrix, sub, ctrl_local, row_sum, row_detected, src
        row0 += n_part
    return parts


def load_part_metadata(parts: list[Path]) -> dict:
    """Small arrays only: shapes, row sums and row maps (never `data`)."""
    meta = {"n_rows": [], "nnz": [], "row_sum": [], "row_ndetected": [],
            "row_well_ix": [], "row_line_ix": [], "row_ncell": [],
            "plate": [], "wells_rows": [], "line_ids_rows": []}
    for part in parts:
        data = np.load(_view_part_path(part), allow_pickle=False)
        meta["n_rows"].append(int(data["indptr"].shape[0]) - 1)
        meta["nnz"].append(int(data["nnz"]))
        meta["row_sum"].append(data["row_sum"])
        meta["row_ndetected"].append(data["row_ndetected"])
        meta["row_well_ix"].append(data["row_well_ix"])
        meta["row_line_ix"].append(data["row_line_ix"])
        meta["row_ncell"].append(data["row_ncell"])
        meta["plate"].append(str(data["plate"]))
        meta["wells_rows"].append(data["wells"][data["row_well_ix"]])
        meta["line_ids_rows"].append(data["line_ids"][data["row_line_ix"]])
    return meta


def compute_control_mean(parts: list[Path], obs: pd.DataFrame, var: pd.DataFrame):
    groups = sorted(view_control_group(obs)[obs["is_control"]].unique())
    sums = np.zeros((len(groups), N_VIEW_GENES), dtype=np.float64)
    counts = np.zeros(len(groups), dtype=np.int64)
    for part in parts:
        data = np.load(_ctrl_part_path(part), allow_pickle=False)
        indptr = data["indptr"]
        values, indices, codes = data["data"], data["indices"], data["groups"]
        for r in range(len(codes)):
            start, stop = indptr[r], indptr[r + 1]
            sums[codes[r], indices[start:stop]] += values[start:stop]
            counts[codes[r]] += 1
    missing = np.flatnonzero(counts == 0)
    if len(missing):
        raise SystemExit(f"control groups without a vehicle control: {[groups[i] for i in missing]}")
    mean = pd.DataFrame(sums / counts[:, None],
                        index=pd.Index(groups, name="control_group"),
                        columns=var.index)
    return mean


# ---------------------------------------------------------------------------
# obs: CM's 24 columns + Tahoe extras + author shims
# ---------------------------------------------------------------------------
def build_obs(obs: pd.DataFrame, var: pd.DataFrame, drugs: pd.DataFrame,
              part_meta: dict) -> pd.DataFrame:
    n = len(obs)
    plate = obs["plate"].astype(str).to_numpy()
    raw_well = obs["well"].astype(str).to_numpy()
    line = obs["cell_line_id"].astype(str).to_numpy()
    is_control = obs["is_control"].to_numpy()
    drug = obs["drug"].astype(str).to_numpy()
    # corpus rows keep the upstream names verbatim ('Erdafitinib ', ...); join stripped
    joined = drugs.set_index("join_name")
    stripped = pd.Series(drug).str.strip().to_numpy()
    missing = sorted(set(stripped) - set(joined.index) - {CONTROL_DRUG})
    if missing:
        raise SystemExit(f"drug names with no drug_metadata row: {missing}")
    parent = pd.Series(stripped).map(joined["parent_inchikey"]).fillna("").to_numpy()
    raw = pd.Series(stripped).map(joined["raw_inchikey"]).fillna("").to_numpy()
    smiles = pd.Series(stripped).map(joined["smiles"]).fillna("").to_numpy()
    cid = pd.Series(stripped).map(joined["pubchem_cid"]).fillna("").to_numpy()
    moa_broad = pd.Series(stripped).map(joined["moa_broad"]).fillna("").to_numpy()
    moa_fine = pd.Series(stripped).map(joined["moa_fine"]).fillna("").to_numpy()
    # controls carry no structure: empty keys and has_structure False
    parent = np.where(is_control, "", parent)
    raw = np.where(is_control, "", raw)
    smiles = np.where(is_control, "", smiles)
    cid = np.where(is_control, "", cid)
    moa_broad = np.where(is_control, "", moa_broad)
    moa_fine = np.where(is_control, "", moa_fine)
    has_structure = np.where(is_control, False, smiles != "")

    index = pd.Index("tahoe::" + pd.Series(plate) + ":" + pd.Series(raw_well)
                     + ":" + pd.Series(line), name="well")
    view_well = index.to_numpy()  # obs["well"] mirrors the index (CM convention)
    if index.duplicated().any():
        raise SystemExit(f"{int(index.duplicated().sum())} duplicate view row ids")
    control_group = view_control_group(obs)

    # within-group replicate index: group = (control_group, drug_key, dose_um), ordered by well
    dose = obs["dose_um"].astype(float).to_numpy()
    group_key = (control_group.to_numpy() + "|" + parent + "|" + pd.Series(dose).astype(str).to_numpy())
    order = np.lexsort((view_well, group_key))
    replicate = np.zeros(n, dtype=np.int32)
    last_key, running = None, 0
    for position in order:
        key = group_key[position]
        if key != last_key:
            running, last_key = 0, key
        replicate[position] = running
        running += 1

    # first control row of each group (ordered by well); control rows point to themselves
    paired = np.empty(n, dtype=object)
    control_positions = np.flatnonzero(is_control)
    control_order = np.lexsort((view_well[control_positions],
                                control_group.to_numpy()[control_positions]))
    first_control = {}
    for position in control_positions[control_order]:
        first_control.setdefault(control_group.iloc[position], index[position])
    paired[:] = [first_control[group] for group in control_group.to_numpy()]
    # Documented shim semantics (and the pretrain fixture encoding): a control row points
    # to its own label; only treated rows point at their group's first control row.
    paired[is_control] = np.asarray(index)[is_control]

    lib_size = np.concatenate(part_meta["row_sum"]).astype(np.int64)
    n_detected = np.concatenate(part_meta["row_ndetected"]).astype(np.int32)
    if len(lib_size) != n:
        raise SystemExit(f"row arrays cover {len(lib_size)} rows, obs has {n}")

    obs_out = pd.DataFrame({
        "well": view_well,
        "dataset_id": pd.Categorical(["tahoe"] * n),
        "compound_name": pd.Categorical(drug),
        "drug_key": pd.Categorical(parent),
        "inchikey": pd.Categorical(raw),
        "smiles": pd.Categorical(smiles),
        "dose_um": dose,
        "dose_source": pd.Categorical(["nominal_atlas"] * n),
        "time_h": np.full(n, np.nan),
        "is_control": is_control,
        "vehicle": pd.Categorical(np.where(is_control, CONTROL_DRUG, "")),
        "plate": pd.Categorical(plate),
        "control_group": pd.Categorical(control_group),
        "cell_type": pd.Categorical([None] * n),
        "donor": pd.Categorical([None] * n),
        "replicate": replicate,
        "atc": pd.Categorical([""] * n),
        "moa": pd.Categorical([""] * n),
        "mesh_classes": pd.Categorical([""] * n),
        "pharm_atc_moa": np.zeros(n, dtype=bool),
        "pharm_mesh": np.zeros(n, dtype=bool),
        "annotation_source": pd.Categorical([""] * n),
        "lib_size": lib_size,
        "n_genes_detected": n_detected,
    }, index=index)
    # Tahoe-native columns (appended after the CM block)
    obs_out["cell_line_id"] = pd.Categorical(line)
    obs_out["cell_name"] = pd.Categorical(obs["cell_name"].astype(str).to_numpy())
    obs_out["n_cells"] = obs["n_cells"].astype(np.int32).to_numpy()
    obs_out["moa_broad"] = pd.Categorical(moa_broad)
    obs_out["moa_fine"] = pd.Categorical(moa_fine)
    obs_out["pubchem_cid"] = cid.astype(object)
    obs_out["has_structure"] = has_structure.astype(bool)
    has_control = control_group.isin(first_control).to_numpy()
    obs_out["has_matched_control"] = has_control
    # author-compatibility shims (marked as shims in qa_report.md, not contract columns)
    obs_out["SMILES"] = smiles.astype(object)
    obs_out["dose"] = dose
    obs_out["paired_control_index"] = paired.astype(object)
    expected = CM_OBS_COLUMNS + TAHOE_OBS_COLUMNS + SHIM_COLUMNS
    assert list(obs_out.columns) == expected, list(obs_out.columns)
    return obs_out


# ---------------------------------------------------------------------------
# pass 2: preallocate the h5ad CSR arrays and fill them part by part
# ---------------------------------------------------------------------------
def _csr_group(parent, name: str, n_rows: int, nnz: int, dtype, n_genes: int = N_VIEW_GENES):
    """An anndata-compatible CSR group (same encoding as build_pseudobulk.py)."""
    group = parent.create_group(name)
    group.attrs["encoding-type"] = "csr_matrix"
    group.attrs["encoding-version"] = "0.1.0"
    group.attrs["shape"] = np.array([n_rows, n_genes], dtype=np.int64)
    data = group.create_dataset("data", shape=(nnz,), dtype=dtype, chunks=(2_000_000,))
    indices = group.create_dataset("indices", shape=(nnz,), dtype=np.int32, chunks=(2_000_000,))
    indptr = group.create_dataset("indptr", shape=(n_rows + 1,), dtype=np.int64)
    return group, data, indices, indptr


def assemble(parts: list[Path], var: pd.DataFrame, obs, control_mean: pd.DataFrame,
             part_meta: dict) -> None:
    import anndata as ad
    import h5py
    import scipy.sparse as sp

    n_rows = sum(part_meta["n_rows"])
    total_nnz = sum(part_meta["nnz"])
    if n_rows != N_VIEW_OBS:
        raise SystemExit(f"assembled rows {n_rows} != {N_VIEW_OBS}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    skeleton = ad.AnnData(
        X=sp.csr_matrix((n_rows, N_VIEW_GENES), dtype=np.int32),
        obs=obs, var=var)
    skeleton.layers["log1p_cp10k"] = sp.csr_matrix((n_rows, N_VIEW_GENES), dtype=np.float32)
    skeleton.uns["control_mean"] = control_mean
    skeleton.write_h5ad(VIEW_H5AD)
    del skeleton
    log(f"  pass2: skeleton written ({time.time() - t0:.0f}s)")

    with h5py.File(VIEW_H5AD, "a") as handle:
        for stale in ("X", "layers"):
            if stale in handle:
                del handle[stale]
        _, x_data, x_indices, x_indptr = _csr_group(handle, "X", n_rows, total_nnz, np.int32)
        layers = handle.create_group("layers")
        _, l_data, l_indices, l_indptr = _csr_group(layers, "log1p_cp10k", n_rows, total_nnz,
                                                    np.float32)
        indptr = np.zeros(n_rows + 1, dtype=np.int64)
        row0 = nz0 = 0
        for part in parts:
            data = np.load(_view_part_path(part), allow_pickle=False)
            part_ptr = data["indptr"]
            part_nnz = int(data["nnz"])
            n_part = int(data["n_rows"])
            values = data["data"]
            indices = data["indices"]
            x_data[nz0:nz0 + part_nnz] = values
            x_indices[nz0:nz0 + part_nnz] = indices
            # same formula as the frozen CM view writer
            # (stage2/_frozen/cm-data-v1/execution/cmdata/views.py): counts / lib_size * 1e4
            lib_per_nnz = np.repeat(data["row_sum"], np.diff(part_ptr))
            log1p = np.log1p(values / lib_per_nnz * 1e4).astype(np.float32)
            l_data[nz0:nz0 + part_nnz] = log1p
            l_indices[nz0:nz0 + part_nnz] = indices
            indptr[row0:row0 + n_part] = nz0 + part_ptr[:-1]
            row0 += n_part
            nz0 += part_nnz
            log(f"  pass2 {part.stem}: {n_part:,} rows, {part_nnz:,} nnz "
                f"({time.time() - t0:.0f}s)")
            del data, values, indices, log1p, lib_per_nnz, part_ptr
        if row0 != n_rows or nz0 != total_nnz:
            raise SystemExit(f"assembled {row0}/{n_rows} rows, {nz0}/{total_nnz} nnz")
        indptr[-1] = nz0
        x_indptr[:] = indptr
        l_indptr[:] = indptr
        handle.attrs["view"] = VIEW_NAME
        handle.attrs["built_from"] = BUILT_FROM
        handle.attrs["pseudobulk_unit"] = PSEUDOBULK_UNIT
    log(f"  pass2: wrote {VIEW_H5AD} ({VIEW_H5AD.stat().st_size / 1e9:.1f} GB, "
        f"{time.time() - t0:.0f}s)")


# ---------------------------------------------------------------------------
# side files: content hash, manifest, compounds, contract expectations
# ---------------------------------------------------------------------------
def content_sha256(path: Path, chunk: int = 32_000_000) -> str:
    """Hash the CSR `data`/`indices`/`indptr` arrays separately, chunked.

    Never densify: `np.ascontiguousarray(X)` would allocate ~9.7 GB
    (67,018 x 18,164 float64) and defeat the streaming assembly.
    """
    import h5py

    digest = hashlib.sha256()
    with h5py.File(path, "r") as handle:
        for name in ("data", "indices", "indptr"):
            dataset = handle[f"X/{name}"]
            for start in range(0, dataset.shape[0], chunk):
                block = dataset[start:start + chunk]
                digest.update(np.ascontiguousarray(block).tobytes())
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# optional `--verify`: replay the projection and recompute the content hash
# ---------------------------------------------------------------------------
def _project_part(part: Path, colmap: np.ndarray):
    """Project one corpus part onto the CM axis in memory (verify only).

    Mirrors the retained-column rule of `build_view_parts` exactly — same
    `colmap`, same zero drop, same CSR row pointers, same `sort_indices` — but
    writes no `_view_parts/` cache.  Returns the CSR matrix and its per-row
    retained nnz.
    """
    import scipy.sparse as sp

    src = np.load(part, allow_pickle=False)
    indptr = src["indptr"].astype(np.int64)
    values, indices = src["data"], src["indices"]
    target = colmap[indices]
    keep = (target >= 0) & (values != 0)
    kept_target = target[keep].astype(np.int32)
    kept_values = values[keep].astype(np.int32)
    counts = np.add.reduceat(keep, indptr[:-1], dtype=np.int32).astype(np.int64)
    newptr = np.zeros(len(indptr), dtype=np.int64)
    np.cumsum(counts, out=newptr[1:])
    matrix = sp.csr_matrix((kept_values, kept_target, newptr),
                           shape=(len(indptr) - 1, N_VIEW_GENES))
    matrix.sort_indices()
    return matrix, counts


def verify_content_hash(chunk: int = 32_000_000) -> int:
    """Recompute the view's `content_sha256` from the corpus + tables and compare it.

    The recorded hash covers only the CSR `X/data`, `X/indices` and `X/indptr`
    arrays (see `content_sha256`), so this replays the pass-1 projection over the
    per-plate corpus parts in the same order and hashes the arrays that a rebuild
    would write — without creating a second 16 GB h5ad.  Two passes keep peak RSS
    at one part: pass A hashes `data` while accumulating `indptr`; pass B
    re-derives and hashes `indices`.  Missing inputs exit with a message, not a
    traceback; a mismatch exits non-zero.
    """
    required = [CORPUS_H5AD, GENES_CSV, TAHOE_META / "gene_metadata.parquet"]
    missing = [path for path in required if not path.exists()]
    if not CORPUS_PARTS.is_dir():
        missing.append(CORPUS_PARTS)
    if missing:
        raise SystemExit("--verify needs inputs that are missing:\n  "
                         + "\n  ".join(str(path) for path in missing))

    parts = input_parts()
    obs = load_corpus_obs()
    verify_row_order(obs, parts)
    colmap = build_colmap(load_cm_var())

    t0 = time.time()
    digest = hashlib.sha256()
    indptr = np.zeros(N_VIEW_OBS + 1, dtype=np.int64)
    row0 = 0
    for part in parts:
        matrix, counts = _project_part(part, colmap)
        n_part = matrix.shape[0]
        if row0 + n_part > N_VIEW_OBS:
            raise SystemExit(f"{part.stem}: {row0 + n_part} rows exceeds {N_VIEW_OBS}")
        base = indptr[row0]
        indptr[row0 + 1:row0 + n_part + 1] = base + np.cumsum(counts)
        data = matrix.data
        for start in range(0, data.shape[0], chunk):
            digest.update(np.ascontiguousarray(data[start:start + chunk]).tobytes())
        row0 += n_part
        log(f"  verify pass A {part.stem}: {n_part:,} rows, {data.shape[0]:,} nnz")
        del matrix, data, counts
    if row0 != N_VIEW_OBS:
        raise SystemExit(f"replayed {row0} rows, the view records {N_VIEW_OBS}")

    for part in parts:
        matrix, _ = _project_part(part, colmap)
        indices = matrix.indices
        for start in range(0, indices.shape[0], chunk):
            digest.update(np.ascontiguousarray(indices[start:start + chunk]).tobytes())
        log(f"  verify pass B {part.stem}: {indices.shape[0]:,} indices")
        del matrix, indices

    for start in range(0, indptr.shape[0], chunk):
        digest.update(np.ascontiguousarray(indptr[start:start + chunk]).tobytes())

    got = digest.hexdigest()
    match = got == RECORDED_CONTENT_SHA256
    log(f"verify: recomputed content_sha256 = {got}")
    log(f"verify: recorded   content_sha256 = {RECORDED_CONTENT_SHA256}")
    log(f"verify: {'MATCH' if match else 'MISMATCH'} ({time.time() - t0:.0f}s)")
    return 0 if match else 1


def build_compounds(drugs: pd.DataFrame) -> pd.DataFrame:
    """379 rows, one per Tahoe drug; `overlap_cm` = desalted parent key in CM's 525."""
    cm = pd.read_csv(CM_COMPOUNDS)
    cm_keys = set(cm["drug_key"].dropna().astype(str))
    if len(cm_keys) != CM_N_DRUG_KEYS:
        raise SystemExit(f"expected {CM_N_DRUG_KEYS} CM drug keys, got {len(cm_keys)}")
    compounds = drugs[["drug", "smiles", "pubchem_cid", "parent_inchikey",
                       "raw_inchikey", "moa_broad", "moa_fine", "has_structure"]].copy()
    compounds["overlap_cm"] = compounds["parent_inchikey"].isin(cm_keys)
    compounds = compounds.sort_values("drug", kind="stable").reset_index(drop=True)
    assert len(compounds) == N_DRUGS
    return compounds


def overlap_numbers(drugs: pd.DataFrame) -> dict:
    """Reproduce the recorded overlap routes on distinct CM drug keys.

    Ticket 02 / the ADR 0029 amendment record 31 desalted, 24 raw (non-desalted),
    36 by PubChem CID and 34 by alphanumerically normalised name, against the CM
    corpus's 525 keys.  The union of the three identity routes (42) is a *derived*
    number, recorded but not attributed to the spec.  A per-row `overlap_cm`
    boolean has 33 True rows (two pairs of Tahoe drug names share one desalted
    parent), so "31" is the distinct-key count, not the row count; both are
    recorded.
    """
    cm = pd.read_csv(CM_COMPOUNDS)
    cm_keys = set(cm["drug_key"].dropna().astype(str))
    cm_cids = {int(x) for x in cm["pubchem_cid"].dropna()}

    def normalise(name: str) -> str:
        # Alphanumeric-only rule: lowercase, every non-`[a-z0-9]` character removed.
        # This is the route that reproduces the recorded 34 name matches (it folds
        # salt/parenthesis/hyphen spellings such as `5-Fluorouracil` vs
        # `5 Fluorouracil`); lower/trim alone finds only 31.
        return re.sub(r"[^a-z0-9]", "", str(name).lower())

    parent_keys = set(drugs["parent_inchikey"]) - {""}
    raw_keys = set(drugs["raw_inchikey"]) - {""}
    drug_cids = {int(x) for x in drugs["pubchem_cid"] if x}
    drug_names = {normalise(x) for x in drugs["drug"]}

    keys_parent = {k for k in cm_keys if k in parent_keys}
    keys_raw = {k for k in cm_keys if k in raw_keys}
    keys_cid = {k for k in cm["drug_key"][cm["pubchem_cid"].fillna(-1).astype("int64")
                                          .isin(drug_cids)]}
    keys_name = {k for k in cm["drug_key"][cm["compound_name"].map(normalise)
                                          .isin(drug_names)]}
    return {
        "cm_drug_keys": len(cm_keys),
        "desalted_parent_keys_matched": len(keys_parent),
        "raw_nondesalted_keys_matched": len(keys_raw),
        "cid_keys_matched": len(keys_cid),
        "normalised_name_keys_matched": len(keys_name),
        "union_keys_matched": len(keys_parent | keys_cid | keys_name),
        "rows_with_overlap_cm": int((drugs["parent_inchikey"].isin(cm_keys)).sum()),
        "drug_rows_sharing_one_parent_key": int(drugs["parent_inchikey"].isin(cm_keys).sum())
        - len(parent_keys & cm_keys),
    }


def check_contract_numbers(overlap: dict) -> None:
    checks = [
        ("desalted parent-key overlap", overlap["desalted_parent_keys_matched"],
         EXPECTED_DESALTED_KEY_OVERLAP),
        ("raw non-desalted overlap", overlap["raw_nondesalted_keys_matched"],
         EXPECTED_RAW_KEY_OVERLAP),
        ("PubChem CID overlap", overlap["cid_keys_matched"],
         EXPECTED_CID_KEY_OVERLAP),
        ("normalised-name overlap", overlap["normalised_name_keys_matched"],
         EXPECTED_NORMALISED_NAME_KEY_OVERLAP),
        # The union is a derived number, not a spec figure; asserted so a route drift
        # cannot silently move it (the report labels it derived).
        ("derived union (parent/CID/name) overlap", overlap["union_keys_matched"],
         EXPECTED_UNION_KEY_OVERLAP),
    ]
    bad = [(n, got, want) for n, got, want in checks if got != want]
    if bad:
        detail = "; ".join(f"{n}: measured {got}, recorded {want}" for n, got, want in bad)
        raise SystemExit("overlap numbers drifted from the recorded routes: " + detail)


def build_manifest(obs, var: pd.DataFrame, control_mean: pd.DataFrame,
                   compounds: pd.DataFrame) -> pd.DataFrame:
    n = len(obs)
    is_control = obs["is_control"].to_numpy()
    manifest = pd.DataFrame([{
        "view": VIEW_NAME,
        "file": VIEW_H5AD.name,
        "n_obs": n,
        "n_perturbed": int((~is_control).sum()),
        "n_controls": int(is_control.sum()),
        "n_drugs": int(obs.loc[~is_control, "compound_name"].nunique()),
        "n_genes": len(var),
        "n_control_groups": len(control_mean),
        "corpus_version": CORPUS_VERSION,
        "source_commit": short_commit(),
        "content_sha256": content_sha256(VIEW_H5AD),
        "bytes": VIEW_H5AD.stat().st_size,
    }])
    return manifest


# ---------------------------------------------------------------------------
# QA: the nine hard gates, run against the written file
# ---------------------------------------------------------------------------
class XReader:
    """Row-block reader over an h5ad CSR matrix, without densifying the whole thing."""

    def __init__(self, path: Path, name: str = "X"):
        import h5py

        self.handle = h5py.File(path, "r")
        self.group = self.handle[name]
        self.indptr = self.group["indptr"][:]
        self.shape = tuple(int(x) for x in self.group.attrs["shape"])
        self.n_rows = self.shape[0]

    def close(self):
        self.handle.close()

    def rows(self, rows: list[int]) -> np.ndarray:
        """Dense float64 copy of the requested rows (used on small samples only)."""
        dense = np.zeros((len(rows), self.shape[1]), dtype=np.float64)
        data, indices = self.group["data"], self.group["indices"]
        for out, row in enumerate(rows):
            start, stop = int(self.indptr[row]), int(self.indptr[row + 1])
            dense[out, indices[start:stop]] = data[start:stop]
        return dense

    def row_sums_and_detected(self):
        """Exact per-row sum and positive count, read in row blocks."""
        data = self.group["data"]
        n = self.n_rows
        sums = np.zeros(n, dtype=np.int64)
        detected = np.zeros(n, dtype=np.int32)
        block = 1024
        for row0 in range(0, n, block):
            row1 = min(row0 + block, n)
            start, stop = int(self.indptr[row0]), int(self.indptr[row1])
            values = data[start:stop]
            starts = self.indptr[row0:row1 + 1] - start
            cumulative = np.concatenate(([0], np.cumsum(values, dtype=np.int64)))
            sums[row0:row1] = cumulative[starts[1:]] - cumulative[starts[:-1]]
            positive = np.concatenate(([0], np.cumsum(values > 0, dtype=np.int64)))
            detected[row0:row1] = (positive[starts[1:]] - positive[starts[:-1]]).astype(np.int32)
        return sums, detected

    def min_max(self):
        data = self.group["data"]
        low, high = None, None
        step = 32_000_000
        for start in range(0, data.shape[0], step):
            block = data[start:start + step]
            low = int(block.min()) if low is None else min(low, int(block.min()))
            high = int(block.max()) if high is None else max(high, int(block.max()))
        return low, high


def old_anndata_roundtrip(path: Path) -> dict:
    """Read `uns['control_mean']` with the anndata 0.10.8 benchmark interpreter.

    The PRnet benchmark env pins anndata 0.10.8, which cannot *write* an
    18k-column DataFrame but reads ours fine.  Recorded (not gated): if the
    interpreter is absent the check is skipped rather than failed.
    """
    import subprocess

    candidates = [Path("/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python"),
                  Path("/home/chang/miniconda3/envs/arc-stage2-bringup-prnet/bin/python")]
    interpreter = next((c for c in candidates if c.exists()), None)
    if interpreter is None:
        return {"skipped": True, "reason": "no anndata 0.10.8 benchmark interpreter on this host"}
    code = (
        "import json, anndata as ad\n"
        f"a = ad.read_h5ad({str(path)!r}, backed='r')\n"
        "m = a.uns['control_mean']\n"
        "print(json.dumps({'anndata': ad.__version__, 'type': type(m).__name__,\n"
        "                  'shape': [int(m.shape[0]), int(m.shape[1])],\n"
        "                  'index_name': m.index.name}))\n"
    )
    try:
        result = subprocess.run([str(interpreter), "-c", code], check=True,
                                capture_output=True, text=True, timeout=600)
        payload = json.loads(result.stdout.strip().splitlines()[-1])
        payload.update({"interpreter": str(interpreter), "ok": True})
        return payload
    except Exception as error:  # noqa: BLE001 - recorded as data, not raised
        return {"ok": False, "interpreter": str(interpreter), "error": repr(error)[:500]}


def run_qa(colmap: np.ndarray) -> dict:
    import h5py

    gates = []
    qa_interop: dict = {}

    def gate(number, name, measured, passed):
        gates.append({"gate": number, "name": name, "pass": bool(passed),
                      "measured": measured})
        return bool(passed)

    if not VIEW_H5AD.exists():
        raise SystemExit(f"no view to QA: {VIEW_H5AD}")
    genes = pd.read_csv(GENES_CSV, dtype=str)
    obs = read_obs_h5py(VIEW_H5AD)
    var = read_var_h5py(VIEW_H5AD)
    means = read_control_mean_h5py(VIEW_H5AD)
    anndata_roundtrip = anndata_roundtrip_check(VIEW_H5AD)
    n_obs, n_genes = obs.shape[0], var.shape[0]
    is_control = obs["is_control"].to_numpy()
    n_controls = int(is_control.sum())
    n_perturbed = int((~is_control).sum())
    groups = obs["control_group"].astype(str)

    # 1. shape / control arithmetic / row identity
    measured = {"n_obs": n_obs, "n_controls": n_controls, "n_perturbed": n_perturbed,
                "n_drugs": int(obs.loc[~is_control, "compound_name"].nunique()),
                "n_genes": n_genes, "n_control_groups": int(groups.nunique()),
                "index_name": obs.index.name,
                "well_equals_index": bool((obs["well"].astype(str) == obs.index.astype(str)).all()),
                "index_unique": bool(obs.index.is_unique)}
    gate(1, "shape, controls, perturbed, control groups, row identity", measured,
         measured["n_obs"] == N_VIEW_OBS and measured["n_controls"] == N_CONTROLS
         and measured["n_perturbed"] == N_VIEW_OBS - N_CONTROLS
         and measured["n_genes"] == N_VIEW_GENES
         and measured["n_drugs"] == N_DRUGS
         and measured["n_control_groups"] == N_CONTROL_GROUPS
         and measured["index_name"] == "well" and measured["well_equals_index"]
         and measured["index_unique"])

    # 2. var is byte-identical to CM's frozen gene axis
    var_columns_ok = list(var.columns) == ["entrez_id", "ensembl_id", "match_type"]
    var_ok = (var_columns_ok
              and list(var.index) == list(genes["symbol"])
              and list(var["entrez_id"].astype(str)) == list(genes["entrez_id"])
              and list(var["ensembl_id"].astype(str)) == list(genes["ensembl_id"])
              and list(var["match_type"].astype(str)) == list(genes["match_type"])
              and list(var["match_type"].cat.categories) == ["alias", "direct", "entrez"]
              and var.index.name == "symbol")
    gate(2, "var identical to corpus/v2/tables/genes.csv (values, order, dtypes)", {
        "rows": n_genes, "index_name": var.index.name,
        "columns": list(var.columns), "columns_ok": var_columns_ok,
        "dtypes": {c: str(t) for c, t in var.dtypes.items()},
        "identical": var_ok}, var_ok)

    # 3 + 5. one streaming pass over X: lib_size, detected, integer/non-negative
    x_read = XReader(VIEW_H5AD, "X")
    try:
        sums, detected = x_read.row_sums_and_detected()
        low, high = x_read.min_max()
        x_dtype = str(x_read.group["data"].dtype)
        measured = {"lib_size_max_abs_diff": int(np.abs(sums - obs["lib_size"].to_numpy()).max()),
                    "n_genes_detected_max_abs_diff":
                        int(np.abs(detected - obs["n_genes_detected"].to_numpy()).max()),
                    "x_dtype": x_dtype,
                    "x_min": low, "x_max": high, "nnz": int(x_read.indptr[-1])}
        gate(3, "lib_size == X.sum(1); n_genes_detected == (X>0).sum(1)", measured,
             measured["lib_size_max_abs_diff"] == 0
             and measured["n_genes_detected_max_abs_diff"] == 0)
        gate(5, "X non-negative int32 counts (subset spot-checked separately)", measured,
             low >= 0 and x_dtype == "int32" and high < 2**31)

        # 4. log1p_cp10k recomputes from X on a random row sample
        rng = np.random.default_rng(20261021)
        sample = np.sort(rng.choice(n_obs, size=200, replace=False))
        dense_x = x_read.rows(sample.tolist())
        dense_l = XReader(VIEW_H5AD, "layers/log1p_cp10k")
        try:
            got = dense_l.rows(sample.tolist())
        finally:
            dense_l.close()
        expected = np.log1p(dense_x / sums[sample][:, None] * 1e4)
        max_diff = float(np.abs(got - expected).max()) if got.size else 0.0
        gate(4, "log1p_cp10k recomputes from X (atol=1e-5, 200 rows)", {
            "rows_sampled": 200, "max_abs_diff": max_diff, "atol": 1e-5},
            max_diff <= 1e-5)

        # 5b. subset spot-check: 3 rows against the corpus `_parts` counts
        spot_rows = np.sort(rng.choice(n_obs, size=3, replace=False))
        parts = input_parts()
        part_rows = np.concatenate([[0], np.cumsum([
            len(np.load(p, allow_pickle=False)["indptr"]) - 1 for p in parts])])
        mismatches = []
        compared = 0
        for row in spot_rows:
            part_index = int(np.searchsorted(part_rows, row, side="right") - 1)
            local = int(row - part_rows[part_index])
            part = np.load(parts[part_index], allow_pickle=False)
            start, stop = int(part["indptr"][local]), int(part["indptr"][local + 1])
            corpus_indices = part["indices"][start:stop]
            target = colmap[corpus_indices]
            keep = (target >= 0) & (part["data"][start:stop] != 0)
            corpus_row = dict(zip(target[keep].tolist(), part["data"][start:stop][keep].tolist()))
            start, stop = int(x_read.indptr[row]), int(x_read.indptr[row + 1])
            view_row = dict(zip(x_read.group["indices"][start:stop].tolist(),
                                x_read.group["data"][start:stop].tolist()))
            if not corpus_row:
                continue
            compared += 1
            if corpus_row != view_row:
                common = set(corpus_row) | set(view_row)
                diff = [(k, corpus_row.get(k), view_row.get(k)) for k in common
                        if corpus_row.get(k) != view_row.get(k)]
                mismatches.append({"row": int(row), "n_diff": len(diff), "example": diff[:3]})
        gate(5.1, "X is a pure subset of the corpus counts (3 random rows)", {
            "rows": [int(r) for r in spot_rows], "compared": compared,
            "mismatches": mismatches},
            compared > 0 and not mismatches)

        # 6. control_group coverage, control_mean exactness, matched controls
        groups_sorted = sorted(groups.unique())
        mean_axis_ok = (isinstance(means, pd.DataFrame) and means.shape == (N_CONTROL_GROUPS, N_VIEW_GENES)
                        and list(means.columns) == list(var.index)
                        and means.index.name == "control_group"
                        and list(means.index) == groups_sorted
                        and means.to_numpy().dtype == np.float64)
        treated_groups = set(groups[~is_control])
        control_groups = set(groups[is_control])
        matched = treated_groups <= control_groups
        sample_group = groups_sorted[123]
        rows = np.flatnonzero((groups == sample_group).to_numpy() & is_control)
        recomputed = x_read.rows(rows.tolist()).mean(axis=0)
        mean_diff = float(np.abs(recomputed - means.loc[sample_group].to_numpy()).max())
        with h5py.File(VIEW_H5AD, "r") as handle:
            encoding = handle["uns/control_mean"].attrs
            mean_encoding = {
                "encoding_type": str(encoding.get("encoding-type", "")),
                "encoding_version": str(encoding.get("encoding-version", "")),
                "_index": str(encoding.get("_index", "")),
            }
        interop = old_anndata_roundtrip(VIEW_H5AD)
        # The declared round-trip contract (§Testing Decisions): a plain consumer read
        # must return the declared 700 x 18,164 frame with the declared index name.
        roundtrip_ok = bool(
            anndata_roundtrip.get("ok")
            and anndata_roundtrip.get("type") == "DataFrame"
            and anndata_roundtrip.get("shape") == [N_CONTROL_GROUPS, N_VIEW_GENES]
            and anndata_roundtrip.get("index_name") == "control_group")
        gate(6, "control_mean covers 700 groups, matches recomputed means, round-trips", {
            "control_mean_shape": list(means.shape) if means is not None else None,
            "index_name": getattr(getattr(means, "index", None), "name", None),
            "dtype": str(means.to_numpy().dtype) if means is not None else None,
            "axis_ok": bool(mean_axis_ok),
            "encoding": mean_encoding,
            "roundtrip_through_anndata": anndata_roundtrip,
            "roundtrip_ok": roundtrip_ok,
            "treated_groups_without_control": len(treated_groups - control_groups),
            "sampled_group": sample_group, "sampled_group_rows": int(len(rows)),
            "max_abs_diff": mean_diff},
            bool(mean_axis_ok and matched and mean_diff <= 1e-9 and roundtrip_ok))
        # recorded interop (not one of the nine gates): can a consumer on anndata 0.10.8 read it?
        qa_interop = interop

        # 7. no forbidden encodings; the D3 empty-value policy and the structure columns
        forbidden = sorted(FORBIDDEN_OBS_COLUMNS & set(obs.columns))
        expected_empty = is_control | (~obs["has_structure"].to_numpy())
        structure_columns = ["drug_key", "inchikey", "smiles", "SMILES", "pubchem_cid"]
        structure_mismatches = {
            column: int(((obs[column].astype(str) == "").to_numpy() != expected_empty).sum())
            for column in structure_columns}
        # D3: every column the atlas cannot supply stays empty — NaN, "" or False as the
        # deviation declares — and the three derived/constant columns carry only their
        # declared value.  A drift to an invented constant fails here, not silently.
        nan_columns = ["time_h", "cell_type", "donor"]
        empty_columns = ["atc", "moa", "mesh_classes", "annotation_source"]
        false_columns = ["pharm_atc_moa", "pharm_mesh"]
        constant_columns = {"dataset_id": "tahoe", "dose_source": "nominal_atlas"}
        d3 = {
            "nan_columns": nan_columns,
            "nan_all_missing": {c: int(obs[c].isna().sum()) for c in nan_columns},
            "empty_columns": empty_columns,
            "non_empty_counts": {c: int((obs[c].astype(str) != "").sum())
                                 for c in empty_columns},
            "false_columns": false_columns,
            "true_counts": {c: int(obs[c].to_numpy(dtype=bool).sum()) for c in false_columns},
            "constant_columns": constant_columns,
            "constant_mismatches": {c: int((obs[c].astype(str) != v).sum())
                                    for c, v in constant_columns.items()},
            # vehicle: DMSO_TF exactly on control rows, "" on every treated row
            "vehicle_mismatches": int(((obs["vehicle"].astype(str) == CONTROL_DRUG)
                                       .to_numpy() != is_control).sum()),
        }
        d3_ok = (all(v == n_obs for v in d3["nan_all_missing"].values())
                 and all(v == 0 for v in d3["non_empty_counts"].values())
                 and all(v == 0 for v in d3["true_counts"].values())
                 and all(v == 0 for v in d3["constant_mismatches"].values())
                 and d3["vehicle_mismatches"] == 0)
        gate(7, "no forbidden obs columns; structure + D3 emptiness policy exact", {
            "forbidden": forbidden,
            "structure_columns": structure_columns,
            "mismatches_by_column": structure_mismatches,
            "rows_expected_empty": int(expected_empty.sum()),
            **d3},
            not forbidden and all(value == 0 for value in structure_mismatches.values())
            and d3_ok)

        # 8. paired_control_index resolves inside the obs index; shim dose marks controls
        paired = obs["paired_control_index"].astype(str).to_numpy()
        known = set(obs.index.astype(str))
        unresolved = int(sum(1 for value in paired if value not in known))
        own = bool((paired[is_control] == obs.index.astype(str).to_numpy()[is_control]).all())
        dose = obs["dose"].to_numpy().astype(float)
        dose_zero_controls = bool((dose[is_control] == 0.0).all())
        dose_nonzero_treated = bool((dose[~is_control] > 0.0).all())
        gate(8, "paired_control_index resolves; shim dose is 0.0 exactly on controls", {
            "unresolved": unresolved,
            "control_rows_point_to_themselves": own,
            "dose_zero_on_controls": dose_zero_controls,
            "dose_nonzero_on_treated": dose_nonzero_treated},
            unresolved == 0 and own and dose_zero_controls and dose_nonzero_treated)
    finally:
        x_read.close()

    # 2b. root attributes (recorded, not one of the nine gates)
    with h5py.File(VIEW_H5AD, "r") as handle:
        root_attrs = {key: str(handle.attrs[key]) for key in handle.attrs}
    all_pass = all(g["pass"] for g in gates)
    return {
        "checked_at": now(),
        "view": VIEW_NAME,
        "file": VIEW_H5AD.name,
        "bytes": VIEW_H5AD.stat().st_size,
        "shape": [int(n_obs), int(n_genes)],
        "root_attrs": root_attrs,
        "gate_count": len(gates),
        "gates_passed": int(sum(1 for g in gates if g["pass"])),
        "gates": sorted(gates, key=lambda g: float(g["gate"])),
        "all_pass": all_pass,
        "interop_anndata_0_10_8": qa_interop,
    }


# ---------------------------------------------------------------------------
# qa.json / qa_report.md
# ---------------------------------------------------------------------------
DEVIATIONS = [
    ("D0", "matrix layout",
     "CM views are dense int32; this view is CSR int32 (+ CSR float32 layer)",
     "mirrors the corpus's per-plate CSR parts and the slice-write assembly; the "
     "plan's 'dense = 4.87 TB' was a units error — dense is 4.87 GB (9.74 GB with "
     "the layer) and the matrix is 83% dense, so this 16.32 GB CSR artifact is "
     "larger than a dense one"),
    ("D1", "control keying",
     "`control_group` = `tahoe::<plate>::<cell_line_id>` (700 groups), "
     "`uns['control_mean']` has 700 rows",
     "a plate pools ~50 cell lines; a plate-level mean is the Simpson's-paradox "
     "mixture ADR 0029's \"Rejected (a)\" rejects, so Δ stays inside one cell line"),
    ("D2", "gene axis",
     "18,164 CM slots; `FAM153B` / ENSG00000182230 is an all-zero row",
     "the gene has no Tahoe measurement; dropping it would break axis identity"),
    ("D3", "provenance gaps",
     "`time_h` NaN; `cell_type`/`donor` NaN; `atc`/`moa`/`mesh_classes`/"
     "`annotation_source` \"\"; `pharm_atc_moa`/`pharm_mesh` False; and the derived or "
     "constant columns `dataset_id` \"tahoe\", `dose_source` \"nominal_atlas\", `vehicle` "
     "`DMSO_TF` on controls / \"\" on treated, `inchikey` derived from the SMILES rather "
     "than looked up",
     "the atlas carries no time column and no CM annotation vocabulary, and it has no "
     "InChIKey; recorded rather than invented"),
]


def write_qa_files(qa: dict, var: pd.DataFrame, drugs: pd.DataFrame,
                   overlap: dict, control_mean: pd.DataFrame) -> None:
    genes = pd.read_csv(GENES_CSV, dtype=str)
    tahoe_genes = pq.read_table(TAHOE_META / "gene_metadata.parquet").to_pandas()
    measured_by_ensembl = len(set(genes["ensembl_id"]) & set(tahoe_genes["ensembl_id"]))
    measured_by_symbol = len(set(genes["symbol"]) & set(tahoe_genes["gene_symbol"]))
    unmeasured = sorted(set(genes["ensembl_id"]) - set(tahoe_genes["ensembl_id"]))
    null_drugs = drugs.loc[~drugs["has_structure"], "drug"].tolist()
    compounds = build_compounds(drugs)
    qa_out = dict(qa)
    qa_out.update({
        "contract_expectations": {
            "desalted_parent_keys_matched": {
                "expected": EXPECTED_DESALTED_KEY_OVERLAP,
                "measured": overlap["desalted_parent_keys_matched"],
                "match": overlap["desalted_parent_keys_matched"] == EXPECTED_DESALTED_KEY_OVERLAP},
            "raw_nondesalted_keys_matched": {
                "expected": EXPECTED_RAW_KEY_OVERLAP,
                "measured": overlap["raw_nondesalted_keys_matched"],
                "match": overlap["raw_nondesalted_keys_matched"] == EXPECTED_RAW_KEY_OVERLAP},
            "cid_keys_matched": {
                "expected": EXPECTED_CID_KEY_OVERLAP,
                "measured": overlap["cid_keys_matched"],
                "match": overlap["cid_keys_matched"] == EXPECTED_CID_KEY_OVERLAP},
            "normalised_name_keys_matched": {
                "expected": EXPECTED_NORMALISED_NAME_KEY_OVERLAP,
                "measured": overlap["normalised_name_keys_matched"],
                "match": overlap["normalised_name_keys_matched"]
                == EXPECTED_NORMALISED_NAME_KEY_OVERLAP},
        },
        "derived_overlap": {
            "union_parent_cid_name_keys_matched": {
                "measured": overlap["union_keys_matched"],
                "note": "union of the three identity routes (parent | CID | name); a "
                        "derived number, not stated in the spec or ticket (which give "
                        "only 31 desalted / 24 raw / 36 CID / 34 name)"},
        },
        "gene_mapping": {
            "cm_axis_rows": len(var),
            "measured_by_ensembl": measured_by_ensembl,
            "measured_by_symbol": measured_by_symbol,
            "unmeasured_ensembl": unmeasured,
            "join": "Tahoe -> CM on ensembl_id only; the symbol route matches "
                    f"{measured_by_symbol} (loses {len(var) - measured_by_symbol} of the "
                    f"{len(var)} CM genes; {measured_by_ensembl - measured_by_symbol} "
                    f"relative to the {measured_by_ensembl} that map by Ensembl)",
        },
        "null_structure_drugs": null_drugs,
        "overlap": dict(overlap, cm_compounds=len(compounds)),
        "shims": {
            "SMILES": "= smiles (author Dataset.py reads 'SMILES')",
            "dose": "= dose_um (0.0 marks controls; pert_model/PRnet/source/data/Dataset.py:39)",
            "paired_control_index": "obs label of the first control row of the row's "
                                    "control_group (own label on control rows)",
        },
        "control_mean": {"shape": list(control_mean.shape),
                         "dtype": str(control_mean.to_numpy().dtype),
                         "groups": len(control_mean),
                         "encoding": next(g["measured"]["encoding"]
                                          for g in qa["gates"] if g["gate"] == 6),
                         "roundtrip_through_anndata":
                             next(g["measured"]["roundtrip_through_anndata"]
                                  for g in qa["gates"] if g["gate"] == 6)},
    })
    (OUT_DIR / "qa.json").write_text(json.dumps(qa_out, indent=2) + "\n")

    rows = []
    # report in gate-number order (qa.json keeps the evaluation order)
    for gate in sorted(qa["gates"], key=lambda g: float(g["gate"])):
        rows.append(f"| {gate['gate']} | {gate['name']} | "
                    f"{'PASS' if gate['pass'] else 'FAIL'} |")
    checks = qa_out["contract_expectations"]
    dev_rows = [f"| {code} | {title} | {actual} | {why} |" for code, title, actual, why in DEVIATIONS]
    interop = qa.get("interop_anndata_0_10_8") or {}
    if interop.get("ok"):
        interop_text = (f"verified — anndata {interop.get('anndata')} read it back as a "
                        f"`{interop.get('type')}` {interop.get('shape')} with index "
                        f"`{interop.get('index_name')}` (interpreter "
                        f"`{interop.get('interpreter')}`)")
    elif interop.get("skipped"):
        interop_text = f"not verified on this host: {interop.get('reason')}"
    else:
        interop_text = (f"**FAILED** on `{interop.get('interpreter')}`: "
                        f"{interop.get('error')}")
    mean_encoding = qa_out["control_mean"]["encoding"]
    import anndata as _anndata

    report = f"""# Pretrain view — QA report

QA run at {qa['checked_at']} by `stage2/05_prnet_pretrain/execution/build_pretrain_view.py`.
Contract: `.scratch/stage2-05-prnet-pretrain/spec.md` §Implementation Decisions
("the view's contract"). Source corpus:
`{BUILT_FROM}` (never rebuilt; the view only projects the frozen rows onto the frozen
18,164-gene CM axis).

## Gates (nine)

| # | gate | result |
|---|---|---|
{chr(10).join(rows)}

Nine gates: **{qa['gates_passed']}/{qa['gate_count']} PASS** (`all_pass={qa['all_pass']}`).
`--qa` exits non-zero when any gate fails.
The machine-readable measurements are in [`qa.json`](qa.json).

## Deviations from the CM contract (D0–D3)

| id | topic | this view | why |
|---|---|---|---|
{chr(10).join(dev_rows)}

## Gene mapping (verbatim CM axis, joined on `ensembl_id`)

- CM axis: **{len(var)}** slots, order and values copied from `corpus/v2/tables/genes.csv`.
- Measured by Ensembl: **{measured_by_ensembl}/{len(var)}**; by symbol:
  **{measured_by_symbol}/{len(var)}** — the symbol route loses
  {len(var) - measured_by_symbol} of the {len(var)} CM genes
  ({measured_by_ensembl - measured_by_symbol} relative to the {measured_by_ensembl} that map
  by Ensembl), so the join is on Ensembl only.
- Unmeasured: {', '.join(f'`{e}`' for e in unmeasured)} (`FAM153B`) — kept as an
  all-zero row so the axis stays byte-identical to every CM view (D2).

## Null-structure drugs

{len(null_drugs)} of {len(drugs)} Tahoe drugs have no `canonical_smiles` (the other
{len(drugs) - len(null_drugs)} resolve to a structure):
{', '.join(f'`{d}`' for d in null_drugs)}. Their `has_structure` is False and
`drug_key`/`inchikey`/`smiles`/`SMILES`/`pubchem_cid` are `""`; the matching says
nothing about biology.

## Drug overlap with the CM corpus

`compounds.csv` has {len(compounds)} rows and `overlap_cm` = *the desalted parent
InChIKey is one of the {overlap['cm_drug_keys']} CM `drug_key`s*.
Raw and desalted keys are both stored (`raw_inchikey`, `parent_inchikey`).

| route | distinct CM drug keys (of {overlap['cm_drug_keys']}) | definition |
|---|---|---|
| desalted parent InChIKey | **{overlap['desalted_parent_keys_matched']}** | `LargestFragmentChooser` + `MolToInchiKey` on `canonical_smiles`; the corpus's own identity rule (ticket/spec: 31) |
| raw, non-desalted InChIKey | **{overlap['raw_nondesalted_keys_matched']}** | `MolToInchiKey` without desalting (ADR 0029's "24/525") |
| PubChem CID | **{overlap['cid_keys_matched']}** | joined on `pubchem_cid` (ticket/spec: 36) |
| alphanumerically normalised name | **{overlap['normalised_name_keys_matched']}** | lowercase, every non-`[a-z0-9]` character removed; joined on `compound_name` (ticket/spec: 34) |
| **union of the three identity routes** | **{overlap['union_keys_matched']}** | derived: parent ∪ CID ∪ name — **not** stated in the spec or ticket |

Counting is on **distinct CM drug keys**, not on Tahoe rows: `{overlap['rows_with_overlap_cm']}`
Tahoe drug rows match by desalted parent key but
{overlap['drug_rows_sharing_one_parent_key']} of them share a key with another row
(e.g. `Trametinib` / `Trametinib (DMSO_TF solvate)`), which is why `overlap_cm`
has {overlap['rows_with_overlap_cm']} True rows and {overlap['desalted_parent_keys_matched']}
distinct keys. The recorded routes are reproduced exactly (expected/measured:
desalted {checks['desalted_parent_keys_matched']['expected']}/{checks['desalted_parent_keys_matched']['measured']},
raw {checks['raw_nondesalted_keys_matched']['expected']}/{checks['raw_nondesalted_keys_matched']['measured']},
CID {checks['cid_keys_matched']['expected']}/{checks['cid_keys_matched']['measured']},
name {checks['normalised_name_keys_matched']['expected']}/{checks['normalised_name_keys_matched']['measured']}).
The union 42 is a *derived* number recorded from those routes; it is not a spec figure.

## Author-compatibility shims (not contract columns)

- `SMILES` = `smiles` (the author `Dataset.py` reads `SMILES`).
- `dose` = `dose_um`, so `0.0` marks a control, matching
  `pert_model/PRnet/source/data/Dataset.py:39`.
- `paired_control_index` = obs label of the first control row of the row's
  `control_group`; on control rows it is the row's own label.

## `uns['control_mean']` encoding (consumer compatibility)

`uns['control_mean']` is an anndata **{mean_encoding['encoding_type']}/{mean_encoding['encoding_version']}**
DataFrame (column-wise storage, `_index = {mean_encoding['_index']}`):
{control_mean.shape[0]} control groups × {control_mean.shape[1]} genes, float64, index named
`control_group`, columns = the 18,164 var symbols in var order. Gate 6 re-reads it through a
plain `anndata.read_h5ad` consumer read, requires the declared 700 x 18,164 shape and
`control_group` index name, and compares it with a group mean recomputed from `X`; the writer
here is anndata {_anndata.__version__}. The PRnet benchmark environment pins anndata 0.10.8,
which cannot *write* an 18k-column DataFrame but reads this one: {interop_text} (recorded,
not one of the nine gates).

## Shape

`{qa['shape'][0]:,} x {qa['shape'][1]:,}` CSR int32 + `layers['log1p_cp10k']` CSR float32,
`uns['control_mean']` float64 {control_mean.shape[0]} x {control_mean.shape[1]} (D1),
{qa['bytes'] / 1e9:.2f} GB on disk. Content SHA-256 in `manifest.csv` covers the
CSR `data`/`indices`/`indptr` arrays separately (chunked; never densified).
"""
    (OUT_DIR / "qa_report.md").write_text(report)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true",
                        help="print resolved paths and exit")
    parser.add_argument("--force", action="store_true",
                        help="rebuild the _view_parts/ cache (a missing part is rebuilt anyway)")
    parser.add_argument("--qa", action="store_true",
                        help="QA the existing view only; exits non-zero if a gate fails")
    parser.add_argument("--verify", action="store_true",
                        help="replay the projection from the corpus + tables, compare "
                             "the recomputed content hash to the recorded value and exit "
                             "(default off; reads the parts, writes nothing)")
    args = parser.parse_args(argv)

    log(f"repo root : {ROOT}")
    log(f"corpus    : {CORPUS_H5AD}")
    log(f"corpus parts: {CORPUS_PARTS}")
    log(f"gene axis : {GENES_CSV}")
    log(f"output    : {OUT_DIR}")
    if args.verify:
        return verify_content_hash()
    if args.check:
        return 0

    drugs = load_drugs()
    overlap = overlap_numbers(drugs)
    var = load_cm_var()
    colmap = build_colmap(var)

    if args.qa:
        qa = run_qa(colmap)
        control_mean = _read_control_mean_for_report()
        write_qa_files(qa, var, drugs, overlap, control_mean)
        _print_gate_table(qa)
        return 0 if qa["all_pass"] else 1

    if not CORPUS_H5AD.exists():
        raise SystemExit(f"corpus missing: {CORPUS_H5AD} (this step never rebuilds it)")
    obs_in = load_corpus_obs()
    parts = input_parts()
    verify_row_order(obs_in, parts)

    t0 = time.time()
    log("pass 1: gather the retained columns")
    build_view_parts(colmap, obs_in, force=args.force)
    part_meta = load_part_metadata(parts)
    control_mean = compute_control_mean(parts, obs_in, var)
    log(f"  control_mean: {control_mean.shape} over {N_CONTROL_GROUPS} (plate, cell_line) groups")
    obs = build_obs(obs_in, var, drugs, part_meta)
    log(f"  obs: {obs.shape[0]:,} rows x {obs.shape[1]} columns")
    log("pass 2: preallocate and fill the h5ad")
    assemble(parts, var, obs, control_mean, part_meta)

    compounds = build_compounds(drugs)
    compounds.to_csv(OUT_DIR / "compounds.csv", index=False)
    check_contract_numbers(overlap)
    manifest = build_manifest(obs, var, control_mean, compounds)
    manifest.to_csv(OUT_DIR / "manifest.csv", index=False)
    log(f"  manifest: {manifest.to_dict('records')[0]}")

    qa = run_qa(colmap)
    write_qa_files(qa, var, drugs, overlap, control_mean)
    _print_gate_table(qa)
    log(f"done in {time.time() - t0:.0f}s; all gates {'PASS' if qa['all_pass'] else 'FAIL'}")
    return 0 if qa["all_pass"] else 1


def _read_control_mean_for_report() -> pd.DataFrame:
    """`--qa` mode: recover control_mean from the written file (column-wise decode)."""
    return read_control_mean_h5py(VIEW_H5AD)


def _print_gate_table(qa: dict) -> None:
    log("QA gates:")
    for gate in sorted(qa["gates"], key=lambda g: float(g["gate"])):
        log(f"  [{ 'PASS' if gate['pass'] else 'FAIL' }] {gate['gate']}: {gate['name']}")
    log(f"  all_pass={qa['all_pass']}")


if __name__ == "__main__":
    sys.exit(main())
