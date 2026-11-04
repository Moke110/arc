#!/usr/bin/env python
"""Build the frozen inputs for the Stage 2 `temp_env_all_5k` benchmark.

Derived from the executed `temp_env_all_18k` builder.  It reads the curated CM
view once, keeps the 5,000 genes selected by `build_gene_set.py`, writes the
Δ/groups/controls/true-set/split artifacts, runs the twelve gates (the ten split
gates plus matrix/true-set integrity) and the five metric sanity cases, and
records every hash in ``manifest.json``.

The rows, Δ, controls, groups and splits must be identical to the 18k archive;
``--source-data`` verifies that column-subset equality before anything is
written.

Usage
-----
    python stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/build_data.py
    python stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/build_data.py --check
    python stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/build_data.py --out /tmp/bench-data-rerun

The default output is ``stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/output/input/``.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import sys
import warnings
import zipfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from importlib.metadata import version as _pkg_version

sys.path.insert(0, str(Path(__file__).resolve().parent))

from repo_root import BENCH, ROOT  # noqa: E402

# The v1 evidence is frozen (ticket #13): the producer is retired, the files
# are read-only in stage2/_frozen/cm-data-v1/ and never rebuilt here.
FROZEN = ROOT / "stage2/_frozen/cm-data-v1"
SOURCE = FROZEN / "datasets-1/temp_env_all-smiles.h5ad"
SOURCE_MANIFEST = FROZEN / "datasets-1/manifest.csv"
GENES = FROZEN / "corpus/v1/genes.csv"
DEFAULT_OUT = BENCH / "output" / "input"
DEFAULT_GENE_SET = DEFAULT_OUT / "gene_set_5k.csv"
DEFAULT_SOURCE_DATA = BENCH.parent / "03.1_cold_drug_18k" / "output" / "input"

SEED = "20261021"
K_LIST = (100, 200, 500)
N_FOLDS = 5
SOURCE_GENES = 18_164
N_GENES = 5_000
N_TREATED = 4_206
N_GROUPS = 1_380
N_CONTROL_GROUPS = 16

# Metric-sanity goldens for the first group on the frozen 5k axis.  They are
# data-dependent regression values (the reverse AUPRC floor and the
# top/bottom-only correlations), recomputed once for the 5k gene set with
# seed 20261021; the 18k builder used the same style of golden.
SANITY_GOLDEN = {
    100: {
        "reverse_auprc": 0.010167334706853206,
        "random_auprc": 0.020753115924147802,
        "masked_pcc": 0.5664336738097762,
        "masked_spearman": 0.3395055309402064,
        "masked_r2": 0.28929768545432377,
    },
    200: {
        "reverse_auprc": 0.02037212470832139,
        "random_auprc": 0.03926800576578986,
        "masked_pcc": 0.6848197324129074,
        "masked_spearman": 0.4704381048630953,
        "masked_r2": 0.44055578822485963,
    },
    500: {
        "reverse_auprc": 0.05185535204252631,
        "random_auprc": 0.09887459775098464,
        "masked_pcc": 0.8433328612998299,
        "masked_spearman": 0.6985699837702515,
        "masked_r2": 0.691457118847044,
    },
}

# Import the metric primitives without requiring the package to be installed.
import metrics  # noqa: E402,F401


class BuildError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# Hashing / serialization helpers


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _content_digest(a: ad.AnnData) -> str:
    """Exact copy of the frozen v1 `_content_digest` recipe.

    The v1 recipe is recorded in `stage2/_frozen/cm-data-v1/README.md` (the
    retired producer has moved there); this copy stays byte-compatible with
    `stage2/_frozen/cm-data-v1/datasets-1/manifest.csv`.

    This is a provenance check only.  `uns["control_mean"]` is read here to
    reproduce the corpus digest recorded in `manifest.csv`; it is never used
    for the Δ baseline (that baseline is recomputed from the DMSO cells in
    `layers["log1p_cp10k"]`).
    """

    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(a.X).tobytes())
    digest.update(np.ascontiguousarray(a.layers["log1p_cp10k"]).tobytes())
    digest.update(a.obs.to_csv(index=False, lineterminator="\n").encode())
    digest.update(a.var.to_csv(lineterminator="\n").encode())
    means = a.uns["control_mean"]
    digest.update(means.to_csv(lineterminator="\n").encode())
    return digest.hexdigest()


def save_npy(path: Path, array: np.ndarray) -> None:
    np.save(path, np.ascontiguousarray(array), allow_pickle=False)


def save_npz(path: Path, arrays: dict[str, np.ndarray]) -> None:
    """Deterministic `.npz` writer.

    `np.savez` stamps entries with the current time, so two identical builds
    differ at the byte level.  A fixed 1980-01-01 timestamp and sorted member
    order make the archive reproducible.
    """

    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_STORED) as archive:
        for name in sorted(arrays):
            buffer = io.BytesIO()
            np.lib.format.write_array(
                buffer, np.ascontiguousarray(arrays[name]), allow_pickle=False
            )
            info = zipfile.ZipInfo(name + ".npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_STORED
            info.external_attr = 0o600 << 16
            archive.writestr(info, buffer.getvalue())


def write_parquet(df: pd.DataFrame, path: Path, schema: pa.Schema) -> None:
    table = pa.Table.from_pandas(df, schema=schema, preserve_index=False)
    table = table.replace_schema_metadata(None)
    pq.write_table(
        table,
        path,
        compression="zstd",
        compression_level=9,
        use_dictionary=False,
        write_statistics=True,
    )


def _str(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value)


def _group_id(drug_key: str, dose: float) -> str:
    return f"{drug_key}|{float(dose):g}"


# ---------------------------------------------------------------------------
# Source loading and verification


def load_gene_set(path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Load the selected 5,000-gene table and check its contract."""

    if not path.exists():
        raise BuildError(f"missing gene set: {path}")
    gene_set = pd.read_csv(path)
    required = {"symbol", "corpus_index", "rank_sum"}
    if not required <= set(gene_set.columns):
        raise BuildError(f"gene set is missing columns {sorted(required - set(gene_set))}")
    if len(gene_set) != N_GENES:
        raise BuildError(f"gene set has {len(gene_set)} rows, expected {N_GENES}")
    if gene_set["symbol"].astype(str).duplicated().any():
        raise BuildError("gene set has duplicate symbols")
    corpus_index = gene_set["corpus_index"].to_numpy(dtype=np.int64)
    if len(np.unique(corpus_index)) != N_GENES:
        raise BuildError("gene set corpus_index is not unique")
    if not np.all(np.diff(corpus_index) > 0):
        raise BuildError("gene set is not in corpus order")
    symbols = gene_set["symbol"].astype(str).to_numpy()
    info = {
        "path": str(path.relative_to(ROOT)),
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
        "n_genes": int(len(gene_set)),
    }
    return symbols, corpus_index, info


def load_source(
    gene_set_path: Path,
) -> tuple[ad.AnnData, pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, dict]:
    if not SOURCE.exists():
        raise BuildError(f"missing source: {SOURCE}")

    selected_symbols, selected_index, gene_set_info = load_gene_set(gene_set_path)

    a = ad.read_h5ad(SOURCE)
    if a.shape != (4_478, SOURCE_GENES):
        raise BuildError(f"unexpected source shape {a.shape}")

    obs = a.obs.copy()
    obs["well"] = obs["well"].astype(str)
    if not (obs.index.astype(str) == obs["well"]).all():
        raise BuildError("obs index is not the well id")
    if obs["well"].duplicated().any():
        raise BuildError("duplicate well ids")

    genes = pd.read_csv(GENES, usecols=["symbol"])
    var_names = np.asarray(a.var_names.astype(str))
    if list(var_names) != genes["symbol"].astype(str).tolist():
        raise BuildError(f"var order does not match {GENES}")
    if len(set(var_names.tolist())) != SOURCE_GENES:
        raise BuildError("duplicate gene symbols")
    if not np.array_equal(var_names[selected_index], selected_symbols):
        raise BuildError("gene set symbols do not match the source var at corpus_index")

    # Provenance is checked on the *full* object before the column subset;
    # `manifest.csv` describes the 18,164-gene corpus view.
    manifest = pd.read_csv(SOURCE_MANIFEST)
    row = manifest[manifest["file"].astype(str) == "temp_env_all-smiles.h5ad"]
    if len(row) != 1:
        raise BuildError("source manifest row missing")
    row = row.iloc[0]
    for key, expected_value in {
        "n_obs": 4_478,
        "n_perturbed": N_TREATED,
        "n_controls": 272,
        "n_drugs": 460,
        "n_genes": SOURCE_GENES,
        "n_control_groups": N_CONTROL_GROUPS,
        "corpus_version": "v1",
    }.items():
        if key in row and str(row[key]) != str(expected_value):
            raise BuildError(f"manifest {key}={row[key]} != {expected_value}")

    content_digest = _content_digest(a)
    if str(row["content_sha256"]) != content_digest:
        raise BuildError(
            "source content_sha256 mismatch: "
            f"{content_digest} != {row['content_sha256']}"
        )

    X_full = np.ascontiguousarray(a.X, dtype=np.int32)
    log1p_full = np.ascontiguousarray(a.layers["log1p_cp10k"], dtype=np.float32)
    if X_full.shape != (4_478, SOURCE_GENES) or log1p_full.shape != (4_478, SOURCE_GENES):
        raise BuildError("matrix shape mismatch")
    if int(X_full.min()) < 0 or not np.isfinite(log1p_full).all():
        raise BuildError("counts/log1p contain invalid values")

    lib_size = X_full.sum(axis=1, dtype=np.int64)
    if not np.array_equal(lib_size, obs["lib_size"].to_numpy(dtype=np.int64)):
        raise BuildError("lib_size does not equal the count row sums")

    # log1p_cp10k must be the documented transform (float32 tolerance).
    check_rows = np.arange(0, X_full.shape[0], max(1, X_full.shape[0] // 64))
    expected = np.log1p(
        X_full[check_rows].astype(np.float64) / lib_size[check_rows, None] * 1e4
    )
    max_diff = float(np.max(np.abs(expected - log1p_full[check_rows].astype(np.float64))))
    if max_diff > 1e-5:
        raise BuildError(f"log1p_cp10k mismatch: max diff {max_diff}")

    X = np.ascontiguousarray(X_full[:, selected_index])
    log1p = np.ascontiguousarray(log1p_full[:, selected_index])
    if X.shape != (4_478, N_GENES) or log1p.shape != (4_478, N_GENES):
        raise BuildError(f"subset matrix shape {X.shape}")

    source_info = {
        "path": str(SOURCE.relative_to(ROOT)),
        "sha256": sha256_file(SOURCE),
        "bytes": SOURCE.stat().st_size,
        "content_sha256": content_digest,
        "shape": list(a.shape),
        "n_genes_used": N_GENES,
        "corpus_version": str(row["corpus_version"]),
        "source_commit": str(row["source_commit"]),
        "gene_set": gene_set_info,
    }
    return a, obs, X, log1p, selected_symbols, source_info


# ---------------------------------------------------------------------------
# Tables


WELLS_SCHEMA = pa.schema(
    [
        ("well", pa.string()),
        ("dataset_id", pa.string()),
        ("compound_name", pa.string()),
        ("drug_key", pa.string()),
        ("inchikey", pa.string()),
        ("smiles", pa.string()),
        ("dose_um", pa.float64()),
        ("dose_source", pa.string()),
        ("time_h", pa.float64()),
        ("is_control", pa.bool_()),
        ("vehicle", pa.string()),
        ("plate", pa.string()),
        ("control_group", pa.string()),
        ("cell_type", pa.string()),
        ("donor", pa.string()),
        ("replicate", pa.int32()),
        ("atc", pa.string()),
        ("moa", pa.string()),
        ("mesh_classes", pa.string()),
        ("pharm_atc_moa", pa.bool_()),
        ("pharm_mesh", pa.bool_()),
        ("annotation_source", pa.string()),
        ("lib_size", pa.int64()),
        ("n_genes_detected", pa.int32()),
        ("delta_row", pa.int32()),
    ]
)

GROUPS_SCHEMA = pa.schema(
    [
        ("group_id", pa.string()),
        ("drug_key", pa.string()),
        ("compound_name", pa.string()),
        ("smiles", pa.string()),
        ("inchikey", pa.string()),
        ("smiles_is_parent", pa.bool_()),
        ("dose_um", pa.float64()),
        ("dose_source", pa.string()),
        ("time_h", pa.float64()),
        ("cell_type", pa.string()),
        ("n_wells", pa.int32()),
        ("n_plates", pa.int32()),
        ("plates", pa.string()),
        ("control_groups", pa.string()),
        ("n_control_groups", pa.int32()),
        ("control_mean_ref", pa.string()),
        ("replicate_min", pa.int32()),
        ("replicate_max", pa.int32()),
        ("pharm_atc_moa", pa.bool_()),
        ("atc", pa.string()),
        ("moa", pa.string()),
    ]
)

SPLITS_SCHEMA = pa.schema(
    [
        ("well", pa.string()),
        ("fold", pa.int8()),
        ("role", pa.string()),
        ("cluster_id", pa.string()),
        ("drug_key", pa.string()),
        ("inchikey", pa.string()),
        ("smiles", pa.string()),
        ("dose_um", pa.float64()),
        ("plate", pa.string()),
        ("control_group", pa.string()),
        ("replicate", pa.int32()),
        ("n_wells", pa.int32()),
        ("n_drugs_in_cluster", pa.int32()),
        ("pharm_atc_moa", pa.bool_()),
        ("test_fold", pa.int8()),
        ("valid_fold", pa.int8()),
    ]
)


def representative_rows(obs_treated: pd.DataFrame) -> pd.DataFrame:
    """One representative structure per drug_key, per ticket 04.

    Prefer the row whose full InChIKey equals the desalted parent key; if none
    matches, take the first row and record `smiles_is_parent = False`.
    """

    records = []
    for drug_key in sorted(obs_treated["drug_key"].astype(str).unique()):
        sub = obs_treated[obs_treated["drug_key"].astype(str) == drug_key]
        parent = sub[sub["inchikey"].astype(str) == drug_key]
        rep = parent.iloc[0] if len(parent) else sub.iloc[0]
        records.append(
            {
                "drug_key": drug_key,
                "compound_name": _str(rep.get("compound_name", "")),
                "smiles": _str(rep.get("smiles", "")),
                "inchikey": _str(rep.get("inchikey", "")),
                "smiles_is_parent": bool(len(parent)),
                "atc": _str(rep.get("atc", "")),
                "moa": _str(rep.get("moa", "")),
                "mesh_classes": _str(rep.get("mesh_classes", "")),
            }
        )
    return pd.DataFrame(records)


def build_control_bundle(
    obs: pd.DataFrame, log1p: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return `(control_groups, cells, mean, wells)` for the 12 DMSO cells."""

    dmso = obs[(obs["is_control"].astype(bool)) & (obs["vehicle"].astype(str) == "DMSO")]
    control_groups = sorted(dmso["control_group"].astype(str).unique())
    if len(control_groups) != N_CONTROL_GROUPS:
        raise BuildError(f"expected {N_CONTROL_GROUPS} DMSO groups, got {len(control_groups)}")

    well_to_row = {str(w): i for i, w in enumerate(obs["well"].astype(str))}
    cells = np.empty((len(control_groups), 12, N_GENES), dtype=np.float32)
    wells = np.empty((len(control_groups), 12), dtype=f"<U{max(len(w) for w in obs['well'].astype(str))}")
    for i, group in enumerate(control_groups):
        sub = dmso[dmso["control_group"].astype(str) == group]
        if len(sub) != 12:
            raise BuildError(f"{group}: expected 12 DMSO cells, got {len(sub)}")
        group_wells = sorted(sub["well"].astype(str).tolist())
        rows = [well_to_row[w] for w in group_wells]
        cells[i] = log1p[rows]
        wells[i] = group_wells

    mean = cells.mean(axis=1, dtype=np.float64).astype(np.float32)
    if not np.isfinite(mean).all():
        raise BuildError("control mean contains non-finite values")
    return (
        np.asarray(control_groups, dtype=str),
        cells,
        mean,
        wells,
    )


def build_delta(
    obs: pd.DataFrame,
    log1p: np.ndarray,
    control_groups: np.ndarray,
    control_mean: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    treated = ~obs["is_control"].astype(bool)
    treated_idx = np.flatnonzero(treated.to_numpy())
    group_to_i = {str(g): i for i, g in enumerate(control_groups)}
    try:
        row_group = np.asarray(
            [group_to_i[str(g)] for g in obs["control_group"].astype(str)], dtype=np.int64
        )
    except KeyError as exc:
        raise BuildError(f"treated well references an unknown control group: {exc}") from exc

    delta = (
        log1p[treated_idx].astype(np.float64)
        - control_mean[row_group[treated_idx]].astype(np.float64)
    ).astype(np.float32)
    if not np.isfinite(delta).all():
        raise BuildError("delta contains non-finite values")

    delta_row = np.full(len(obs), -1, dtype=np.int32)
    delta_row[treated_idx] = np.arange(len(treated_idx), dtype=np.int32)
    return delta, delta_row, treated_idx


def build_group_delta(
    obs_treated: pd.DataFrame, delta: np.ndarray, groups: pd.DataFrame
) -> np.ndarray:
    group_of_row = {
        (str(r.drug_key), float(r.dose_um)): i
        for i, r in groups.reset_index(drop=True).iterrows()
    }
    well_to_delta = {
        str(w): i for i, w in enumerate(obs_treated["well"].astype(str))
    }
    sums = np.zeros((len(groups), N_GENES), dtype=np.float64)
    counts = np.zeros(len(groups), dtype=np.int64)
    for row in obs_treated.itertuples(index=False):
        gi = group_of_row[(str(row.drug_key), float(row.dose_um))]
        sums[gi] += delta[well_to_delta[str(row.well)]].astype(np.float64)
        counts[gi] += 1
    if (counts == 0).any():
        raise BuildError("empty group in group_delta")
    group_delta = (sums / counts[:, None]).astype(np.float32)
    if not np.isfinite(group_delta).all():
        raise BuildError("group_delta contains non-finite values")
    return group_delta


def build_splits(
    obs_treated: pd.DataFrame, reps: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """458-cluster greedy 6:2:2 split, exactly as ticket 06."""

    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    pairs = obs_treated[["drug_key", "inchikey"]].drop_duplicates()
    for drug_key, inchikey in pairs.itertuples(index=False):
        union("K:" + str(drug_key), "B:" + str(inchikey)[:14])

    members: dict[str, list[str]] = {}
    for drug_key in sorted(obs_treated["drug_key"].astype(str).unique()):
        members.setdefault(find("K:" + drug_key), []).append(drug_key)

    cluster_of_drug: dict[str, str] = {}
    for root, keys in members.items():
        cluster_id = "|".join(sorted(keys))
        for key in keys:
            cluster_of_drug[key] = cluster_id

    obs_treated = obs_treated.copy()
    obs_treated["cluster_id"] = obs_treated["drug_key"].astype(str).map(cluster_of_drug)

    cluster = (
        obs_treated.groupby("cluster_id", observed=True)
        .agg(n_wells=("well", "size"), n_drugs=("drug_key", "nunique"))
        .reset_index()
    )
    cluster["order_key"] = cluster["cluster_id"].map(
        lambda cid: hashlib.sha256(f"{SEED}:{cid}".encode()).hexdigest()
    )
    cluster = cluster.sort_values(
        ["n_wells", "order_key"], ascending=[False, True], kind="mergesort"
    ).reset_index(drop=True)

    load_wells = [0] * N_FOLDS
    load_drugs = [0] * N_FOLDS
    fold_of: dict[str, int] = {}
    for row in cluster.itertuples(index=False):
        candidates = [
            k for k in range(N_FOLDS) if load_drugs[k] + int(row.n_drugs) <= 92
        ]
        if not candidates:
            raise BuildError(f"no fold can take cluster {row.cluster_id}")
        k = min(candidates, key=lambda k: (load_wells[k], load_drugs[k], k))
        load_wells[k] += int(row.n_wells)
        load_drugs[k] += int(row.n_drugs)
        fold_of[str(row.cluster_id)] = k

    obs_treated["test_fold"] = obs_treated["cluster_id"].map(fold_of).astype(np.int8)
    obs_treated["valid_fold"] = ((obs_treated["test_fold"].astype(int) + 4) % 5).astype(np.int8)

    rep = reps.set_index("drug_key")
    frames = []
    for k in range(N_FOLDS):
        fold = obs_treated.copy()
        fold["fold"] = np.int8(k)
        fold["role"] = np.where(
            fold["test_fold"] == k, "test", np.where(fold["valid_fold"] == k, "valid", "train")
        )
        fold = fold.merge(
            rep[["smiles", "inchikey"]].rename(
                columns={"smiles": "rep_smiles", "inchikey": "rep_inchikey"}
            ),
            left_on="drug_key",
            right_index=True,
            how="left",
            validate="many_to_one",
        )
        fold["smiles"] = fold["rep_smiles"].fillna("").astype(str)
        fold["inchikey"] = fold["rep_inchikey"].fillna("").astype(str)
        frames.append(
            pd.DataFrame(
                {
                    "well": fold["well"].astype(str),
                    "fold": fold["fold"].astype(np.int8),
                    "role": fold["role"].astype(str),
                    "cluster_id": fold["cluster_id"].astype(str),
                    "drug_key": fold["drug_key"].astype(str),
                    "inchikey": fold["inchikey"].astype(str),
                    "smiles": fold["smiles"].astype(str),
                    "dose_um": fold["dose_um"].astype(np.float64),
                    "plate": fold["plate"].astype(str),
                    "control_group": fold["control_group"].astype(str),
                    "replicate": fold["replicate"].astype(np.int32),
                    "n_wells": fold["cluster_id"].map(
                        cluster.set_index("cluster_id")["n_wells"]
                    ).astype(np.int32),
                    "n_drugs_in_cluster": fold["cluster_id"].map(
                        cluster.set_index("cluster_id")["n_drugs"]
                    ).astype(np.int32),
                    "pharm_atc_moa": fold["pharm_atc_moa"].astype(bool),
                    "test_fold": fold["test_fold"].astype(np.int8),
                    "valid_fold": fold["valid_fold"].astype(np.int8),
                }
            )
        )
    splits = pd.concat(frames, ignore_index=True)
    splits = splits.sort_values(["fold", "well"], kind="mergesort").reset_index(drop=True)
    return splits, cluster


def build_groups(obs_treated: pd.DataFrame, reps: pd.DataFrame) -> pd.DataFrame:
    rep = reps.set_index("drug_key")
    rows = []
    for (drug_key, dose), sub in obs_treated.groupby(
        ["drug_key", "dose_um"], observed=True, sort=True
    ):
        drug_key = str(drug_key)
        dose = float(dose)
        r = rep.loc[drug_key]
        plates = sorted(sub["plate"].astype(str).unique())
        cgs = sorted(sub["control_group"].astype(str).unique())
        dose_sources = sorted(sub["dose_source"].astype(str).unique())
        rows.append(
            {
                "group_id": _group_id(drug_key, dose),
                "drug_key": drug_key,
                "compound_name": _str(r["compound_name"]),
                "smiles": _str(r["smiles"]),
                "inchikey": _str(r["inchikey"]),
                "smiles_is_parent": bool(r["smiles_is_parent"]),
                "dose_um": dose,
                "dose_source": ";".join(dose_sources),
                "time_h": float(sub["time_h"].iloc[0]),
                "cell_type": _str(sub["cell_type"].iloc[0]),
                "n_wells": int(len(sub)),
                "n_plates": int(len(plates)),
                "plates": ";".join(plates),
                "control_groups": ";".join(cgs),
                "n_control_groups": int(len(cgs)),
                "control_mean_ref": ";".join(cgs),
                "replicate_min": int(sub["replicate"].min()),
                "replicate_max": int(sub["replicate"].max()),
                "pharm_atc_moa": bool(sub["pharm_atc_moa"].any()),
                "atc": _str(r["atc"]),
                "moa": _str(r["moa"]),
            }
        )
    groups = pd.DataFrame(rows).sort_values(
        ["drug_key", "dose_um"], kind="mergesort"
    ).reset_index(drop=True)
    return groups


# ---------------------------------------------------------------------------
# True sets


def build_true_sets(
    group_delta: np.ndarray, group_ids: np.ndarray
) -> dict[int, dict[str, np.ndarray]]:
    out: dict[int, dict[str, np.ndarray]] = {}
    n = group_delta.shape[1]
    for k in K_LIST:
        labels = np.zeros((len(group_ids), n), dtype=np.int8)
        n_up = np.zeros(len(group_ids), dtype=np.int32)
        n_down = np.zeros(len(group_ids), dtype=np.int32)
        tie_up = np.zeros(len(group_ids), dtype=np.int32)
        tie_down = np.zeros(len(group_ids), dtype=np.int32)
        for i, delta in enumerate(group_delta):
            up, down = metrics.true_set_indices(delta, k)
            labels[i, up] = 1
            labels[i, down] = -1
            n_up[i] = len(up)
            n_down[i] = len(down)

            desc = np.sort(delta.astype(np.float64))[::-1]
            tie_up[i] = int(np.sum(delta == desc[k - 1])) - 1
            asc = np.sort(delta.astype(np.float64))
            tie_down[i] = int(np.sum(delta == asc[k - 1])) - 1
        if not np.array_equal(np.count_nonzero(labels == 1, axis=1), n_up):
            raise BuildError(f"k={k}: up label count mismatch")
        if not np.array_equal(np.count_nonzero(labels == -1, axis=1), n_down):
            raise BuildError(f"k={k}: down label count mismatch")
        out[k] = {
            "labels": labels,
            "group_ids": group_ids.astype(str),
            "n_up": n_up,
            "n_down": n_down,
            "boundary_ties_up": tie_up,
            "boundary_ties_down": tie_down,
            "k": np.asarray(k, dtype=np.int32),
        }
    return out


# ---------------------------------------------------------------------------
# Gates


def _role_counts(splits: pd.DataFrame, fold: int, role: str) -> tuple[int, int, int]:
    sub = splits[(splits["fold"] == fold) & (splits["role"] == role)]
    return (
        int(sub["drug_key"].nunique()),
        int(sub.groupby(["drug_key", "dose_um"], observed=True).ngroups),
        int(len(sub)),
    )


def run_gates(
    obs: pd.DataFrame,
    X: np.ndarray,
    log1p: np.ndarray,
    var_names: np.ndarray,
    groups: pd.DataFrame,
    splits: pd.DataFrame,
    control_groups: np.ndarray,
    cells: np.ndarray,
    control_mean: np.ndarray,
    delta: np.ndarray,
    delta_row: np.ndarray,
    group_delta: np.ndarray,
    true_sets: dict[int, dict[str, np.ndarray]],
    gene_set: np.ndarray | None = None,
) -> dict[str, bool]:
    gates: dict[str, bool] = {}
    treated = obs[~obs["is_control"].astype(bool)].copy()
    all_wells = set(treated["well"].astype(str))

    gates["01_coverage"] = (
        len(splits) == N_TREATED * N_FOLDS
        and not splits.duplicated(["well", "fold"]).any()
    )

    partition_ok = True
    for k in range(N_FOLDS):
        fold = splits[splits["fold"] == k]
        roles = {r: set(fold.loc[fold["role"] == r, "well"].astype(str)) for r in ("train", "valid", "test")}
        partition_ok &= (roles["train"] | roles["valid"] | roles["test"]) == all_wells
        partition_ok &= not (roles["train"] & roles["valid"])
        partition_ok &= not (roles["train"] & roles["test"])
        partition_ok &= not (roles["valid"] & roles["test"])
    gates["02_partition"] = bool(partition_ok)

    counts_ok = True
    for k in range(N_FOLDS):
        counts_ok &= _role_counts(splits, k, "test")[:2] == (92, 276)
        counts_ok &= _role_counts(splits, k, "valid")[:2] == (92, 276)
        counts_ok &= _role_counts(splits, k, "train")[:2] == (276, 828)
    gates["03_counts"] = bool(counts_ok and len(groups) == N_GROUPS)

    cold_ok = True
    for k in range(N_FOLDS):
        fold = splits[splits["fold"] == k]
        test = set(fold.loc[fold["role"] == "test", "drug_key"].astype(str))
        train = set(fold.loc[fold["role"] == "train", "drug_key"].astype(str))
        valid = set(fold.loc[fold["role"] == "valid", "drug_key"].astype(str))
        cold_ok &= not (test & (train | valid)) and not (valid & train)
    gates["04_cold_drug"] = bool(cold_ok)

    identity_ok = splits.groupby("cluster_id")["test_fold"].nunique().max() == 1
    identity_ok &= splits.groupby("cluster_id")["valid_fold"].nunique().max() == 1
    identity_ok &= splits.groupby("drug_key")["cluster_id"].nunique().max() == 1
    identity_ok &= treated.groupby("inchikey", observed=True)["drug_key"].nunique().max() == 1
    gates["05_identity"] = bool(identity_ok)

    test_once = splits[splits["role"] == "test"].groupby("drug_key")["fold"].nunique()
    valid_once = splits[splits["role"] == "valid"].groupby("drug_key")["fold"].nunique()
    train_three = splits[splits["role"] == "train"].groupby("drug_key")["fold"].nunique()
    gates["06_global_roles"] = bool(
        len(test_once) == 460
        and len(valid_once) == 460
        and len(train_three) == 460
        and (test_once == 1).all()
        and (valid_once == 1).all()
        and (train_three == 3).all()
    )

    cg = set(map(str, control_groups))
    controls_ok = len(cg) == N_CONTROL_GROUPS
    controls_ok &= cells.shape == (N_CONTROL_GROUPS, 12, N_GENES)
    controls_ok &= control_mean.shape == (N_CONTROL_GROUPS, N_GENES)
    controls_ok &= set(treated["control_group"].astype(str)) <= cg
    pfoa = splits[
        (splits["role"] == "test")
        & (splits["drug_key"].astype(str) == "SNGREZUHAYWORS-UHFFFAOYSA-N")
    ]
    controls_ok &= bool(len(pfoa) and pfoa.groupby("fold")["plate"].nunique().eq(2).all())
    gates["07_controls"] = bool(controls_ok)

    support_ok = True
    for k in range(N_FOLDS):
        train = splits[(splits["fold"] == k) & (splits["role"] == "train")]
        support_ok &= bool(
            (train.groupby(["drug_key", "dose_um"], observed=True).size() >= 1).all()
        )
    gates["08_min_support"] = bool(support_ok)

    test_wells = splits[splits["role"] == "test"].groupby("fold").size()
    balance_ok = int(test_wells.max() - test_wells.min()) <= 30
    balance_ok &= all(
        splits[(splits["fold"] == k) & (splits["role"] == "test")]["plate"].nunique() == 16
        for k in range(N_FOLDS)
    )
    balance_ok &= all(
        splits[(splits["fold"] == k) & (splits["role"] == "test")]
        .drop_duplicates("drug_key")["pharm_atc_moa"]
        .sum()
        >= 1
        for k in range(N_FOLDS)
    )
    gates["09_balance"] = bool(balance_ok)

    # Gate 10 is the deterministic rerun check; the in-build part asserts the
    # source/gene ordering and that no 978 path is present.
    if gene_set is None:
        gene_set = pd.read_csv(DEFAULT_GENE_SET, usecols=["symbol"])["symbol"].astype(str).to_numpy()
    gates["10_determinism_inputs"] = bool(
        len(var_names) == N_GENES
        and len(set(map(str, var_names))) == N_GENES
        and list(map(str, var_names)) == list(map(str, gene_set))
        and not any("978" in str(p) for p in (SOURCE, GENES, DEFAULT_GENE_SET))
    )

    # Additional integrity checks that must also hold.
    gates["11_matrix_integrity"] = bool(
        X.shape == (4_478, N_GENES)
        and log1p.shape == (4_478, N_GENES)
        and delta.shape == (N_TREATED, N_GENES)
        and delta_row.shape == (4_478,)
        and group_delta.shape == (N_GROUPS, N_GENES)
        and np.isfinite(delta).all()
        and np.isfinite(group_delta).all()
    )
    gates["12_true_sets"] = bool(
        all(
            v["labels"].shape == (N_GROUPS, N_GENES)
            and np.array_equal(v["n_up"], np.full(N_GROUPS, k, dtype=np.int32))
            and np.array_equal(v["n_down"], np.full(N_GROUPS, k, dtype=np.int32))
            for k, v in true_sets.items()
        )
        and all(
            np.array_equal(np.count_nonzero(v["labels"] == 1, axis=1), v["n_up"])
            and np.array_equal(np.count_nonzero(v["labels"] == -1, axis=1), v["n_down"])
            for v in true_sets.values()
        )
    )

    failed = [name for name, ok in gates.items() if not ok]
    if failed:
        raise BuildError(f"gates failed: {failed}")
    return gates


# ---------------------------------------------------------------------------
# Metric sanity cases


def run_metric_sanity(
    group_delta: np.ndarray,
    group_ids: np.ndarray,
    true_sets: dict[int, dict[str, np.ndarray]],
    first_group_id: str,
) -> dict[str, dict]:
    """Five fixed cases against the first treated well's group (ticket 07)."""

    matches = np.flatnonzero(group_ids == first_group_id)
    if not len(matches):
        raise BuildError(f"first treated group {first_group_id} not in groups.parquet")
    i = int(matches[0])
    true = group_delta[i].astype(np.float64)
    n = true.size
    report: dict[str, dict] = {}

    for k in K_LIST:
        labels = true_sets[k]["labels"][i]
        up = np.flatnonzero(labels == 1)
        down = np.flatnonzero(labels == -1)
        up_ref, down_ref = metrics.true_set_indices(true, k)
        if not np.array_equal(np.sort(up), np.sort(up_ref)):
            raise BuildError(f"k={k}: persisted up set != recomputed set")
        if not np.array_equal(np.sort(down), np.sort(down_ref)):
            raise BuildError(f"k={k}: persisted down set != recomputed set")

        prevalence = k / n
        cases: dict[str, dict] = {}

        perfect = metrics.group_metrics(true, true, k)
        if not (
            math.isclose(perfect["auprc"], 1.0, abs_tol=1e-12)
            and math.isclose(perfect["pcc"], 1.0, abs_tol=1e-12)
            and math.isclose(perfect["spearman"], 1.0, abs_tol=1e-12)
            and math.isclose(perfect["r2"], 1.0, abs_tol=1e-12)
            and math.isclose(perfect["rmse"], 0.0, abs_tol=1e-12)
        ):
            raise BuildError(f"k={k}: perfect-prediction sanity failed: {perfect}")
        cases["perfect"] = perfect

        constant = metrics.group_metrics(true, np.zeros_like(true), k)
        if not math.isclose(constant["auprc"], prevalence, abs_tol=1e-12):
            raise BuildError(f"k={k}: constant AUPRC {constant['auprc']} != {prevalence}")
        if not (math.isnan(constant["pcc"]) and math.isnan(constant["spearman"])):
            raise BuildError(f"k={k}: constant prediction should have NaN PCC/Spearman")
        expected_r2 = 1.0 - float(np.mean(true ** 2)) / float(np.var(true))
        expected_rmse = float(np.sqrt(np.mean(true ** 2)))
        if not math.isclose(constant["r2"], expected_r2, abs_tol=1e-12):
            raise BuildError(f"k={k}: constant R2 {constant['r2']} != {expected_r2}")
        if not math.isclose(constant["rmse"], expected_rmse, abs_tol=1e-12):
            raise BuildError(f"k={k}: constant RMSE {constant['rmse']} != {expected_rmse}")
        cases["constant"] = constant

        reverse = metrics.group_metrics(true, -true, k)
        if not math.isclose(reverse["auprc_up"], reverse["auprc_down"], abs_tol=1e-12):
            raise BuildError(f"k={k}: reverse up/down AP differ")
        if not math.isclose(
            reverse["auprc"], SANITY_GOLDEN[k]["reverse_auprc"], abs_tol=1e-9
        ):
            raise BuildError(
                f"k={k}: reverse floor {reverse['auprc']} != 5k golden"
            )
        if not (
            math.isclose(reverse["pcc"], -1.0, abs_tol=1e-12)
            and math.isclose(reverse["spearman"], -1.0, abs_tol=1e-12)
        ):
            raise BuildError(f"k={k}: reverse PCC/Spearman not -1")
        cases["reverse"] = reverse

        rng = np.random.default_rng(20261021)
        random_pred = rng.normal(size=n)
        random = metrics.group_metrics(true, random_pred, k)
        if not math.isclose(
            random["auprc"], SANITY_GOLDEN[k]["random_auprc"], abs_tol=1e-9
        ):
            raise BuildError(f"k={k}: random AUPRC {random['auprc']} != 5k golden")
        cases["random"] = random

        masked_pred = np.where(labels != 0, true, 0.0)
        masked = metrics.group_metrics(true, masked_pred, k)
        if not math.isclose(masked["auprc"], 1.0, abs_tol=1e-12):
            raise BuildError(f"k={k}: top/bottom-only AUPRC {masked['auprc']} != 1")
        for metric in ("pcc", "spearman", "r2"):
            if not math.isclose(
                masked[metric], SANITY_GOLDEN[k][f"masked_{metric}"], abs_tol=1e-9
            ):
                raise BuildError(
                    f"k={k}: top/bottom-only {metric} {masked[metric]} != 5k golden"
                )
        cases["top_bottom_only"] = masked

        report[str(k)] = {
            "group_id": str(first_group_id),
            "prevalence": prevalence,
            "cases": cases,
        }
    return report


# ---------------------------------------------------------------------------
# Manifest / README


def build_manifest(
    out: Path,
    source_info: dict,
    var_names: np.ndarray,
    groups: pd.DataFrame,
    splits: pd.DataFrame,
    true_sets: dict[int, dict[str, np.ndarray]],
    gates: dict[str, bool],
    sanity: dict,
    source_parity: dict | None,
) -> dict:
    artifacts = {}
    for path in sorted(out.iterdir()):
        if path.name in {"manifest.json", "README.md"}:
            continue
        if path.suffix == ".npy":
            arr = np.load(path, mmap_mode="r")
            artifacts[path.name] = {
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "shape": list(arr.shape),
                "dtype": str(arr.dtype),
            }
        elif path.suffix == ".npz":
            with np.load(path, allow_pickle=False) as npz:
                arrays = {
                    name: {"shape": list(npz[name].shape), "dtype": str(npz[name].dtype)}
                    for name in npz.files
                }
                shape = "; ".join(
                    f"{name} {'×'.join(map(str, info['shape']))} {info['dtype']}"
                    for name, info in arrays.items()
                )
                artifacts[path.name] = {
                    "sha256": sha256_file(path),
                    "bytes": path.stat().st_size,
                    "shape": shape,
                    "arrays": arrays,
                }
        elif path.suffix == ".parquet":
            artifacts[path.name] = {
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "rows": int(len(pq.read_table(path))),
            }
        else:
            artifacts[path.name] = {
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }

    split_counts = []
    for k in range(N_FOLDS):
        row = {"fold": k}
        for role in ("test", "valid", "train"):
            sub = splits[(splits["fold"] == k) & (splits["role"] == role)]
            row[f"{role}_drugs"] = int(sub["drug_key"].nunique())
            row[f"{role}_groups"] = int(sub.groupby(["drug_key", "dose_um"], observed=True).ngroups)
            row[f"{role}_wells"] = int(len(sub))
        row["test_plates"] = int(
            splits[(splits["fold"] == k) & (splits["role"] == "test")]["plate"].nunique()
        )
        row["test_atc_moa_drugs"] = int(
            splits[(splits["fold"] == k) & (splits["role"] == "test")]
            .drop_duplicates("drug_key")["pharm_atc_moa"]
            .sum()
        )
        split_counts.append(row)

    return {
        "corpus_version": source_info["corpus_version"],
        "benchmark": "temp_env_all_5k",
        "gene_set": source_info["gene_set"],
        "source": source_info,
        "gene_order_sha256": hashlib.sha256("\n".join(map(str, var_names)).encode()).hexdigest(),
        "seed": SEED,
        "algorithm": "458 connectivity-block clusters; n_wells desc then sha256(seed:cluster) asc; capacity-constrained greedy into fold minimizing (treated wells, drugs, fold); valid=(test+4)%5",
        "build_command": "python stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/build_data.py",
        "n_genes": N_GENES,
        "n_groups": N_GROUPS,
        "n_treated": N_TREATED,
        "groups_sha256": sha256_file(out / "groups.parquet"),
        "splits_sha256": sha256_file(out / "splits.parquet"),
        "split_counts": split_counts,
        "true_sets": {
            str(k): {
                "prevalence": k / N_GENES,
                "max_boundary_tie_up": int(np.max(v["boundary_ties_up"])),
                "max_boundary_tie_down": int(np.max(v["boundary_ties_down"])),
                "min_n_up": int(np.min(v["n_up"])),
                "min_n_down": int(np.min(v["n_down"])),
            }
            for k, v in true_sets.items()
        },
        "sanity": sanity,
        "gates": gates,
        "source_parity": source_parity,
        "artifacts": artifacts,
        "environment": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "pyarrow": pa.__version__,
            "anndata": _pkg_version("anndata"),
        },
    }


def write_readme(out: Path, manifest: dict) -> None:
    rows = []
    for name, info in manifest["artifacts"].items():
        if "shape" in info:
            if isinstance(info["shape"], str):
                shape = info["shape"]
            else:
                shape = "×".join(map(str, info["shape"]))
        elif "rows" in info:
            shape = f"{info['rows']} rows"
        else:
            shape = ""
        rows.append(
            f"| `{name}` | {shape} | `{info['sha256']}` | {info['bytes']:,} |"
        )
    split_rows = []
    for row in manifest["split_counts"]:
        split_rows.append(
            "| {fold} | {test_drugs} / {test_groups} / {test_wells} | "
            "{valid_drugs} / {valid_groups} / {valid_wells} | "
            "{train_drugs} / {train_groups} / {train_wells} | {test_plates} | {test_atc_moa_drugs} |".format(**row)
        )
    text = f"""# Benchmark data — Stage 2 `temp_env_all_5k`

This directory is the frozen input of the 5,000-gene sensitivity benchmark.
Its wells, Δ, controls, groups and splits are the same as the executed
`temp_env_all_18k` archive; only the gene axis differs.  The axis is the 5,000
genes with the smallest sum of per-source HVG ranks, selected by
`scripts/build_gene_set.py` (provenance in `gene_set_5k.manifest.json`).

Built by `python stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/build_data.py` from
`{manifest['source']['path']}` (sha256 `{manifest['source']['sha256']}`,
{manifest['source']['bytes']:,} bytes; corpus `{manifest['source']['corpus_version']}`,
source commit `{manifest['source']['source_commit']}`), column-subset to the
gene set `{manifest['gene_set']['path']}` (sha256
`{manifest['gene_set']['sha256']}`).

This directory is the frozen input of the benchmark.  Every adapter and the
evaluator read these files; nothing here is recomputed during a model run.

## Gene axis

- The selection is the sum of three per-source HVG ranks: `dtoxs_cm`,
  `screen_ctox` and `temp_env_all-smiles`.  Each source is ranked with
  `scanpy.pp.highly_variable_genes(flavor="seurat", n_top_genes=None)` on
  `layers["log1p_cp10k"]` over all wells, no `batch_key`.
- Genes with mean in `[0.0125, 3]` are ranked 1..M by descending
  `dispersions_norm`; every other gene (including zero-mean/NaN) gets the worst
  rank 18,165.  The 5,000 smallest rank sums are kept, ties by corpus order.
- `gene_set_5k.csv` carries the symbol, corpus index, per-source means and
  dispersions, per-source rank and the rank sum.

## Δ and controls

- `delta_log1p_cp10k.npy` is `log1p_cp10k(treated well) − mean(log1p_cp10k of
  the 12 DMSO cells in the same control_group)`, float32, treated wells only.
- `wells.parquet.delta_row` maps a treated row to its Δ row (`-1` for controls),
  so any adapter can recover the raw input and the matched baseline.
- `controls.npz` holds the 16 DMSO control groups: `control_groups`, `cells`
  (16×12×5,000 float32), `mean` (16×5,000 float32), and `wells` (16×12).
  The corpus `uns["control_mean"]` is **not** used; it pools DMSO + Media.
  The only access to it is the provenance hash check below.

## Groups and true sets

- One group = `(drug_key, dose_um)`; `group_id` = `drug_key|dose` (for example
  `SNGREZUHAYWORS-UHFFFAOYSA-N|10`). 1,380 groups, all with 3–6 wells.
- `group_delta.npy` is the true group-mean Δ, row order = `groups.parquet`.
- `true_sets_k{{100,200,500}}.npz` store `labels` (`int8` 1,380×5,000; +1 up,
  −1 down, 0 neutral), `group_ids`, `n_up`, `n_down`,
  `boundary_ties_up/down`, and `k`. Ties are broken by gene-axis order
  (stable argsort); each direction is exactly *k* genes on this axis.
  Random prevalence is 2.0% / 4.0% / 10.0% for k = 100 / 200 / 500, five times
  the 18,164-gene prevalence.

## Splits

`{manifest['algorithm']}`.

| fold | test (drugs/groups/wells) | valid (drugs/groups/wells) | train (drugs/groups/wells) | test plates | test atc-moa drugs |
|---|---|---|---|---|---|
{chr(10).join(split_rows)}

`{manifest['splits_sha256']}` is the split file hash.  Controls are never held
out and are not rows in `splits.parquet`; every test drug's Δ uses its own
plate's DMSO cells.  In `splits.parquet`, `n_wells` and `n_drugs_in_cluster`
are **cluster-level** metadata; per-group sizes are in `groups.parquet`.

## Artifacts

| file | shape | sha256 | bytes |
|---|---|---|---|
{chr(10).join(rows)}

## Build and verification

```bash
python stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/build_gene_set.py
python stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/build_data.py
python stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/build_data.py --check
python -m unittest discover -s stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/tests -v
```

`--check` verifies every artifact's sha256 and size against `manifest.json`;
the unit tests re-run the twelve gates (the ten split gates plus
matrix/true-set integrity) and the five metric sanity cases
without rebuilding.  `build_data.py` also verifies that every rebuilt array is
exactly the 18k array at the same gene columns, and that `wells.parquet`,
`groups.parquet` and `splits.parquet` are byte-identical to the 18k archive.
A second build to a fresh `--out` directory is byte-identical: `.npy` files are
plain arrays, `.npz` archives use a fixed 1980 timestamp, and Parquet is
written with a fixed schema and no pandas metadata.

## Notes

- The corpus actually has **two** `drug_key`s with multiple structures/names
  (PFOA `SNGREZUHAYWORS-UHFFFAOYSA-N`, Tamoxifen `NKANXQFJJICGDU-QPLCGJKRSA-N`),
  not the seven stated in the older ticket text.  Both representative rows are
  the parent (free) form; 27 keys have no exact parent row and use the first
  row with `smiles_is_parent = false`.
- The 21-well cluster is Dieldrin (12 wells) + Endrin (9); the 18-well cluster
  is Allethrin (9) + S-Bioallethrin (9).  The older ticket parenthetical
  swapped these names; the cluster-size distribution (439×9, 15×12, 3×18,
  1×21) and the split are unchanged.
- `counts.npy` is `int32`; `log1p_cp10k.npy`, Δ, and group Δ are `float32`.
  Gene order is the 5,000-row subset of `corpus/v1/genes.csv`; the order hash
  is in `manifest.json` and the full selection is in `gene_set_5k.csv`.
- The HVG selection is unsupervised and uses all wells, including wells of
  drugs that are test in some fold.  This matches the authors' pipelines and
  is fixed across folds; no fold-specific gene axis exists.
- No 978-landmark file is read.  No Stage 1 data or code is imported.
- Environment: Python {manifest['environment']['python']}, numpy
  {manifest['environment']['numpy']}, pandas {manifest['environment']['pandas']},
  pyarrow {manifest['environment']['pyarrow']}, anndata
  {manifest['environment']['anndata']}.
"""
    (out / "README.md").write_text(text, encoding="utf-8")


# ---------------------------------------------------------------------------
# Main


def _array_parity(
    new_path: Path, source_path: Path, index: np.ndarray | None, name: str
) -> int:
    """Assert `new[:, j] == source[:, index[j]]` and return rows checked."""

    new = np.load(new_path, mmap_mode="r")
    source = np.load(source_path, mmap_mode="r")
    if new.shape[1] != N_GENES:
        raise BuildError(f"{name}: new width {new.shape[1]} != {N_GENES}")
    if index is not None and source.shape[1] != SOURCE_GENES:
        raise BuildError(f"{name}: source width {source.shape[1]}")
    if new.shape[0] != source.shape[0]:
        raise BuildError(f"{name}: row mismatch {new.shape[0]} != {source.shape[0]}")
    for start in range(0, new.shape[0], 256):
        stop = min(start + 256, new.shape[0])
        got = np.asarray(new[start:stop])
        want = np.asarray(source[start:stop])
        if index is not None:
            want = want[:, index]
        if not np.array_equal(got, want):
            raise BuildError(f"{name}: rows {start}:{stop} differ from the 18k source")
    return int(new.shape[0])


def verify_source_parity(
    out: Path, source_data: Path, selected_index: np.ndarray
) -> dict:
    """Prove the 5k build is the 18k build on the selected gene columns."""

    if not (source_data / "manifest.json").exists():
        raise BuildError(f"missing 18k source data: {source_data}")
    checks = {}
    for name in ("counts.npy", "log1p_cp10k.npy", "delta_log1p_cp10k.npy", "group_delta.npy"):
        rows = _array_parity(out / name, source_data / name, selected_index, name)
        checks[name] = {"rows": rows, "columns": N_GENES, "match": True}

    with np.load(out / "controls.npz", allow_pickle=False) as new, np.load(
        source_data / "controls.npz", allow_pickle=False
    ) as source:
        if not np.array_equal(new["control_groups"], source["control_groups"]):
            raise BuildError("controls: control_groups differ")
        if not np.array_equal(new["wells"], source["wells"]):
            raise BuildError("controls: wells differ")
        if not np.array_equal(new["cells"], source["cells"][:, :, selected_index]):
            raise BuildError("controls: cells differ from the 18k subset")
        if not np.array_equal(new["mean"], source["mean"][:, selected_index]):
            raise BuildError("controls: mean differs from the 18k subset")
    checks["controls.npz"] = {"groups": N_CONTROL_GROUPS, "match": True}

    for name in ("wells.parquet", "groups.parquet", "splits.parquet"):
        new = pd.read_parquet(out / name)
        source = pd.read_parquet(source_data / name)
        try:
            pd.testing.assert_frame_equal(new, source, check_dtype=True)
        except AssertionError as exc:
            raise BuildError(f"{name} differs from the 18k table: {exc}") from exc
        checks[name] = {"rows": int(len(new)), "match": True}
    return checks


def build(out: Path, gene_set_path: Path, source_data: Path | None) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    for artifact in (
        gene_set_path,
        gene_set_path.with_suffix(".manifest.json"),
        gene_set_path.with_name(gene_set_path.stem + "_ranks.csv"),
    ):
        target = out / artifact.name
        if artifact.resolve() != target.resolve():
            target.write_bytes(artifact.read_bytes())
    a, obs, X, log1p, var_names, source_info = load_source(gene_set_path)
    treated = obs[~obs["is_control"].astype(bool)].copy()
    selected_index = pd.read_csv(gene_set_path)["corpus_index"].to_numpy(dtype=np.int64)
    del a

    reps = representative_rows(treated)
    groups = build_groups(treated, reps)

    control_groups, cells, control_mean, control_wells = build_control_bundle(obs, log1p)
    delta, delta_row, treated_idx = build_delta(obs, log1p, control_groups, control_mean)
    group_delta = build_group_delta(treated, delta, groups)

    group_ids = groups["group_id"].astype(str).to_numpy()
    true_sets = build_true_sets(group_delta, group_ids)

    splits, cluster = build_splits(treated, reps)

    wells = pd.DataFrame(
        {
            "well": obs["well"].astype(str),
            "dataset_id": obs["dataset_id"].astype(str),
            "compound_name": obs["compound_name"].map(_str),
            "drug_key": obs["drug_key"].map(_str),
            "inchikey": obs["inchikey"].map(_str),
            "smiles": obs["smiles"].map(_str),
            "dose_um": obs["dose_um"].astype(np.float64),
            "dose_source": obs["dose_source"].map(_str),
            "time_h": obs["time_h"].astype(np.float64),
            "is_control": obs["is_control"].astype(bool),
            "vehicle": obs["vehicle"].map(_str),
            "plate": obs["plate"].astype(str),
            "control_group": obs["control_group"].astype(str),
            "cell_type": obs["cell_type"].astype(str),
            "donor": obs["donor"].map(_str),
            "replicate": obs["replicate"].astype(np.int32),
            "atc": obs["atc"].map(_str),
            "moa": obs["moa"].map(_str),
            "mesh_classes": obs["mesh_classes"].map(_str),
            "pharm_atc_moa": obs["pharm_atc_moa"].astype(bool),
            "pharm_mesh": obs["pharm_mesh"].astype(bool),
            "annotation_source": obs["annotation_source"].map(_str),
            "lib_size": obs["lib_size"].astype(np.int64),
            "n_genes_detected": obs["n_genes_detected"].astype(np.int32),
            "delta_row": delta_row,
        }
    )
    if len(wells) != 4_478 or len(groups) != N_GROUPS:
        raise BuildError(f"table sizes: wells={len(wells)} groups={len(groups)}")
    if len(splits) != N_TREATED * N_FOLDS:
        raise BuildError(f"splits rows {len(splits)}")

    # Run the fail-closed gates before writing any manifest.
    gates = run_gates(
        obs, X, log1p, var_names, groups, splits, control_groups, cells,
        control_mean, delta, delta_row, group_delta, true_sets,
        gene_set=var_names,
    )
    sanity = run_metric_sanity(
        group_delta,
        group_ids,
        true_sets,
        str(group_ids[0]),
    )

    save_npy(out / "counts.npy", X)
    save_npy(out / "log1p_cp10k.npy", log1p)
    save_npy(out / "delta_log1p_cp10k.npy", delta)
    save_npy(out / "group_delta.npy", group_delta)
    save_npz(out / "controls.npz", {
        "control_groups": control_groups,
        "cells": cells,
        "mean": control_mean,
        "wells": control_wells,
    })
    for k, arrays in true_sets.items():
        save_npz(out / f"true_sets_k{k}.npz", arrays)

    write_parquet(wells, out / "wells.parquet", WELLS_SCHEMA)
    write_parquet(groups, out / "groups.parquet", GROUPS_SCHEMA)
    write_parquet(splits, out / "splits.parquet", SPLITS_SCHEMA)

    # splits.md is written before the manifest so its hash is recorded.
    split_rows = []
    for k in range(N_FOLDS):
        row = {"fold": k}
        for role in ("test", "valid", "train"):
            sub = splits[(splits["fold"] == k) & (splits["role"] == role)]
            row[f"{role}_drugs"] = int(sub["drug_key"].nunique())
            row[f"{role}_groups"] = int(sub.groupby(["drug_key", "dose_um"], observed=True).ngroups)
            row[f"{role}_wells"] = int(len(sub))
        row["test_plates"] = int(
            splits[(splits["fold"] == k) & (splits["role"] == "test")]["plate"].nunique()
        )
        row["test_atc_moa_drugs"] = int(
            splits[(splits["fold"] == k) & (splits["role"] == "test")]
            .drop_duplicates("drug_key")["pharm_atc_moa"].sum()
        )
        split_rows.append(row)
    md = [
        "# Cold-drug split — 5 folds (ticket 06)",
        "",
        f"Seed: `{SEED}`.  Assignment: {json.dumps(SPLIT_ALGORITHM)}",
        "",
        "| fold | test drugs | test groups | test wells | valid drugs | valid groups | valid wells | train drugs | train groups | train wells | test plates | test atc-moa drugs |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for row in split_rows:
        md.append(
            "| {fold} | {test_drugs} | {test_groups} | {test_wells} | {valid_drugs} | {valid_groups} | {valid_wells} | {train_drugs} | {train_groups} | {train_wells} | {test_plates} | {test_atc_moa_drugs} |".format(**row)
        )
    md += [
        "",
        f"`splits.parquet` is long format: {N_TREATED} treated wells × {N_FOLDS} folds = {N_TREATED * N_FOLDS:,} rows, primary key `(well, fold)`.",
        "Controls are never held out and are not split rows; each tested well uses its own plate's DMSO cells.",
        "PFOA (`SNGREZUHAYWORS-UHFFFAOYSA-N`) is one cluster and its 18 wells stay on one side; its test fold references both Plate5 and Plate9.",
    ]
    (out / "splits.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    parity = None
    if source_data is not None:
        parity = verify_source_parity(out, source_data, selected_index)

    manifest = build_manifest(
        out, source_info, var_names, groups, splits, true_sets, gates, sanity, parity
    )
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_readme(out, manifest)
    return manifest


SPLIT_ALGORITHM = (
    "458 connectivity-block clusters; order by n_wells desc then "
    "sha256(seed:cluster_id) asc; capacity-constrained greedy into the fold "
    "minimizing (treated wells, drugs, fold index) subject to <=92 drugs; "
    "valid_fold=(test_fold+4)%5"
)


def check(out: Path) -> None:
    manifest_path = out / "manifest.json"
    if not manifest_path.exists():
        raise BuildError(f"missing {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    for name, info in manifest["artifacts"].items():
        path = out / name
        if not path.exists():
            raise BuildError(f"missing artifact {name}")
        if sha256_file(path) != info["sha256"]:
            raise BuildError(f"hash mismatch: {name}")
        if path.stat().st_size != info["bytes"]:
            raise BuildError(f"size mismatch: {name}")
    print(f"OK: {len(manifest['artifacts'])} artifacts match manifest.json")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--gene-set", type=Path, default=DEFAULT_GENE_SET)
    parser.add_argument("--source-data", type=Path, default=DEFAULT_SOURCE_DATA)
    parser.add_argument(
        "--skip-source-check",
        action="store_true",
        help="do not compare the rebuilt arrays with the 18k archive",
    )
    parser.add_argument("--check", action="store_true", help="verify existing artifacts")
    args = parser.parse_args(argv)
    if args.check:
        check(args.out)
        return 0
    manifest = build(
        args.out,
        args.gene_set,
        None if args.skip_source_check else args.source_data,
    )
    print(json.dumps({k: manifest[k] for k in ("n_genes", "n_groups", "n_treated")}, indent=1))
    print(f"wrote {len(manifest['artifacts'])} artifacts to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
