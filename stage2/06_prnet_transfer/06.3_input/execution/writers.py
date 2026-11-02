"""Observation assembly and h5ad writers for the frozen 06.3 input.

Two artifacts, one contract:

* `cm_transfer.h5ad` — dense `int32` raw counts + dense `float32`
  `log1p_cp10k`, 18,164-gene axis, the 24 inherited CM columns plus the identity,
  raw-design and eligibility/QC flags.
* `tahoe_transfer.h5ad` — the 04 pretrain view unchanged in `X`/`layers`
  (CSR, streamed by h5py so the 16 GB matrix is never materialised), with new
  identity/QC/split columns.

Neither artifact carries `uns['control_mean']`: the 02 cache mixes DMSO+Media
**raw counts** inside one group and the 04 cache is raw-count means against
log-scale targets (the ADR 0031 recorded deviation).  The frozen sources keep
their caches untouched; the frozen input ships the correctly-scaled caches
instead — per-fold `c*` (with row membership) for CM, and control-group
membership for Tahoe whose mean is taken on the log layer at load time.

Building requires the `arc` environment (anndata 0.12.16) because it writes the
`dataframe/0.2.0` obs/var encoding the pinned trainer (anndata 0.10.8) reads.
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from contract import (
    CM_OBS, CM_OBS_INHERITED, ContractError, N_GENES, TAHOE_OBS,
    TAHOE_OBS_INHERITED, TAHOE_VALID_PLATES, TIME_SOURCE_TAHOE, VAR_COLUMNS,
)

_WELL_RE = re.compile(r"^([A-Z]+)0*(\d+)$")


# --------------------------------------------------------------------------- #
# Shared helpers
# --------------------------------------------------------------------------- #

def var_frame(axis: pd.DataFrame, gene_measured: pd.Series,
              n_nonzero: np.ndarray) -> pd.DataFrame:
    """`var` with the corpus identity columns, measurement mask and nnz audit."""
    if len(axis) != N_GENES:
        raise ContractError(f"gene axis has {len(axis)} rows, expected {N_GENES}")
    if len(n_nonzero) != N_GENES:
        raise ContractError(f"n_nonzero has {len(n_nonzero)} entries, expected {N_GENES}")
    measured = gene_measured.reindex(axis.index)
    if measured.isna().any():
        raise ContractError("gene_measured does not cover the gene axis")
    var = pd.DataFrame({
        "hgnc_id": axis["hgnc_id"].astype(str),
        "entrez_id": axis["entrez_id"].astype(str),
        "ensembl_id": axis["ensembl_id"].astype(str),
        "match_type": axis["match_type"].astype(str),
        "gene_measured": measured.to_numpy(dtype=bool),
        "n_nonzero": np.asarray(n_nonzero, dtype=np.int32),
    }, index=axis.index)
    var.index.name = "symbol"
    return var[VAR_COLUMNS]


def _as_categorical(series: pd.Series) -> pd.Series:
    values = series.astype(str)
    if (values == "").all():
        return values.astype(object)
    return values.astype("category")


def _assert_obs_columns(obs: pd.DataFrame, expected: list[str], artifact: str) -> None:
    if list(obs.columns) != expected:
        raise ContractError(
            f"{artifact}: obs columns are not the frozen contract.\n"
            f"  expected: {expected}\n  got:      {list(obs.columns)}")
    if obs.index.name != "well":
        raise ContractError(f"{artifact}: obs index must be named 'well'")
    if not obs.index.is_unique:
        raise ContractError(f"{artifact}: obs index is not unique")
    if obs.isna().all(axis=1).any() or (obs.index.astype(str) == "").any():
        raise ContractError(f"{artifact}: empty observation id")


# --------------------------------------------------------------------------- #
# CM: raw-design re-join and obs assembly
# --------------------------------------------------------------------------- #

def parse_plate_well(well_ids) -> pd.DataFrame:
    """`temp_env::Plate8-B02` -> plate `Plate8`, plate_no 8, normalised well `B2`.

    The normalisation matches the raw GSE262419 `Well_ID` spelling (zero padding
    stripped); it is the join key the 02 loader itself used.
    """
    rows = []
    for well in well_ids:
        text = str(well)
        tail = text.split("-")[-1]
        match = _WELL_RE.match(tail)
        if not match or text.count("-") < 1:
            raise ContractError(f"cannot parse CM well id {text!r}")
        plate = text.split("::")[-1].split("-")[0]
        if not plate.startswith("Plate"):
            raise ContractError(f"cannot parse plate from CM well id {text!r}")
        rows.append({"well": text, "plate": plate,
                     "plate_no": int(plate.replace("Plate", "")),
                     "physical_well": f"{match.group(1)}{match.group(2)}"})
    return pd.DataFrame(rows).set_index("well")


def join_raw_design(obs: pd.DataFrame, hash_table: pd.DataFrame) -> pd.DataFrame:
    """Re-join `Replicate_Number`, `Solvent`, `Chemical_Index` from the raw table.

    The 02 loader drops all three (`cmdata/sources.py:480-555`), and
    `obs.replicate` is a generated index, so the raw design value must be
    recovered by `(plate, normalised well)`.  Verified 1:1 on all rows.
    """
    needed = {"Plate_ID", "Well_ID", "Replicate_Number", "Solvent",
              "Chemical_Index", "Chemical_name", "Chemical_Concentration_uM"}
    missing = needed - set(hash_table.columns)
    if missing:
        raise ContractError(f"GSE262419 hash table missing {sorted(missing)}")
    parsed = parse_plate_well(obs.index.to_numpy())
    joined = parsed.join(
        hash_table.set_index(["Plate_ID", "Well_ID"]),
        on=["plate_no", "physical_well"], how="left", validate="many_to_one")
    if joined[["Replicate_Number", "Solvent", "Chemical_Index"]].isna().any().any():
        orphans = joined.index[joined["Replicate_Number"].isna()][:5].tolist()
        raise ContractError(f"raw-design join left orphan CM wells: {orphans}")
    out = obs.copy()
    out["plate_no"] = joined["plate_no"].astype(np.int16)
    out["physical_well"] = joined["physical_well"].astype(str)
    out["replicate_raw"] = joined["Replicate_Number"].astype(np.int16)
    out["solvent"] = _as_categorical(joined["Solvent"].astype(str))
    out["chemical_index"] = _as_categorical(joined["Chemical_Index"].astype(str))
    return out


def assemble_cm_obs(obs: pd.DataFrame, hash_table: pd.DataFrame,
                    eligibility: pd.DataFrame, qc: pd.DataFrame) -> pd.DataFrame:
    """Attach the frozen flags to the inherited 24 CM columns."""
    missing = [c for c in CM_OBS_INHERITED if c not in obs.columns]
    if missing:
        raise ContractError(f"CM view is missing inherited columns {missing}")
    out = obs[CM_OBS_INHERITED].copy()
    out.index = pd.Index(obs.index.astype(str), name="well")
    out = join_raw_design(out, hash_table)
    out["observation_id"] = out.index.astype(str)

    drug_key = out["drug_key"].astype(str)
    decision = drug_key.map(eligibility.set_index("drug_key")["eligibility_decision"])
    elig_class = drug_key.map(eligibility.set_index("drug_key")["elig_class"])

    qc_map = qc.set_index("well")
    qc_pass = out.index.to_series().map(qc_map["qc_pass"])
    qc_reason = out.index.to_series().map(qc_map["qc_reason"])
    if qc_pass.isna().any():
        missing_ids = out.index[qc_pass.isna()][:5].tolist()
        raise ContractError(
            f"QC handoff does not cover {int(qc_pass.isna().sum())} CM rows, "
            f"e.g. {missing_ids}; the frozen input must not invent QC decisions")

    out["eligibility_decision"] = decision.astype("category")
    out["elig_class"] = elig_class.astype("category")
    out["qc_pass"] = qc_pass.astype(bool)
    out["qc_reason"] = _as_categorical(qc_reason.fillna("").astype(str))

    is_control = out["is_control"].astype(bool)
    is_treated = ~is_control
    resolvable = is_treated & (drug_key != "")
    if (is_treated & ~resolvable).any():
        raise ContractError(
            f"{int((is_treated & ~resolvable).sum())} treated CM rows have no drug_key")
    unresolved = is_treated & decision.isna()
    if unresolved.any():
        keys = sorted(drug_key[unresolved].unique())[:5]
        raise ContractError(
            f"treated CM rows reference drug_key(s) absent from the #1 eligibility "
            f"table: {keys}; #1 must cover every resolvable key")
    out["eligibility_decision"] = out["eligibility_decision"].cat.set_categories(
        ["include", "exclude", "uncertain"])
    main_set = (out["eligibility_decision"].astype(str) == "include") & out["qc_pass"]
    out["main_set"] = (main_set & is_treated).astype(bool)

    return out[CM_OBS]


def cm_keep_mask(obs: pd.DataFrame) -> pd.Series:
    """Rows the frozen CM artifact keeps: QC-pass treated (resolvable) + QC-pass controls."""
    is_control = obs["is_control"].astype(bool)
    qc_pass = obs["qc_pass"].astype(bool)
    resolvable = obs["drug_key"].astype(str) != ""
    return qc_pass & (is_control | (~is_control & resolvable))


# --------------------------------------------------------------------------- #
# Tahoe: identity / split / usability and obs assembly
# --------------------------------------------------------------------------- #

def parse_tahoe_well(well_ids) -> pd.DataFrame:
    """`tahoe::plate1:smp_1495:CVCL_0399` -> plate, physical well, cell line."""
    rows = []
    for well in well_ids:
        text = str(well)
        # tahoe::<plate>:<physical well>:<cell line>
        parts = text.split(":")
        if len(parts) != 5 or parts[0] != "tahoe" or parts[1] != "":
            raise ContractError(f"cannot parse Tahoe observation id {text!r}")
        _, _, plate, physical, cell_line = parts
        rows.append({"well": text, "plate_parsed": plate,
                     "physical_well": physical, "cell_line_parsed": cell_line})
    return pd.DataFrame(rows).set_index("well")


def tahoe_split(plate: pd.Series) -> pd.Series:
    """`train` | `valid` by the existing plate13/plate14 technical default."""
    plates = plate.astype(str)
    return pd.Series(np.where(plates.isin(TAHOE_VALID_PLATES), "valid", "train"),
                     index=plate.index, dtype=object)


def assemble_tahoe_obs(obs: pd.DataFrame, qc: pd.DataFrame) -> pd.DataFrame:
    """Attach identity, split, QC and usability flags to the 35 inherited columns."""
    missing = [c for c in TAHOE_OBS_INHERITED if c not in obs.columns]
    if missing:
        raise ContractError(f"Tahoe view is missing inherited columns {missing}")
    out = obs[TAHOE_OBS_INHERITED].copy()
    out.index = pd.Index(obs.index.astype(str), name="well")
    out["observation_id"] = out.index.astype(str)

    parsed = parse_tahoe_well(out.index.to_numpy())
    if not (parsed["plate_parsed"].to_numpy() == out["plate"].astype(str).to_numpy()).all():
        raise ContractError("Tahoe obs.well plate token disagrees with obs.plate")
    if not (parsed["cell_line_parsed"].to_numpy()
            == out["cell_line_id"].astype(str).to_numpy()).all():
        raise ContractError("Tahoe obs.well cell-line token disagrees with obs.cell_line_id")
    out["physical_well"] = parsed["physical_well"].astype(str)
    out["tahoe_split"] = _as_categorical(tahoe_split(out["plate"]))

    is_control = out["is_control"].astype(bool)
    dose = out["dose_um"].astype(float)
    line = out["cell_line_id"].astype(str)
    key = out["drug_key"].astype(str)
    name = out["compound_name"].astype(str)
    plate = out["plate"].astype(str)
    treated_group = np.where(key.to_numpy() != "", key.to_numpy(),
                            "unresolved:" + name.to_numpy())
    treated_key = (pd.Series(treated_group, index=out.index).astype(str)
                   + "|" + dose.map(lambda d: f"{d:g}").astype(str)
                   + "|" + line)
    control_key = "control|" + plate + "|" + line
    out["tahoe_valid_group"] = np.where(is_control.to_numpy(),
                                        control_key.to_numpy(),
                                        treated_key.to_numpy())

    qc_map = qc.set_index("well")
    qc_pass = out.index.to_series().map(qc_map["qc_pass"])
    qc_reason = out.index.to_series().map(qc_map["qc_reason"])
    if qc_pass.isna().any():
        missing_ids = out.index[qc_pass.isna()][:5].tolist()
        raise ContractError(
            f"QC handoff does not cover {int(qc_pass.isna().sum())} Tahoe rows, "
            f"e.g. {missing_ids}")
    out["qc_pass"] = qc_pass.astype(bool)
    out["qc_reason"] = _as_categorical(qc_reason.fillna("").astype(str))

    group_pass = (out.loc[out["is_control"].astype(bool)]
                  .groupby("control_group", observed=True)["qc_pass"].sum())
    has_control = out["control_group"].astype(str).map(group_pass).fillna(0) > 0
    usable = (~is_control) & out["has_structure"].astype(bool) & out["qc_pass"] & has_control
    out["tahoe_usable"] = usable.astype(bool)

    # time_h stays NaN: no atlas table carries a time column (ADR 0031 D3).  The
    # 24 h claim is literature provenance only, never a written constant.
    if out["time_h"].notna().any():
        raise ContractError("Tahoe time_h is not NaN; refusing to overwrite a source value")
    out["time_source"] = TIME_SOURCE_TAHOE
    return out[TAHOE_OBS]


def tahoe_control_groups(obs: pd.DataFrame) -> pd.DataFrame:
    """`control_group` -> the control row positions the matched mean is taken over.

    One row per `(plate, cell_line)` group; `control_rows` holds the row
    positions (as a list, parquet-friendly) so the pipeline can compute the mean
    of that group's `log1p_cp10k` rows directly.  Never materialise 67k x 18k
    control vectors: 700 groups x 2 rows is all the membership needed.
    """
    needed = {"control_group", "is_control", "plate", "cell_line_id", "qc_pass"}
    missing = needed - set(obs.columns)
    if missing:
        raise ContractError(f"tahoe_control_groups: obs missing {sorted(missing)}")
    control = obs["is_control"].astype(bool).to_numpy()
    positions = np.flatnonzero(control)
    group = obs["control_group"].astype(str).to_numpy()
    qc_pass = obs["qc_pass"].astype(bool).to_numpy()
    rows = []
    for label in sorted(set(group[positions])):
        idx = positions[group[positions] == label]
        first = int(idx[0])
        rows.append({
            "control_group": label,
            "plate": str(obs["plate"].iloc[first]),
            "cell_line_id": str(obs["cell_line_id"].iloc[first]),
            "control_rows": [int(i) for i in idx],
            "n_controls": int(len(idx)),
            "n_qc_pass": int(qc_pass[idx].sum()),
        })
    table = pd.DataFrame(rows)
    if len(table) == 0:
        raise ContractError("tahoe_control_groups: no control groups")
    if (table["n_qc_pass"] == 0).any():
        empty = table.loc[table["n_qc_pass"] == 0, "control_group"].tolist()[:5]
        raise ContractError(f"control group(s) lost every control to QC: {empty}")
    return table


def tahoe_pairs(obs: pd.DataFrame) -> dict:
    """Treated row -> its control group, plus the split/role arrays."""
    treated = ~obs["is_control"].astype(bool).to_numpy()
    groups = obs["control_group"].astype(str).to_numpy()
    return {
        "treated_row": np.flatnonzero(treated).astype(np.int32),
        "control_group": groups[treated],
        "usable": obs["tahoe_usable"].astype(bool).to_numpy()[treated],
        "split": obs["tahoe_split"].astype(str).to_numpy()[treated],
    }


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #

def write_cm_h5ad(path: Path, obs: pd.DataFrame, X: np.ndarray,
                  log1p: np.ndarray, var: pd.DataFrame, uns: dict) -> None:
    """Write `cm_transfer.h5ad` (dense int32 X + dense float32 log layer)."""
    import anndata as ad
    _assert_obs_columns(obs, CM_OBS, "cm_transfer")
    if X.dtype != np.int32:
        raise ContractError(f"cm X must be int32, got {X.dtype}")
    if log1p.dtype != np.float32:
        raise ContractError(f"cm log layer must be float32, got {log1p.dtype}")
    if X.shape != (len(obs), len(var)) or log1p.shape != X.shape:
        raise ContractError(
            f"cm shapes disagree: X {X.shape} layer {log1p.shape} "
            f"obs {len(obs)} var {len(var)}")
    if (X < 0).any():
        raise ContractError("cm X has negative counts")
    if not np.isfinite(log1p).all():
        raise ContractError("cm log layer has non-finite values")
    adata = ad.AnnData(X=X, obs=obs, var=var,
                       layers={"log1p_cp10k": log1p}, uns=uns)
    adata.write_h5ad(path, compression="gzip")


def write_tahoe_h5ad(src_path: Path, dst_path: Path, obs: pd.DataFrame,
                     var: pd.DataFrame, uns: dict, *,
                     row_limit: int | None = None) -> None:
    """Write `tahoe_transfer.h5ad`, streaming the CSR `X`/`layers` from the view.

    Full build: h5py copies the `X` and `layers` groups between files without
    materialising either matrix, then the new `obs`/`var`/`uns` groups are
    grafted in place of the source's.  `row_limit` (smoke only) materialises a
    leading slice instead, so the writer path is exercised without 16 GB.
    """
    import anndata as ad
    import h5py

    _assert_obs_columns(obs, TAHOE_OBS, "tahoe_transfer")
    if len(var) != N_GENES:
        raise ContractError(f"tahoe var has {len(var)} rows, expected {N_GENES}")

    if row_limit is not None:
        src = ad.read_h5ad(src_path, backed="r")
        try:
            sub = src[: int(row_limit)].to_memory()
            sub_obs = obs.iloc[: int(row_limit)]
            sub_obs.index = pd.Index(sub_obs.index.astype(str), name="well")
            adata = ad.AnnData(X=sub.X, obs=sub_obs, var=var,
                               layers={"log1p_cp10k": sub.layers["log1p_cp10k"]},
                               uns=uns)
            adata.write_h5ad(dst_path, compression="gzip")
        finally:
            src.file.close()
        return

    with tempfile.TemporaryDirectory(prefix="tahoe_shell_") as tmp:
        shell_path = Path(tmp) / "shell.h5ad"
        shell = ad.AnnData(obs=obs, var=var, uns=uns)
        shell.write_h5ad(shell_path, compression="gzip")

        with h5py.File(dst_path, "w") as out, h5py.File(src_path, "r") as src:
            for key, value in src.attrs.items():
                out.attrs[key] = value
            out.attrs["input_contract"] = "06.3_input/1.0.0"
            for group in ("X", "layers", "obsm", "obsp", "varm", "varp"):
                if group in src:
                    src.copy(group, out)
            with h5py.File(shell_path, "r") as sh:
                for group in ("obs", "var", "uns"):
                    if group in out:
                        del out[group]
                    sh.copy(group, out)


def verify_tahoe_h5ad(path: Path, *, expected_rows: int | None = None) -> dict:
    """Read back the written Tahoe artifact and assert the frozen schema."""
    import anndata as ad
    adata = ad.read_h5ad(path, backed="r")
    try:
        shape = tuple(adata.shape)
        if expected_rows is not None and shape[0] != expected_rows:
            raise ContractError(f"{path}: {shape[0]} rows, expected {expected_rows}")
        if shape[1] != N_GENES:
            raise ContractError(f"{path}: {shape[1]} genes, expected {N_GENES}")
        if tuple(adata.var.index) != tuple(adata.var_names):
            raise ContractError(f"{path}: var index is not the gene axis")
        missing = [c for c in TAHOE_OBS if c not in adata.obs.columns]
        if missing:
            raise ContractError(f"{path}: obs missing columns {missing}")
        if "control_mean" in adata.uns:
            raise ContractError(f"{path}: inherited uns['control_mean'] must not be carried")
        X = adata.X
        if getattr(X, "dtype", None) != np.int32:
            raise ContractError(f"{path}: X is {getattr(X, 'dtype', None)}, expected int32")
        nnz = getattr(X, "nnz", None)
        if nnz is None:                       # CSRDataset exposes only shape/dtype
            import h5py
            with h5py.File(path, "r") as fh:
                nnz = int(fh["X/indptr"][-1])
        return {"rows": int(shape[0]), "cols": int(shape[1]),
                "nnz": int(nnz), "X_dtype": str(X.dtype),
                "layers": sorted(map(str, adata.layers.keys()))}
    finally:
        adata.file.close()
