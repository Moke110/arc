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
    read_h5ad_header,
)

_WELL_RE = re.compile(r"^([A-Z]+)0*(\d+)$")


# --------------------------------------------------------------------------- #
# Shared helpers
# --------------------------------------------------------------------------- #


def _matched_control_columns(obs: pd.DataFrame, qc_map: pd.DataFrame,
                             is_control: pd.Series, group_key: pd.Series,
                             group_has: pd.Series) -> tuple[pd.Series, pd.Series, str]:
    """Carry #2's matched-control decision (or derive it when absent).

    `matched_control_ok` answers "does this observation have an available
    matched control?" (#2's own decision, keyed by plate/vehicle for CM and by
    plate/cell line for Tahoe).  When the handoff does not carry it, the rule in
    the README is applied and the manifest records `derived`.  A control row is
    its own matched control, so its flags are `True` / empty key by
    construction.  A treated row flagged `True` by #2 must actually have a
    QC-pass control in its group — a contradiction fails the build.
    """
    controls = is_control.to_numpy(dtype=bool)
    source = "derived"
    if "matched_control_ok" in qc_map.columns:
        source = "handoff"
        ok = obs.index.to_series().map(qc_map["matched_control_ok"]).to_numpy(dtype=object)
        key = (obs.index.to_series().map(qc_map["matched_control_key"]).to_numpy(dtype=object)
               if "matched_control_key" in qc_map.columns else None)
    else:
        ok = group_has.to_numpy(dtype=object)
        key = None
    has = group_has.to_numpy(dtype=bool)
    ok_filled = np.where(controls, True,
                         np.where(pd.isna(ok), False, np.asarray(ok, dtype=object)))
    ok_filled = np.asarray(ok_filled, dtype=bool)
    if source == "handoff":
        contradicted = (~controls) & ok_filled & (~has)
        if contradicted.any():
            bad = obs.index[contradicted][:5].tolist()
            raise ContractError(
                "the QC handoff says these treated rows have a matched control, "
                f"but their group has no QC-pass control: {bad}")
    if key is None:
        key_values = group_key.astype(str).to_numpy()
    else:
        key_values = np.where(pd.isna(key), group_key.astype(str).to_numpy(),
                              np.asarray(key, dtype=object))
    key_values = np.where(controls, "", np.asarray(key_values, dtype=object)).astype(str)
    return (pd.Series(ok_filled, index=obs.index, name="matched_control_ok"),
            pd.Series(key_values, index=obs.index, name="matched_control_key"),
            source)


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
    pass_dmso = (out["is_control"].astype(bool)
                 & (out["vehicle"].astype(str) == "DMSO")
                 & out["qc_pass"].astype(bool))
    group_has = out["plate"].astype(str).map(
        out.loc[pass_dmso].groupby("plate", observed=True).size()).fillna(0) > 0
    matched_ok, matched_key, matched_source = _matched_control_columns(
        out, qc_map, is_control, out["plate"].astype(str), group_has)
    out["matched_control_ok"] = matched_ok.to_numpy(dtype=bool)
    out["matched_control_key"] = _as_categorical(matched_key.astype(str))
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
    decision_text = out["eligibility_decision"].astype(str)
    out["main_set"] = ((decision_text == "include") & out["qc_pass"]
                       & is_treated).astype(bool)
    # Ticket #11's expansion arm adds the QC-pass non-pharmacological compounds
    # (decision == exclude) to training; they are never scored.  `uncertain`
    # drugs are in neither set.
    out["expansion_set"] = ((decision_text == "exclude") & out["qc_pass"]
                            & is_treated).astype(bool)

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
    matched_ok, matched_key, matched_source = _matched_control_columns(
        out, qc_map, is_control, out["control_group"].astype(str), has_control)
    out["matched_control_ok"] = matched_ok.to_numpy(dtype=bool)
    out["matched_control_key"] = _as_categorical(matched_key.astype(str))
    # #2's matched-control decision is authoritative; `has_structure` still
    # excludes the two structure-less drugs from any use.
    usable = (~is_control) & out["has_structure"].astype(bool) & out["qc_pass"] \
        & out["matched_control_ok"].astype(bool)
    out["tahoe_usable"] = usable.astype(bool)

    # time_h stays NaN: no atlas table carries a time column (ADR 0031 D3).  The
    # 24 h claim is literature provenance only, never a written constant.
    if out["time_h"].notna().any():
        raise ContractError("Tahoe time_h is not NaN; refusing to overwrite a source value")
    out["time_source"] = TIME_SOURCE_TAHOE
    return out[TAHOE_OBS]


def tahoe_control_groups(obs: pd.DataFrame) -> pd.DataFrame:
    """Matched-control membership in the contracted **long** form.

    One row per ``(group, member observation)``: ``group_index``, ``group_key``
    (``plate|cell_line_id``), ``plate``, ``cell_line_id``, ``row`` (position in
    the artifact).  Only **QC-pass control** rows are members, so the mean of a
    group's ``log1p_cp10k`` rows is exactly the matched control.  A group that
    lost both controls to QC keeps its treated rows in the artifact but is not a
    group here (they are ``tahoe_usable=False`` and have no pair); controls are
    never averaged across cell lines and 67k x 18k control vectors are never
    materialised.
    """
    needed = {"control_group", "is_control", "plate", "cell_line_id", "qc_pass"}
    missing = needed - set(obs.columns)
    if missing:
        raise ContractError(f"tahoe_control_groups: obs missing {sorted(missing)}")
    control = obs["is_control"].astype(bool).to_numpy()
    qc_pass = obs["qc_pass"].astype(bool).to_numpy()
    positions = np.flatnonzero(control)
    group = obs["control_group"].astype(str).to_numpy()
    plate = obs["plate"].astype(str).to_numpy()
    line = obs["cell_line_id"].astype(str).to_numpy()

    rows = []
    group_index = 0
    for label in sorted(set(group[positions])):
        members = positions[group[positions] == label]
        passed = members[qc_pass[members]]
        key = f"{plate[members[0]]}|{line[members[0]]}"
        if len({(plate[m], line[m]) for m in members}) != 1:
            raise ContractError(
                f"control group {label!r} spans more than one (plate, cell_line)")
        for member in passed:
            rows.append({
                "group_index": group_index,
                "group_key": key,
                "plate": plate[member],
                "cell_line_id": line[member],
                "row": int(member),
                "n_controls": int(len(members)),
                "n_qc_pass": int(len(passed)),
                "qc_pass": True,
            })
        if passed.size:
            group_index += 1
    table = pd.DataFrame(rows)
    if len(table) == 0:
        raise ContractError("tahoe_control_groups: no QC-pass control group")
    for column, dtype in (("group_index", np.int32), ("row", np.int32),
                          ("n_controls", np.int16), ("n_qc_pass", np.int16)):
        table[column] = table[column].astype(dtype)
    # Every usable treated row must point at a group that has a member.
    usable = obs["tahoe_usable"].astype(bool).to_numpy()
    present = set(table["group_key"].astype(str))
    key_of_row = np.array([f"{p}|{l}" for p, l in zip(plate, line)], dtype=object)
    orphan = usable & ~np.isin(key_of_row, list(present))
    if orphan.any():
        bad = obs.index[orphan][:5].tolist()
        raise ContractError(
            "usable Tahoe rows reference a control group with no QC-pass control: "
            f"{bad}")
    return table


def tahoe_pairs(obs: pd.DataFrame, groups: pd.DataFrame) -> dict:
    """Treated row -> control ``group_index`` (ticket #4 reads `row`/`group_index`).

    Only treated rows whose ``(plate, cell_line)`` group has ≥ 1 QC-pass control
    are listed: the 204 kept observations without a matched control are reported
    in the manifest and never re-paired across cell lines.  `usable` and `split`
    are carried for audit; they do not change the mapping.
    """
    treated = ~obs["is_control"].astype(bool).to_numpy()
    plate = obs["plate"].astype(str).to_numpy()
    line = obs["cell_line_id"].astype(str).to_numpy()
    key_of_row = np.array([f"{p}|{l}" for p, l in zip(plate, line)], dtype=object)
    index_of_key = dict(zip(groups["group_key"].astype(str),
                            groups["group_index"].astype(np.int64)))
    treated_rows = np.flatnonzero(treated)
    rows = []
    group_index = []
    usable = []
    split = []
    skipped_no_control = 0
    for row in treated_rows:
        gid = index_of_key.get(key_of_row[row])
        if gid is None:
            skipped_no_control += 1
            continue
        rows.append(int(row))
        group_index.append(int(gid))
        usable.append(bool(obs["tahoe_usable"].to_numpy()[row]))
        split.append(str(obs["tahoe_split"].astype(str).to_numpy()[row]))
    if not rows:
        raise ContractError("tahoe_pairs: no treated row has a matched control")
    return {
        "row": np.asarray(rows, dtype=np.int32),
        "group_index": np.asarray(group_index, dtype=np.int32),
        "usable": np.asarray(usable, dtype=bool),
        "split": np.asarray(split, dtype=str),
        "skipped_no_control": skipped_no_control,
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


def _csr_row_slice(src_path: Path, n_rows: int, group: str = "X"):
    """Leading `n_rows` of a CSR group, read straight from the h5py datasets.

    Only the requested rows' `data`/`indices` are read, so `--tahoe-rows N`
    costs O(N x nnz-per-row) instead of the whole 8 GB matrix.
    """
    import h5py
    import scipy.sparse as sp
    with h5py.File(src_path, "r") as handle:
        node = handle[group]
        shape = tuple(int(v) for v in node.attrs["shape"])
        if n_rows > shape[0]:
            raise ContractError(
                f"{src_path}: row limit {n_rows} exceeds the {shape[0]} source rows")
        indptr = np.asarray(node["indptr"][: int(n_rows) + 1], dtype=np.int64)
        nnz = int(indptr[-1])
        data = np.asarray(node["data"][:nnz], dtype=node["data"].dtype)
        indices = np.asarray(node["indices"][:nnz], dtype=node["indices"].dtype)
    return sp.csr_matrix((data, indices, indptr), shape=(int(n_rows), shape[1]))


def write_tahoe_h5ad(src_path: Path, dst_path: Path, obs: pd.DataFrame,
                     var: pd.DataFrame, uns: dict, *,
                     row_limit: int | None = None) -> None:
    """Write `tahoe_transfer.h5ad`, streaming the CSR `X`/`layers` from the view.

    Full build: h5py copies the `X` and `layers` groups between files without
    materialising either matrix, then the new `obs`/`var`/`uns` groups are
    grafted in place of the source's.  `row_limit` (smoke only) reads a leading
    CSR slice instead, so the writer path is exercised without 16 GB.

    Neither path may call `anndata.read_h5ad(..., backed="r")`: backed mode
    materialises `layers` (8.1 GB here), which is outside the memory envelope.
    """
    import anndata as ad
    import h5py

    _assert_obs_columns(obs, TAHOE_OBS, "tahoe_transfer")
    if len(var) != N_GENES:
        raise ContractError(f"tahoe var has {len(var)} rows, expected {N_GENES}")

    if row_limit is not None:
        limit = int(row_limit)
        sub_obs = obs.iloc[:limit].copy()
        sub_obs.index = pd.Index(sub_obs.index.astype(str), name="well")
        adata = ad.AnnData(
            X=_csr_row_slice(src_path, limit, "X"),
            obs=sub_obs, var=var, uns=uns,
            layers={"log1p_cp10k": _csr_row_slice(src_path, limit,
                                                  "layers/log1p_cp10k")},
        )
        adata.write_h5ad(dst_path, compression="gzip")
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
    """Read back the written Tahoe artifact's header and assert the frozen schema."""
    header = read_h5ad_header(path)
    shape = tuple(header.shape)
    if expected_rows is not None and shape[0] != expected_rows:
        raise ContractError(f"{path}: {shape[0]} rows, expected {expected_rows}")
    if shape[1] != N_GENES:
        raise ContractError(f"{path}: {shape[1]} genes, expected {N_GENES}")
    if tuple(header.var.index.astype(str)) != tuple(header.var.index.astype(str).unique()) \
            or len(header.var) != N_GENES:
        raise ContractError(f"{path}: var does not hold {N_GENES} unique gene symbols")
    if list(header.var.columns) != VAR_COLUMNS:
        raise ContractError(
            f"{path}: var columns are {list(header.var.columns)}, expected {VAR_COLUMNS}")
    missing = [c for c in TAHOE_OBS if c not in header.obs.columns]
    if missing:
        raise ContractError(f"{path}: obs missing columns {missing}")
    if "control_mean" in header.uns:
        raise ContractError(f"{path}: inherited uns['control_mean'] must not be carried")
    if header.x_dtype != "int32":
        raise ContractError(f"{path}: X is {header.x_dtype}, expected int32")
    if header.layer_names != ["log1p_cp10k"]:
        raise ContractError(
            f"{path}: layers are {header.layer_names}, expected ['log1p_cp10k']")
    return {"rows": int(shape[0]), "cols": int(shape[1]),
            "nnz": int(header.x_nnz), "X_dtype": str(header.x_dtype),
            "layers": list(header.layer_names)}
