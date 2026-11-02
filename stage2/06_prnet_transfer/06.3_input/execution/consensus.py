"""Per-fold CM consensus control `c*` (ADR 0032 decision 4, spec §控制与模型).

Deployment semantics frozen here: training sees matched same-plate DMSO
controls; validation, test and inference receive **one fixed `c*` per fold**;
the experimental truth Δ is computed separately at scoring time from matched
controls and never enters a model input.

`c*` is built only from control rows — never from a treated row of any drug —
and only from plates that carry at least one *training* drug of the fold.  A
plate whose fold drugs are all test/valid is a held-out environment and
contributes nothing.  Per plate the pool is averaged on the log scale, then
plates are weighted equally; spelling the weight per well
(`1 / (n_plates * n_wells_in_plate)`) makes the pool auditable and reproducible.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from contract import ContractError, N_GENES, N_FOLDS

CONSENSUS_SCHEMA = "06.3_input/consensus/1"

MEMBERSHIP_COLUMNS = [
    # `row` is the 0-based position in cm_transfer.h5ad obs/X/layers (the name
    # ticket #4's pipeline reads; keep it exact).
    "fold", "plate", "observation_id", "row", "weight", "plate_weight",
    "n_wells_plate", "qc_pass", "vehicle",
]


def allowed_plates(fold_table: pd.DataFrame, fold: int) -> list[str]:
    """Plates carrying at least one training drug of `fold`, sorted.

    Accepts both the long `(drug_key, fold, role, plate)` table and the
    collapsed export where a multi-plate drug's plates are ``;``-joined.
    """
    sub = fold_table.loc[(fold_table["fold"].astype(int) == int(fold))
                         & (fold_table["role"] == "train"), "plate"]
    plates: set[str] = set()
    for value in sub.astype(str):
        plates.update(part for part in value.split(";") if part)
    plates = sorted(plates)
    if not plates:
        raise ContractError(f"fold {fold} has no training-drug plate")
    return plates


def consensus_pool(obs: pd.DataFrame, fold_table: pd.DataFrame,
                   fold: int) -> pd.DataFrame:
    """The candidate control rows for one fold, with per-well weights."""
    for col in ("is_control", "vehicle", "plate", "qc_pass", "observation_id"):
        if col not in obs.columns:
            raise ContractError(f"consensus_pool: obs is missing {col!r}")
    plates = set(allowed_plates(fold_table, fold))
    is_control = obs["is_control"].astype(bool).to_numpy()
    is_dmso = obs["vehicle"].astype(str).to_numpy() == "DMSO"
    qc_pass = obs["qc_pass"].astype(bool).to_numpy()
    plate = obs["plate"].astype(str).to_numpy()
    in_scope = np.array([p in plates for p in plate])
    mask = is_control & is_dmso & qc_pass & in_scope
    rows = np.flatnonzero(mask)
    if rows.size == 0:
        raise ContractError(
            f"fold {fold}: no QC-pass DMSO control rows on the training plates "
            f"({sorted(plates)}) — cannot build c*")
    pool = pd.DataFrame({
        "fold": np.int16(fold),
        "plate": plate[rows],
        "observation_id": obs["observation_id"].astype(str).to_numpy()[rows],
        "row": rows.astype(np.int32),
    })
    counts = pool.groupby("plate", observed=True).size().rename("n_wells_plate")
    pool = pool.merge(counts, left_on="plate", right_index=True, how="left")
    n_plates = int(pool["plate"].nunique())
    pool["plate_weight"] = np.float32(1.0 / n_plates)
    pool["weight"] = (pool["plate_weight"].astype(np.float64)
                      / pool["n_wells_plate"].astype(np.float64)).astype(np.float32)
    pool["qc_pass"] = True
    pool["vehicle"] = "DMSO"
    return pool[MEMBERSHIP_COLUMNS].reset_index(drop=True)


def build_consensus(log1p, obs: pd.DataFrame, fold_table: pd.DataFrame, *,
                    n_folds: int = N_FOLDS, n_genes: int = N_GENES):
    """Compute `c*` per fold plus the auditable membership table.

    Returns `(c_star, membership, summary, plate_names, plate_weights)`.
    """
    if log1p.shape != (len(obs), n_genes):
        raise ContractError(
            f"consensus: log layer shape {log1p.shape} != {(len(obs), n_genes)}")
    folds = sorted(fold_table["fold"].astype(int).unique())
    if folds != list(range(n_folds)):
        raise ContractError(f"consensus: fold table covers {folds}, expected 0..{n_folds - 1}")

    c_star = np.zeros((n_folds, n_genes), dtype=np.float32)
    memberships = []
    summaries = []
    used_plates: set[str] = set()
    for fold in range(n_folds):
        pool = consensus_pool(obs, fold_table, fold)
        used_plates.update(pool["plate"].astype(str).unique())
        weighted = np.zeros(n_genes, dtype=np.float64)
        for plate, sub in pool.groupby("plate", observed=True, sort=True):
            idx = sub["row"].to_numpy()
            plate_mean = np.asarray(log1p[idx], dtype=np.float64).mean(axis=0)
            weighted += float(sub["plate_weight"].iloc[0]) * plate_mean
        if not np.isfinite(weighted).all():
            raise ContractError(f"fold {fold}: c* has non-finite values")
        c_star[fold] = weighted.astype(np.float32)
        memberships.append(pool)
        summaries.append({
            "fold": fold,
            "n_plates": int(pool["plate"].nunique()),
            "n_wells": int(len(pool)),
            "plates": ";".join(sorted(pool["plate"].astype(str).unique())),
            "plate_weight": float(1.0 / pool["plate"].nunique()),
        })

    membership = pd.concat(memberships, ignore_index=True)
    membership["fold"] = membership["fold"].astype(np.int16)
    membership["row"] = membership["row"].astype(np.int32)
    membership["weight"] = membership["weight"].astype(np.float32)
    membership["plate_weight"] = membership["plate_weight"].astype(np.float32)
    membership["n_wells_plate"] = membership["n_wells_plate"].astype(np.int32)
    summary = pd.DataFrame(summaries)

    plate_names = sorted(used_plates)
    plate_weights = np.zeros((n_folds, len(plate_names)), dtype=np.float32)
    for row in summary.itertuples(index=False):
        for plate in str(row.plates).split(";"):
            plate_weights[int(row.fold), plate_names.index(plate)] = float(row.plate_weight)

    return (c_star, membership, summary,
            np.array(plate_names, dtype=object), plate_weights)


def recompute_c_star(log1p, membership: pd.DataFrame, *,
                     n_folds: int = N_FOLDS, n_genes: int = N_GENES) -> np.ndarray:
    """Rebuild `c*` from the membership table alone (used by `--verify`)."""
    c_star = np.zeros((n_folds, n_genes), dtype=np.float32)
    if membership.empty:
        raise ContractError("consensus membership is empty")
    for fold in range(n_folds):
        sub = membership.loc[membership["fold"].astype(int) == fold]
        if sub.empty:
            raise ContractError(f"consensus membership has no fold {fold}")
        idx = sub["row"].to_numpy()
        weights = sub["weight"].to_numpy(dtype=np.float64)
        if abs(weights.sum() - 1.0) > 1e-5:
            raise ContractError(
                f"fold {fold}: consensus weights sum to {weights.sum():.8f}, not 1")
        weighted = (np.asarray(log1p[idx], dtype=np.float64) * weights[:, None]).sum(axis=0)
        if not np.isfinite(weighted).all():
            raise ContractError(f"fold {fold}: recomputed c* is non-finite")
        c_star[fold] = weighted.astype(np.float32)
    return c_star


def assert_consensus_provenance(membership: pd.DataFrame, obs: pd.DataFrame,
                                fold_table: pd.DataFrame, *,
                                n_folds: int = N_FOLDS) -> dict:
    """Fail closed unless every `c*` member is a permitted control row.

    The three leak assertions the ticket names:
      1. every member row is a control row (`is_control`), never a treated well;
      2. every member row is `vehicle == 'DMSO'` and QC-pass;
      3. every member row's plate carries at least one *training* drug of that
         fold (a held-out plate contributes nothing) — and the treated
         expression of test/valid drugs is therefore structurally absent from
         the pool.
    """
    for col in MEMBERSHIP_COLUMNS:
        if col not in membership.columns:
            raise ContractError(f"consensus membership is missing {col!r}")
    row_ix = membership["row"].to_numpy(dtype=np.int64)
    if row_ix.min() < 0 or row_ix.max() >= len(obs):
        raise ContractError("consensus membership row index outside obs")
    is_control = obs["is_control"].astype(bool).to_numpy()
    vehicle = obs["vehicle"].astype(str).to_numpy()
    qc_pass = obs["qc_pass"].astype(bool).to_numpy()
    plate = obs["plate"].astype(str).to_numpy()

    if not is_control[row_ix].all():
        bad = row_ix[~is_control[row_ix]][:5]
        raise ContractError(f"consensus pool contains TREATED rows: {bad.tolist()}")
    if not (vehicle[row_ix] == "DMSO").all():
        bad = row_ix[vehicle[row_ix] != "DMSO"][:5]
        raise ContractError(f"consensus pool contains non-DMSO rows: {bad.tolist()}")
    if not qc_pass[row_ix].all():
        bad = row_ix[~qc_pass[row_ix]][:5]
        raise ContractError(f"consensus pool contains QC-failed rows: {bad.tolist()}")
    if not (plate[row_ix] == membership["plate"].astype(str).to_numpy()).all():
        raise ContractError("consensus membership plate disagrees with obs.plate")

    facts = {"per_fold": [], "total_members": int(len(membership))}
    for fold in range(n_folds):
        sub = membership.loc[membership["fold"].astype(int) == fold]
        if sub.empty:
            raise ContractError(f"fold {fold} has no consensus membership")
        train_plates = set(allowed_plates(fold_table, fold))
        fold_plates = set(sub["plate"].astype(str))
        held_out = fold_plates - train_plates
        if held_out:
            raise ContractError(
                f"fold {fold}: consensus pool draws on held-out plate(s) "
                f"{sorted(held_out)} with no training drug")
        if abs(float(sub["weight"].sum()) - 1.0) > 1e-5:
            raise ContractError(
                f"fold {fold}: consensus weights sum to {sub['weight'].sum():.8f}")
        per_plate = sub.groupby("plate", observed=True)["weight"].sum()
        if not np.allclose(per_plate.to_numpy(), 1.0 / len(per_plate), atol=1e-6):
            raise ContractError(f"fold {fold}: plates are not equally weighted")
        facts["per_fold"].append({
            "fold": fold, "n_plates": int(len(fold_plates)),
            "n_wells": int(len(sub)),
            "plates": sorted(fold_plates),
        })
    return facts
