"""Five drug-key folds for the frozen 06.3 input (ADR 0032 decision 2).

The split unit is the **drug key**: every dose and every replicate of a drug
moves together, every main-set drug is test in exactly one fold and valid in
exactly one fold, and salt/alias pairs cannot straddle a boundary because
`drug_key` already merges them while a supplementary connectivity-block union
(`inchikey[:14]`) keeps structures that share a skeleton in one cluster.

Pure functions on plain frames — no h5ad, no I/O — so the algorithm is testable
without the 16 GB artifacts.  `build_input.py` supplies `drugs` from the frozen
CM view.
"""

from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd

from contract import ContractError, FOLD_ROLES, N_FOLDS, FOLD_SEED

# Columns a caller must supply per drug.
DRUG_COLUMNS = ("drug_key", "inchikey")


def _digest(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def cluster_drugs(drugs: pd.DataFrame, *, n_folds: int = N_FOLDS,
                  seed: int = FOLD_SEED) -> pd.DataFrame:
    """Assign each main-set drug to a cluster and to a test/valid fold.

    Parameters
    ----------
    drugs
        One row per drug with `drug_key`, `inchikey` and `n_wells` (main-set
        wells).  `inchikey` may be empty; then the drug is its own cluster.
    n_folds, seed
        Frozen contract values (`contract.N_FOLDS`, `contract.FOLD_SEED`).

    Returns
    -------
    DataFrame with `drug_key`, `cluster_id`, `test_fold`, `valid_fold`,
    `n_wells`.
    """
    missing = [c for c in DRUG_COLUMNS if c not in drugs.columns]
    if missing:
        raise ContractError(f"cluster_drugs: missing column(s) {missing}")
    if "n_wells" not in drugs.columns:
        raise ContractError("cluster_drugs: missing column 'n_wells'")
    table = drugs[["drug_key", "inchikey", "n_wells"]].copy()
    table["drug_key"] = table["drug_key"].astype(str)
    table["inchikey"] = table["inchikey"].fillna("").astype(str)
    if table["drug_key"].duplicated().any():
        raise ContractError("cluster_drugs: duplicated drug_key rows")
    if len(table) == 0:
        raise ContractError("cluster_drugs: no drugs")
    if len(table) < n_folds:
        raise ContractError(
            f"cluster_drugs: {len(table)} drugs cannot fill {n_folds} folds")

    # Union-find over drug_key and the 14-character InChIKey connectivity block.
    parent: dict[str, str] = {}

    def find(node: str) -> str:
        parent.setdefault(node, node)
        root = node
        while parent[root] != root:
            root = parent[root]
        while parent[node] != root:
            parent[node], node = root, parent[node]
        return root

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for drug_key, inchikey in zip(table["drug_key"], table["inchikey"]):
        key_node = "K:" + drug_key
        parent.setdefault(key_node, key_node)
        block = inchikey[:14]
        if block:
            union(key_node, "B:" + block)

    members: dict[str, list[str]] = {}
    for drug_key in sorted(table["drug_key"]):
        members.setdefault(find("K:" + drug_key), []).append(drug_key)

    cluster_of: dict[str, str] = {}
    for root, keys in members.items():
        cluster_id = "|".join(sorted(keys))
        for key in keys:
            cluster_of[key] = cluster_id

    table["cluster_id"] = table["drug_key"].map(cluster_of)
    cluster = (table.groupby("cluster_id", observed=True)
               .agg(n_drugs=("drug_key", "nunique"), n_wells=("n_wells", "sum"))
               .reset_index())
    cluster["order_key"] = cluster["cluster_id"].map(
        lambda cid: _digest(f"{seed}:{cid}"))
    # Multi-drug clusters first (they are the ones a capacity pass can strand),
    # then by wells and by a seeded hash order — fully deterministic.
    cluster = cluster.sort_values(["n_drugs", "n_wells", "order_key"],
                                 ascending=[False, False, True],
                                 kind="mergesort").reset_index(drop=True)

    # Capacity keeps the drug split ~6:2:2 (test:valid:train) by construction.
    per_fold = -(-len(table) // n_folds)          # ceil(n_drugs / n_folds)
    max_cluster = int(cluster["n_drugs"].max())
    if max_cluster > per_fold:
        raise ContractError(
            f"one identity cluster holds {max_cluster} of {len(table)} drugs "
            f"(capacity {per_fold}/fold): a five-fold drug split cannot keep the "
            "cluster together; the build fails closed instead of silently "
            "degrading the split")

    def greedy(capacity: int | None) -> dict[str, int] | None:
        load_drugs = [0] * n_folds
        load_wells = [0] * n_folds
        assignment: dict[str, int] = {}
        for row in cluster.itertuples(index=False):
            candidates = [k for k in range(n_folds)
                          if capacity is None
                          or load_drugs[k] + int(row.n_drugs) <= capacity]
            if not candidates:
                return None
            k = min(candidates, key=lambda f: (load_drugs[f], load_wells[f], f))
            load_drugs[k] += int(row.n_drugs)
            load_wells[k] += int(row.n_wells)
            assignment[str(row.cluster_id)] = k
        return assignment

    fold_of = greedy(per_fold)
    if fold_of is None:                           # pragma: no cover - defended above
        raise ContractError(
            "the capacity-constrained greedy pass could not place every cluster "
            f"(capacity {per_fold}/fold, largest cluster {max_cluster} drugs)")
    realised = [0] * n_folds
    for cluster_id, fold in fold_of.items():
        realised[fold] += int(
            cluster.loc[cluster["cluster_id"] == cluster_id, "n_drugs"].iloc[0])
    if max(realised) > per_fold:                  # pragma: no cover - greedy guarantees
        raise ContractError(f"fold sizes {realised} exceed {per_fold} drugs/fold")

    table["test_fold"] = table["cluster_id"].map(fold_of).astype(np.int16)
    # Every drug is valid in exactly one fold; offset -1 spreads valid sets.
    table["valid_fold"] = ((table["test_fold"].astype(int) + n_folds - 1) % n_folds
                           ).astype(np.int16)
    return table.sort_values("drug_key", kind="mergesort").reset_index(drop=True)


def expand_fold_table(assign: pd.DataFrame, drug_plates: pd.DataFrame, *,
                      n_folds: int = N_FOLDS) -> pd.DataFrame:
    """One row per (drug_key, fold, role, plate) — the frozen `folds.parquet`.

    `drug_plates` has `drug_key`, `plate`, `n_wells` and `n_doses`.
    """
    required = {"drug_key", "plate", "n_wells", "n_doses"}
    missing = required - set(drug_plates.columns)
    if missing:
        raise ContractError(f"expand_fold_table: missing column(s) {sorted(missing)}")
    base = assign[["drug_key", "cluster_id", "test_fold", "valid_fold"]]
    merged = drug_plates.merge(base, on="drug_key", how="left", validate="many_to_one")
    if merged["test_fold"].isna().any():
        unknown = merged.loc[merged["test_fold"].isna(), "drug_key"].unique()[:5]
        raise ContractError(f"expand_fold_table: drugs without assignment: {list(unknown)}")

    frames = []
    for fold in range(n_folds):
        frame = merged.copy()
        frame["fold"] = np.int16(fold)
        frame["role"] = np.where(frame["test_fold"] == fold, "test",
                                 np.where(frame["valid_fold"] == fold, "valid",
                                          "train"))
        frames.append(frame)
    table = pd.concat(frames, ignore_index=True)
    table = table[[
        "drug_key", "cluster_id", "fold", "role", "plate", "n_wells", "n_doses",
        "test_fold", "valid_fold",
    ]].sort_values(["fold", "role", "drug_key", "plate"], kind="mergesort"
                   ).reset_index(drop=True)
    table["fold"] = table["fold"].astype(np.int16)
    table["test_fold"] = table["test_fold"].astype(np.int16)
    table["valid_fold"] = table["valid_fold"].astype(np.int16)
    table["n_wells"] = table["n_wells"].astype(np.int32)
    table["n_doses"] = table["n_doses"].astype(np.int32)
    for col in ("drug_key", "cluster_id", "role", "plate"):
        table[col] = table[col].astype(str)
    return table


def collapse_fold_table(table: pd.DataFrame) -> pd.DataFrame:
    """The exported `folds.parquet`: exactly one row per ``(drug_key, fold)``.

    Consumers key the split by ``drug_key`` only, so the per-plate expansion is
    collapsed: a drug's plates are ``;``-joined (one key spans two plates) and
    its wells/doses are summed/maxed.  The per-plate long form is recomputed on
    demand from the CM artifact, so no information is lost — the collapse keeps
    the table's primary key unambiguous.
    """
    required = {"drug_key", "cluster_id", "fold", "role", "test_fold",
                "valid_fold", "plate", "n_wells", "n_doses"}
    missing = required - set(table.columns)
    if missing:
        raise ContractError(f"collapse_fold_table: missing column(s) {sorted(missing)}")
    grouped = (table.groupby(["drug_key", "fold"], observed=True)
               .agg(cluster_id=("cluster_id", "first"),
                    role=("role", "first"),
                    plates=("plate", lambda values: ";".join(sorted(set(values)))),
                    n_wells=("n_wells", "sum"),
                    n_doses=("n_doses", "max"),
                    test_fold=("test_fold", "first"),
                    valid_fold=("valid_fold", "first"))
               .reset_index())
    grouped["n_plates"] = grouped["plates"].str.count(";") + 1
    grouped = grouped.rename(columns={"plates": "plate"})
    out = grouped[["drug_key", "cluster_id", "fold", "role", "plate", "n_wells",
                   "n_doses", "n_plates", "test_fold", "valid_fold"]]
    out = out.sort_values(["fold", "role", "drug_key"], kind="mergesort"
                         ).reset_index(drop=True)
    out["fold"] = out["fold"].astype(np.int16)
    out["test_fold"] = out["test_fold"].astype(np.int16)
    out["valid_fold"] = out["valid_fold"].astype(np.int16)
    out["n_wells"] = out["n_wells"].astype(np.int32)
    out["n_doses"] = out["n_doses"].astype(np.int32)
    out["n_plates"] = out["n_plates"].astype(np.int8)
    for column in ("drug_key", "cluster_id", "role", "plate"):
        out[column] = out[column].astype(str)
    if out.duplicated(["drug_key", "fold"]).any():
        raise ContractError("folds.parquet would still duplicate a (drug_key, fold)")
    return out


def assert_fold_invariants(table: pd.DataFrame, assign: pd.DataFrame, *,
                           expected_drugs: set[str] | None = None,
                           n_folds: int = N_FOLDS) -> dict:
    """Fail closed on any violation of the frozen fold contract.

    Returns a small dict of measured facts for the manifest.
    """
    need = {"drug_key", "fold", "role", "plate", "cluster_id"}
    missing = need - set(table.columns)
    if missing:
        raise ContractError(f"fold table missing column(s) {sorted(missing)}")
    bad_roles = set(table["role"]) - set(FOLD_ROLES)
    if bad_roles:
        raise ContractError(f"fold table has roles outside {FOLD_ROLES}: {sorted(bad_roles)}")
    if set(table["fold"].astype(int)) != set(range(n_folds)):
        raise ContractError("fold table does not cover exactly folds 0..n-1")

    # One role per (drug, fold); roles repeat only across plates.
    per_drug_fold = table.groupby(["drug_key", "fold"], observed=True)["role"].nunique()
    if (per_drug_fold != 1).any():
        bad = per_drug_fold[per_drug_fold != 1].index.tolist()[:5]
        raise ContractError(f"a drug has more than one role inside a fold: {bad}")

    roles = table.drop_duplicates(["drug_key", "fold"])[["drug_key", "fold", "role"]]
    counts = roles.groupby(["drug_key", "role"], observed=True).size().unstack(fill_value=0)
    for role in FOLD_ROLES:
        if role not in counts.columns:
            counts[role] = 0
    for drug, row in counts.iterrows():
        if int(row["test"]) != 1 or int(row["valid"]) != 1 or int(row["train"]) != n_folds - 2:
            raise ContractError(
                f"drug {drug} has test={int(row['test'])}, valid={int(row['valid'])}, "
                f"train={int(row['train'])}; expected 1/1/{n_folds - 2}")
    if len(counts) != len(assign):
        raise ContractError("fold table drugs != fold assignment drugs")

    # Every drug's (test_fold, valid_fold) must equal the assignment, and a
    # cluster may never be split across test/valid folds.
    cluster_folds = table.groupby("cluster_id", observed=True)[["test_fold", "valid_fold"]].nunique()
    if (cluster_folds > 1).any().any():
        bad = cluster_folds[(cluster_folds > 1).any(axis=1)].index.tolist()[:5]
        raise ContractError(f"cluster split across folds: {bad}")

    if expected_drugs is not None:
        got = set(table["drug_key"].astype(str))
        missing_drugs = sorted(expected_drugs - got)
        extra = sorted(got - expected_drugs)
        if missing_drugs or extra:
            raise ContractError(
                f"fold table drug set mismatch: missing={missing_drugs[:5]} extra={extra[:5]}")

    sizes = fold_sizes(table)
    return {
        "n_drugs": int(table["drug_key"].nunique()),
        "n_rows": int(len(table)),
        "per_fold": sizes.to_dict(orient="records"),
    }


def fold_sizes(table: pd.DataFrame) -> pd.DataFrame:
    """Per (fold, role): drugs, wells, plates, doses — the tracked `fold_sizes.csv`."""
    grouped = (table.groupby(["fold", "role"], observed=True)
               .agg(n_drugs=("drug_key", "nunique"),
                    n_wells=("n_wells", "sum"),
                    n_plates=("plate", "nunique"),
                    n_doses=("n_doses", "sum"))
               .reset_index())
    grouped["fold"] = grouped["fold"].astype(int)
    grouped["role"] = grouped["role"].astype(str)
    full = (grouped.pivot(index="fold", columns="role", values="n_drugs")
            .reindex(columns=list(FOLD_ROLES)).fillna(0).astype(int))
    wells = (grouped.pivot(index="fold", columns="role", values="n_wells")
             .reindex(columns=list(FOLD_ROLES)).fillna(0).astype(int))
    out = pd.DataFrame({
        "fold": full.index,
        "test_drugs": full["test"].values,
        "valid_drugs": full["valid"].values,
        "train_drugs": full["train"].values,
        "test_wells": wells["test"].values,
        "valid_wells": wells["valid"].values,
        "train_wells": wells["train"].values,
    })
    plates = (table.groupby(["fold", "role"], observed=True)["plate"].nunique()
              .unstack(fill_value=0).reindex(columns=list(FOLD_ROLES)).fillna(0).astype(int))
    out["test_plates"] = plates["test"].values
    out["valid_plates"] = plates["valid"].values
    out["train_plates"] = plates["train"].values
    return out.sort_values("fold").reset_index(drop=True)


def _split_plates(values) -> list[str]:
    out: set[str] = set()
    for value in values:
        out.update(part for part in str(value).split(";") if part)
    return sorted(out)


def fold_plates(table: pd.DataFrame, fold: int, role: str) -> list[str]:
    """Sorted plates carrying a drug of `role` in `fold` (the consensus scope)."""
    sub = table.loc[(table["fold"].astype(int) == int(fold))
                    & (table["role"] == role), "plate"]
    return _split_plates(sub)
