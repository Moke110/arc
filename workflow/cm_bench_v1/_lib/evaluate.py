"""
Adapter evaluation over the two generalization axes.

Cold-drug (primary): for each of the 5 folds, finetune on the fold's train
drugs' wells and predict the held-out drugs' wells; aggregate well predictions
to (dataset, drug, dose) groups by mean and score direction-aware AUPRC against
the true sets. Each eligible drug is held out exactly once.

Cold-dose (secondary): leave-dose-level-out on gse244740 only (interpolation =
one interior level, extrapolation = the top level).

Aggregation views (reported side by side):
  macro           per-dataset mean of group AUPRC, then equally-weighted mean
                  across datasets (the headline number)
  micro           one AUPRC over all test groups' predicted/label vectors pooled
  per_drug_equal  mean over drugs of that drug's mean group AUPRC
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from adapter import Cases  # noqa: E402
from metric import direction_aware_auprc, prevalence  # noqa: E402


def _label_map(groups, labels):
    return {k: labels[i] for i, k in enumerate(groups["group_key"])}


def _reliable(groups):
    return set(groups.loc[groups["reliable"], "group_key"])


def _score_group(rows, pooled, sub, idx, pred, key, lab, axis, tag):
    idx = np.asarray(idx)
    p = pred[idx].mean(axis=0)
    a = direction_aware_auprc(p, lab)
    first = sub.iloc[idx[0]]
    rows.append({"axis": axis, "split": tag, "group_key": key,
                 "dataset": first["dataset"], "drug": first["drug"],
                 "dose_level": int(first["dose_level"]), "n_wells": int(len(idx)),
                 "auprc": a, "prevalence": prevalence(lab)})
    if np.isfinite(a):
        pooled.append((p, lab))


def evaluate_cold_drug(adapter, meta, baseline, delta, groups, labels, fold_df):
    lab_of = _label_map(groups, labels)
    reliable = _reliable(groups)
    rows, pooled = [], []
    folds = sorted(fold_df.loc[fold_df["eligible"], "fold"].unique())
    for f in folds:
        held = set(fold_df.loc[(fold_df["fold"] == f) & fold_df["eligible"], "drug"])
        test_keys = set(groups.loc[groups["drug"].isin(held)
                                   & groups["group_key"].isin(reliable), "group_key"])
        if not test_keys:
            continue
        test_mask = meta["group_key"].isin(test_keys).to_numpy()
        train_mask = ~meta["drug"].isin(held).to_numpy()
        train = Cases.from_frame(meta.loc[train_mask], baseline[train_mask],
                                 delta[train_mask])
        test = Cases.from_frame(meta.loc[test_mask], baseline[test_mask])
        pred = adapter.predict(train, test)
        sub = meta.loc[test_mask]
        for key, idx in sub.groupby("group_key").indices.items():
            _score_group(rows, pooled, sub, idx, pred, key, lab_of[key],
                         "cold-drug", f"fold{f}")
    return pd.DataFrame(rows), pooled


def evaluate_cold_dose(adapter, meta, baseline, delta, groups, labels, dose_df):
    if not getattr(adapter, "supports_dose", False):
        return pd.DataFrame(), None  # N/A (no dose input)
    lab_of = _label_map(groups, labels)
    # cold-dose true sets rest on ~3 wells per (drug, dose): keep every group
    # with >= 2 wells (a t-test needs 2) rather than the 4-well cold-drug bar,
    # and report the weakness with the numbers (spec: cold-dose caveat).
    scoreable = set(groups.loc[groups["n_wells"] >= 2, "group_key"])
    g = meta[meta["dataset"] == "gse244740"].reset_index(drop=True)
    g_base = baseline[meta["dataset"].to_numpy() == "gse244740"]
    g_delta = delta[meta["dataset"].to_numpy() == "gse244740"]
    rows, pooled = [], []
    for split in ("interpolation", "extrapolation"):
        sel = dose_df[(dose_df["split"] == split) & dose_df["held_out"]]
        if sel.empty:
            continue
        level = int(sel["dose_level"].iloc[0])
        test_keys = set(groups.loc[groups["group_key"].isin(scoreable)
                                   & (groups["dataset"] == "gse244740")
                                   & (groups["dose_level"] == level), "group_key"])
        if not test_keys:
            continue
        test_mask = g["group_key"].isin(test_keys).to_numpy()
        train_mask = (g["dose_level"] != level).to_numpy()
        train = Cases.from_frame(g.loc[train_mask], g_base[train_mask],
                                 g_delta[train_mask])
        test = Cases.from_frame(g.loc[test_mask], g_base[test_mask])
        pred = adapter.predict(train, test)
        sub = g.loc[test_mask]
        for key, idx in sub.groupby("group_key").indices.items():
            _score_group(rows, pooled, sub, idx, pred, key, lab_of[key],
                         "cold-dose", split)
    return pd.DataFrame(rows), pooled


def _pooled_micro(pooled):
    if not pooled:
        return {"macro": None, "micro": None, "per_drug_equal": None}
    preds = np.concatenate([p for p, _ in pooled])
    labs = np.concatenate([l for _, l in pooled])
    return {"micro": float(direction_aware_auprc(preds, labs))}


def summarize(df, pooled=None):
    """Aggregate a per-group result frame into macro / micro / per-drug-equal."""
    if df.empty:
        return {"macro": None, "micro": None, "per_drug_equal": None,
                "n_groups": 0, "n_datasets": 0, "n_drugs": 0}
    valid = df.dropna(subset=["auprc"])
    if valid.empty:
        return {"macro": None, "micro": None, "per_drug_equal": None,
                "n_groups": 0, "n_datasets": 0, "n_drugs": 0}
    macro = float(valid.groupby("dataset")["auprc"].mean().mean())
    per_drug = float(valid.groupby("drug")["auprc"].mean().mean())
    micro = _pooled_micro(pooled)["micro"] if pooled else macro
    return {"macro": macro, "micro": micro, "per_drug_equal": per_drug,
            "n_groups": int(len(valid)), "n_datasets": int(valid["dataset"].nunique()),
            "n_drugs": int(valid["drug"].nunique())}
