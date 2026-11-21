"""Cross the three tenv held doses within one frozen development population (#264)."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STEP / "output/api"))
from tsg_transfer import development, fit, tenv, transfer

PARENT = "drug_parent_inchikey"
IDENTITY = [PARENT, "dose_um"]
DOSES = tenv.DOSES
METRICS = ("pcc", "rmse", "direction_error", "amplitude_error")
SEEDS = (42, 17)
EPOCHS = 20
LEARNING_RATE = 1e-4


def complete_rows(obs, rows, *, balance):
    subset = obs.iloc[rows][IDENTITY].assign(row=rows)
    by_dose = subset.groupby(IDENTITY)["row"].apply(list).unstack()
    complete = by_dose.reindex(columns=DOSES).dropna()
    kept = []
    for profiles in complete.itertuples(index=False, name=None):
        n = min(map(len, profiles))
        for group in profiles:
            order = sorted(group, key=lambda row: hashlib.sha256(str(obs.index[row]).encode()).digest())
            kept.extend(order[:n] if balance else order)
    return np.asarray(sorted(kept), dtype=int)


def crossing_splits(data, train, validation):
    """Derive three disjoint fit, selection and assessment splits from outer development."""
    obs = data.obs
    if (not obs.index.is_unique or len(set(train)) != len(train)
            or len(set(validation)) != len(validation)):
        raise ValueError("Expected unique frozen development identities")
    validation = complete_rows(obs, validation, balance=False)
    balanced_train = complete_rows(obs, train, balance=True)
    base = development.split_validation(data, balanced_train, validation)
    base["selection"] = complete_rows(obs, base["selection"], balance=True)
    splits = {}
    for dose in DOSES:
        roles = {role: rows[obs.iloc[rows].dose_um.eq(dose).to_numpy()
                            if role == "assessment" else obs.iloc[rows].dose_um.ne(dose).to_numpy()]
                 for role, rows in base.items()}
        if not all(map(len, roles.values())):
            raise ValueError("Incomplete dose crossing")
        splits[dose] = roles
    if len({len(roles["train"]) for roles in splits.values()}) != 1:
        raise ValueError("Training profile budgets differ across doses")
    return splits


def assess_profiles(native, model, data, rows):
    prediction = fit.predict(native, model, data, rows, seed=fit.VALIDATION_SEED)
    truth = data.expression[rows]
    control = data.control[rows]
    pcc, rmse = development.delta_scores(truth, prediction, control)
    measured, predicted = truth - control, prediction - control
    direction = (np.sign(measured) != np.sign(predicted)).mean(axis=1)
    amplitude = np.abs(np.abs(measured) - np.abs(predicted)).mean(axis=1)
    return data.obs.iloc[rows][IDENTITY].assign(pcc=pcc, rmse=rmse, direction_error=direction,
                                                amplitude_error=amplitude)


def pair_scores(random, initialized):
    if (not random.index.is_unique or not initialized.index.is_unique
            or not random[IDENTITY].equals(initialized[IDENTITY])):
        raise ValueError("Initialization arms must have identical ordered assessment profiles")
    contrasts = random[IDENTITY].copy()
    for metric in METRICS:
        defined = np.isfinite(random[metric]) & np.isfinite(initialized[metric])
        contrasts[f"{metric}_random"] = random[metric].where(defined)
        contrasts[f"{metric}_transfer"] = initialized[metric].where(defined)
        contrasts[f"{metric}_change"] = (initialized[metric] - random[metric]).where(defined)
    groups = contrasts.groupby(IDENTITY)
    summary = groups.mean()
    summary["n_profiles"] = groups.size()
    for metric in METRICS:
        summary[f"{metric}_n"] = groups[f"{metric}_change"].count()
    return summary.reset_index()


def digest_profiles(profile_ids):
    ordered = json.dumps(profile_ids.tolist(), separators=(",", ":")).encode()
    return hashlib.sha256(ordered).hexdigest()


def design_report(data, splits, train, validation):
    obs = data.obs
    development_rows = np.asarray(sorted(set(train) | set(validation)))
    rows = []
    for held, split in splits.items():
        active = set(np.concatenate(list(split.values())).tolist())
        excluded = np.asarray(sorted(set(development_rows) - active))
        held_parents = set(obs.iloc[split["assessment"]][PARENT])
        omitted = obs.iloc[excluded]
        is_parent = omitted[PARENT].isin(held_parents).to_numpy()
        is_dose = omitted.dose_um.eq(held).to_numpy()
        if (is_parent & is_dose).any():
            raise ValueError("Matched assessment well was excluded")
        groups = dict(**split, excluded_held_parent=excluded[is_parent],
                      excluded_held_dose=excluded[is_dose],
                      excluded_incomplete_or_unbalanced=excluded[~(is_parent | is_dose)])
        for role, indices in groups.items():
            part = obs.iloc[indices]
            rows.append(dict(held_dose_um=held, role=role, n_profiles=len(part),
                             n_parents=part[PARENT].nunique(), profile_sha256=digest_profiles(part.index),
                             dose_counts={str(dose): int(part.dose_um.eq(dose).sum()) for dose in DOSES}))
    return rows


def bootstrap_summary(paired):
    rng = np.random.default_rng(264)
    outputs = []
    # Averaging seeds within parent gives a descriptive two-seed estimate.
    for metric in METRICS:
        scores = paired.pivot(index=PARENT, columns=["dose_um", "seed"], values=f"{metric}_change")
        scores = scores.dropna()
        if len(scores) < 2:
            raise ValueError(f"Fewer than two complete parent clusters for {metric}")
        values = scores.to_numpy()
        samples = rng.integers(len(scores), size=(4000, len(scores)))
        draws = values[samples].mean(axis=1)
        draws = pd.DataFrame(draws, columns=scores.columns)
        for dose in DOSES:
            dose_draws = draws[dose].mean(axis=1)
            dose_seeds = scores[dose].mean(axis=0)
            low, high = np.quantile(dose_draws, [.025, .975])
            outputs.append(dict(comparison="transfer-minus-random", dose_um=dose, metric=metric,
                                n_parents=len(scores), average=dose_seeds.mean(), ci_low=low, ci_high=high,
                                seed_42=dose_seeds[42], seed_17=dose_seeds[17]))
            if dose != 1.0:
                contrast = dose_draws - draws[1.0].mean(axis=1)
                seed_contrasts = dose_seeds - scores[1.0].mean(axis=0)
                low, high = np.quantile(contrast, [.025, .975])
                outputs.append(dict(comparison="gain-minus-1uM-gain", dose_um=dose, metric=metric,
                                    n_parents=len(scores), average=seed_contrasts.mean(),
                                    ci_low=low, ci_high=high, seed_42=seed_contrasts[42],
                                    seed_17=seed_contrasts[17]))
    return pd.DataFrame(outputs)


def run(args):
    data = fit.load_inputs(args.inputs, "tenv")
    # Cold-drug fold 0 withholds the same outer parents at every monitored dose.
    train = data.rows("cold-drug", 0, "train")
    validation = data.rows("cold-drug", 0, "validation")
    splits = crossing_splits(data, train, validation)
    report = {
        "schema": "stage2-opt-264-dose-crossing-v1",
        "tenv_manifest_sha256": fit.sha256(args.inputs / "manifest.json"),
        "outer_scenario": "cold-drug", "outer_fold": 0, "outer_split_seed": 226,
        "development_counts": design_report(data, splits, train, validation),
        "learning_rate": LEARNING_RATE, "epochs": EPOCHS, "seeds": list(SEEDS),
        "selection": "earliest maximum parent-equal Delta PCC",
        "assessment_sampling_seed": fit.VALIDATION_SEED, "bootstrap_seed": 264,
        "bootstrap_draws": 4000,
    }
    if args.design_only:
        fit.write_json(args.output, report)
        return
    source_fit, winner, source_doses = development.load_tahoe_winner(
        args.tahoe_inputs, args.experiment / "tahoe", args.source)
    report["winner_sha256"] = winner["checkpoint_sha256"]
    report["fits"] = []
    paired = []
    budgets = set()
    for dose, roles in splits.items():
        def initialize(model, transform):
            return transfer.initialize_from_tahoe(source_fit, model, winner["transform"], transform,
                                                   source_doses + DOSES, held_dose=dose)

        for seed in SEEDS:
            scores = []
            for arm, initializer in (("random", None), ("transfer", initialize)):
                print(f"Held {dose:g} uM, seed {seed}, arm {arm}", file=sys.stderr, flush=True)
                result = development.run_fit(data, args.source, roles, learning_rate=LEARNING_RATE,
                                             seed=seed, epochs=EPOCHS, initialize=initializer,
                                             assessment_evaluator=assess_profiles)
                chosen = result["selected"]["pcc"]
                updates = result["history"][-1]["updates"]
                budgets.add(updates)
                scores.append(chosen["scores"])
                report["fits"].append(dict(held_dose_um=dose, seed=seed, arm=arm,
                                           selected=chosen["selected"], assessment=chosen["assessment"],
                                           transform=result["transform"], updates=updates,
                                           history=result["history"]))
            paired.append(pair_scores(*scores).assign(seed=seed))
    if len(budgets) != 1:
        raise ValueError("Optimizer budgets differed despite equal profile counts")
    parent = pd.concat(paired, ignore_index=True)
    summary = bootstrap_summary(parent)
    args.output.mkdir()
    parent.to_csv(args.output / "parent_pairs.csv", index=False, float_format="%.12g")
    summary.to_csv(args.output / "summary.csv", index=False, float_format="%.12g")
    report["tables_sha256"] = {name: fit.sha256(args.output / name)
                               for name in ("parent_pairs.csv", "summary.csv")}
    fit.write_json(args.output / "run.json", report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--design-only", action="store_true")
    for name in ("tahoe-inputs", "experiment", "source"):
        parser.add_argument(f"--{name}", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("--output must be new")
    if not args.design_only and any(getattr(args, key) is None for key in
                                    ("tahoe_inputs", "experiment", "source")):
        parser.error("fitting requires --tahoe-inputs, --experiment and --source")
    run(args)


if __name__ == "__main__":
    main()
