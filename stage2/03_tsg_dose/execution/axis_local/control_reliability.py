"""Development-only repeat and residual audits for the accepted Step 03 inputs."""
import argparse
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from baselines import SOURCE_POPULATIONS, STEP, frozen_pair, mean_cell, sha256
from tsg_round1.inputs import Inputs

SCENARIOS = ("cold-drug", "cold-drug-dose")
IDENTITY = ["profile_id", "drug_parent_inchikey", "dose_um"]
CONDITION = ["drug_parent_inchikey", "dose_um", "cell_line"]
SEED = 262
DRAWS = 2000


def pcc(x, y):
    x = x - x.mean()
    y = y - y.mean()
    norm = np.linalg.norm(x) * np.linalg.norm(y)
    return float(np.dot(x, y) / norm) if norm else np.nan


def paired_control_scores(x, y, pools):
    """Compare both assignments of equally sized, disjoint same-batch pools."""
    a, b = pools
    scores = [pcc(x - a, y - a), pcc(x - b, y - b),
              pcc(x - a, y - b), pcc(x - b, y - a)]
    if not np.isfinite(scores).all():
        return np.nan, np.nan
    return np.mean(scores[:2]), np.mean(scores[2:])


def control_pools(data, development):
    pools = {}
    for row in development:
        indices = data.controls[row]
        treated = data.all_obs.iloc[row]
        group = treated.matched_control_group
        if group in pools:
            continue
        vehicles = data.all_obs.iloc[indices].sort_values(["source_item", "source_row"])
        if (vehicles.empty or not vehicles.role.eq("vehicle").all()
                or not vehicles.plate.eq(treated.plate).all()
                or not vehicles.cell_line.eq(treated.cell_line).all()
                or not vehicles.source_item.eq(treated.source_item).all()
                or vehicles[["source_item", "source_row"]].duplicated().any()):
            raise ValueError(f"control identity mismatch for {group}")
        if pd.notna(treated.vehicle) and not vehicles.vehicle.eq(treated.vehicle).all():
            raise ValueError(f"vehicle mismatch for {group}")
        ordered = vehicles.index
        vehicle_rows = data.all_obs.index.get_indexer(ordered)
        n = len(vehicle_rows) // 2
        pools[group] = (set(vehicle_rows), None if n == 0 else np.stack([
            data.expression[vehicle_rows[:2*n:2]].mean(axis=0, dtype=np.float64),
            data.expression[vehicle_rows[1:2*n:2]].mean(axis=0, dtype=np.float64),
        ]))
    return pools


def repeat_scores(data, rows, pools):
    """One row per condition, averaging all eligible treated pairs."""
    records = []
    for (parent, dose, line), condition in data.obs.iloc[rows].groupby(CONDITION):
        group_rows = data.obs.index.get_indexer(condition.index)
        same, paired, cross = [], [], []
        reasons = Counter()
        for a, b in itertools.combinations(group_rows, 2):
            key_a = data.obs.iloc[a].matched_control_group
            key_b = data.obs.iloc[b].matched_control_group
            vehicle_a, halves = pools[key_a]
            vehicle_b, _ = pools[key_b]
            x, y = data.expression[a], data.expression[b]
            original = pcc(x - data.control[a], y - data.control[b])
            if key_a == key_b:
                same.append(original)
                if halves is not None:
                    scores = paired_control_scores(x, y, halves)
                    if np.isfinite(scores).all():
                        paired.append(scores)
                    else:
                        reasons["undefined_paired"] += 1
                else:
                    reasons["single_vehicle"] += 1
            elif vehicle_a.isdisjoint(vehicle_b):
                cross.append(original)
                if data.obs.iloc[a].plate == data.obs.iloc[b].plate:
                    reasons["distinct_group_same_plate"] += 1
            else:
                reasons["partially_shared"] += 1
        valid_same = np.asarray(same)[np.isfinite(same)]
        valid_cross = np.asarray(cross)[np.isfinite(cross)]
        shared, independent = np.mean(paired, axis=0) if paired else (np.nan, np.nan)
        records.append(dict(
            parent=parent, dose=dose, cell_line=line, profiles=len(condition),
            original_shared=np.mean(valid_same) if len(valid_same) else np.nan,
            original_disjoint=np.mean(valid_cross) if len(valid_cross) else np.nan,
            matched_pool_shared=shared, matched_pool_disjoint=independent,
            shared_minus_disjoint=shared - independent,
            original_shared_pairs=len(same), original_disjoint_pairs=len(cross),
            paired_defined=len(paired), undefined_original=len(same)+len(cross)-len(valid_same)-len(valid_cross),
            **{name: reasons[name] for name in (
                "undefined_paired", "single_vehicle", "distinct_group_same_plate", "partially_shared")},
        ))
    return pd.DataFrame(records)


def estimate(frame, metrics):
    """Equal condition/profile weight within parent, then equal parent weight."""
    result = []
    rng = np.random.default_rng(SEED)
    for dose, group in frame.groupby("dose", sort=True):
        for metric in metrics:
            defined = group.dropna(subset=[metric])
            parents = defined.groupby("parent")[metric].mean().to_numpy()
            lower = upper = score = np.nan
            if len(parents):
                score = float(parents.mean())
            if len(parents) > 1:
                samples = rng.choice(parents, size=(DRAWS, len(parents)), replace=True)
                lower, upper = np.quantile(samples.mean(axis=1), [.025, .975])
            result.append(dict(dose=dose, metric=metric, score=score, lower=lower, upper=upper,
                               defined=len(defined), n_parent=len(parents), total=len(group),
                               total_parent=group["parent"].nunique()))
    return pd.DataFrame(result)


def residual_scores(data, run, train, validation, accepted_manifest):
    fit_path = run / "manifest.json"
    validation_path = run / "validation_manifest.json"
    if sha256(fit_path) != accepted_manifest:
        raise ValueError(f"unaccepted Step 03 fit manifest: {run}")
    fit = json.loads(fit_path.read_text())
    if sha256(validation_path) != fit["outputs"]["validation_manifest.json"]:
        raise ValueError(f"validation manifest differs from accepted fit: {run}")
    validation_manifest = json.loads(validation_path.read_text())
    if validation_manifest["input_manifest_sha256"] != sha256(data.directory / "manifest.json"):
        raise ValueError(f"validation prediction source mismatch: {run}")
    for name in ("validation_prediction.npy", "validation_identity.csv"):
        if sha256(run / name) != validation_manifest["outputs"][name]:
            raise ValueError(f"validation digest mismatch: {run / name}")
    target = data.obs.iloc[validation]
    identity = pd.read_csv(run / "validation_identity.csv")
    if not identity.equals(target[IDENTITY].reset_index(drop=True)):
        raise ValueError(f"validation row order mismatch: {run}")
    prediction = np.load(run / "validation_prediction.npy", mmap_mode="r")
    if prediction.shape != (len(validation), data.expression.shape[1]):
        raise ValueError(f"validation prediction shape mismatch: {run}")

    representative = target.loc[target.cell_line.isin(data.obs.iloc[train].cell_line)]
    representative = representative.drop_duplicates("cell_line")
    rows = data.obs.index.get_indexer(representative.index)
    reference = dict(zip(representative.cell_line, mean_cell(data, train, rows), strict=True))
    records = []
    for i, row in enumerate(validation):
        obs = target.iloc[i]
        truth = data.expression[row]
        control = data.control[row]
        pred = prediction[i]
        original = pcc(pred - control, truth - control)
        mean = reference.get(obs.cell_line)
        residual = pcc(pred - mean, truth - mean) if mean is not None else np.nan
        records.append(dict(parent=obs.drug_parent_inchikey, dose=obs.dose_um,
                            original_delta=original, residual=residual,
                            residual_minus_delta=residual - original,
                            no_training_cell=mean is None))
    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()
    pair = frozen_pair(args.inputs)
    acceptance = STEP / "log/round1-baselines-acceptance.json"
    accepted = json.loads(acceptance.read_text())["scientific_artifacts"]
    parts = {name: [] for name in ("repeat_summary", "repeat_coverage", "residual_summary")}
    provenance = dict(input_pair_sha256=pair, source_manifests={}, validation_manifests={},
                      acceptance_sha256=sha256(acceptance),
                      training_seed=42, split_seed=226, bootstrap_seed=SEED, bootstrap_draws=DRAWS,
                      checkpoint="Step 03 latent-embedding best validation MSE")
    for source in SOURCE_POPULATIONS:
        data = Inputs(args.inputs / source)
        provenance["source_manifests"][source] = sha256(data.directory / "manifest.json")
        development = sorted(set(data.splits.loc[data.splits.partition.eq("validation"), "row"]))
        pools = control_pools(data, development)
        for scenario, fold in itertools.product(SCENARIOS, range(3)):
            val = data.rows(scenario, fold, "validation")
            train = data.rows(scenario, fold, "train")
            run = args.matrix / "runs" / source / scenario / f"fold-{fold}" / "latent-embedding"
            provenance["validation_manifests"][f"{source}/{scenario}/{fold}"] = sha256(
                run / "validation_manifest.json")
            repeat = repeat_scores(data, val, pools)
            residual = residual_scores(data, run, train, val,
                                       accepted[str(run.relative_to(args.matrix) / "manifest.json")])
            repeat_metrics = ["original_shared", "original_disjoint", "matched_pool_shared",
                              "matched_pool_disjoint", "shared_minus_disjoint"]
            residual_metrics = ["original_delta", "residual", "residual_minus_delta"]
            coverage = repeat.assign(single_profile=repeat.profiles == 1)
            for name, frame in (
                    ("repeat_summary", estimate(repeat, repeat_metrics)),
                    ("repeat_coverage", coverage.groupby("dose")[
                        ["profiles", "single_profile", "original_shared_pairs", "original_disjoint_pairs",
                         "paired_defined", "undefined_original", "undefined_paired", "single_vehicle",
                         "distinct_group_same_plate", "partially_shared"]].sum().reset_index()),
                    ("residual_summary", estimate(residual, residual_metrics))):
                parts[name].append(frame.assign(source=source, scenario=scenario, fold=fold))
    args.destination.mkdir(parents=True)
    for name, frames in parts.items():
        pd.concat(frames, ignore_index=True).to_csv(
            args.destination / f"{name}.csv", index=False, float_format="%.8f")
    (args.destination / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
