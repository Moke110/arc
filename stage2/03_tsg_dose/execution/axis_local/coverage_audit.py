"""Audit #261 from frozen metadata and accepted per-profile Delta PCC scores."""
import argparse
from hashlib import file_digest
import json
from pathlib import Path

import numpy as np
import pandas as pd

SOURCES = ("tahoe_bulk", "tenv")
SCENARIOS = ("cold-drug", "cold-drug-dose")
METHODS = ("Mean_cell", "no-dose", "latent-embedding")
IDENTITY = ["profile_id", "drug_parent", "dose_um", "cell_line"]
SIMILARITY_BOUNDS = [0, .25, .50, .75, 1]
SIMILARITY_LABELS = ["0–0.25", ">0.25–0.50", ">0.50–0.75", ">0.75–1"]


def file_digest_in_path(path):
    with path.open("rb") as stream:
        return file_digest(stream, "sha256").hexdigest()


def verified(path, expected):
    if file_digest_in_path(path) != expected:
        raise ValueError(f"Digest mismatch: {path}")
    return path


def read_json(path):
    return json.loads(path.read_text())


def paired_scores(scores):
    finite = np.isfinite(scores[list(METHODS)]).all(axis=1)
    table = scores.loc[finite, ["drug_parent", *METHODS]].copy()
    for method in METHODS[1:]:
        table[f"{method}_minus_Mean_cell"] = table[method] - table["Mean_cell"]
    table["latent-embedding_minus_no-dose"] = table["latent-embedding"] - table["no-dose"]
    return table


def summarize_scores(scores, seed, draws):
    table = paired_scores(scores)
    counts = table.groupby("drug_parent").size()
    parent = table.groupby("drug_parent").sum()
    profile_counts = counts.to_numpy()
    totals = parent.to_numpy()
    parent_means = totals / profile_counts[:, None]
    rng = np.random.default_rng(seed)
    profile_draws, parent_draws = [], []
    n = len(parent)
    if n > 1:
        for start in range(0, draws, 128):
            chosen = rng.integers(n, size=(min(128, draws - start), n))
            profile_draws.append(totals[chosen].sum(axis=1) / profile_counts[chosen].sum(axis=1)[:, None])
            parent_draws.append(parent_means[chosen].mean(axis=1))
    distributions = {}
    rows = []
    points = {
        "profile": totals.sum(axis=0) / max(1, profile_counts.sum()),
        "parent": parent_means.mean(axis=0) if n else np.full(len(parent.columns), np.nan),
    }
    for weight, samples in (("profile", profile_draws), ("parent", parent_draws)):
        replicates = np.concatenate(samples) if samples else np.full((draws, len(parent.columns)), np.nan)
        for col, metric in enumerate(parent.columns):
            distributions[weight, metric] = replicates[:, col]
            interval = np.quantile(replicates[:, col], [.025, .975])
            rows.append(dict(weight=weight, metric=metric, value=points[weight][col] if n else np.nan,
                             ci_low=interval[0], ci_high=interval[1], n_profiles=len(table),
                             n_parents=n, n_excluded=len(scores) - len(table)))
    return pd.DataFrame(rows), distributions


def nearest_similarity(fingerprints, training_parents):
    values = fingerprints.to_numpy(dtype=np.float32)
    train = fingerprints.loc[sorted(training_parents)].to_numpy(dtype=np.float32)
    if not len(train) or not np.isin(values, [0, 1]).all() or np.any(values.sum(axis=1) == 0):
        raise ValueError("Need binary, nonempty ECFP4 vectors and at least one training parent")
    intersection = values @ train.T
    union = values.sum(axis=1)[:, None] + train.sum(axis=1) - intersection
    return pd.Series((intersection / union).max(axis=1), index=fingerprints.index)


def inspect_coverage(rows):
    conditions = rows.groupby(["drug_parent", "dose_um", "cell_line"]).size()
    similarity = rows.drop_duplicates("drug_parent").nearest_fit_ecfp4
    return dict(n_profiles=len(rows), n_parents=rows.drug_parent.nunique(),
                n_doses=rows.dose_um.nunique(), n_cells=rows.cell_line.nunique(),
                n_plates=rows.plate.nunique(), n_conditions=len(conditions),
                n_replicated_conditions=int((conditions > 1).sum()),
                profiles_per_condition_p50=conditions.median(),
                profiles_per_condition_max=conditions.max(),
                nearest_fit_ecfp4_p10=similarity.quantile(.1),
                nearest_fit_ecfp4_p50=similarity.quantile(.5),
                nearest_fit_ecfp4_p90=similarity.quantile(.9))


def check_scores(path, expected, pcc_rows):
    frame = pd.read_csv(verified(path, expected), dtype={"profile_id": str, "cell_line": str},
                        float_precision="round_trip", usecols=[*IDENTITY, "pcc"])
    pd.testing.assert_frame_equal(frame[IDENTITY].reset_index(drop=True),
                                  pcc_rows[IDENTITY].reset_index(drop=True), check_dtype=False)
    return frame.pcc.to_numpy()


def assess(args):
    step = args.step
    checkpoints = step / "execution/checkpoints"
    inputs = checkpoints / "round1-axis-local-inputs"
    baselines = checkpoints / "round1-baselines-f32/runs"
    formal = checkpoints / "round1-axis-local-formal"
    accepted = read_json(step / "log/round1-baselines-acceptance.json")
    gray = read_json(step / "log/baselines-gray-acceptance.json")
    pooled = pd.read_csv(verified(step / "figure/data/pooled.csv", gray["files"]["data/pooled.csv"]))
    pair = read_json(verified(inputs / "manifest.json",
                              read_json(step / "log/247-axis-local-smoke.json")["pair_sha256"]))
    coverage, replicate_counts, conditions, development, results, per_parent = [], [], [], [], [], []
    bootstraps = {}
    current_seed = np.random.SeedSequence(args.seed).spawn(len(SOURCES) * len(SCENARIOS))
    for i, source in enumerate(SOURCES):
        directory = inputs / source
        manifest = read_json(verified(directory / "manifest.json", pair["sources"][source]))
        files = {name: verified(directory / name, manifest["outputs"][name]) for name in
                 ("all_obs.csv", "splits.csv", "parents.csv", "morgan.npy")}
        obs = pd.read_csv(files["all_obs.csv"], dtype={"profile_id": str, "cell_line": str},
                          keep_default_na=False)
        obs = obs.loc[obs.role.eq("treated")].rename(columns={"drug_parent_inchikey": "original_parent",
                                                            "dose_um": "original_dose"})
        if obs.profile_id.duplicated().any() or obs[["source_item", "source_row"]].duplicated().any():
            raise ValueError("Original treated profile/source-row identity is not unique")
        ledger = pd.read_csv(files["splits.csv"]).merge(
            obs[["profile_id", "original_parent", "original_dose", "cell_line", "plate"]],
            on="profile_id", how="left", validate="many_to_one")
        if (not ledger.drug_parent.eq(ledger.original_parent).all()
                or not ledger.dose_um.eq(ledger.original_dose).all()
                or ledger.duplicated(["scenario", "fold", "profile_id"]).any()):
            raise ValueError("Split rows differ from original treated identities")
        parents = pd.read_csv(files["parents.csv"]).drug_parent_inchikey
        fingerprints = pd.DataFrame(np.load(files["morgan.npy"], allow_pickle=False), index=parents)
        for j, scenario in enumerate(SCENARIOS):
            tests = []
            for fold in range(3):
                split = ledger.loc[ledger.scenario.eq(scenario) & ledger.fold.eq(fold)].copy()
                train_parents = set(split.loc[split.partition.eq("train"), "drug_parent"])
                split["nearest_fit_ecfp4"] = split.drug_parent.map(nearest_similarity(fingerprints, train_parents))
                if len(split) != len(obs) or not {"train", "test", "validation"}.issubset(set(split.partition)):
                    raise ValueError("Incomplete frozen split")
                for partition, part in split.groupby("partition", sort=True):
                    key = dict(source=source, scenario=scenario, fold=fold, partition=partition)
                    coverage.append(dict(**key, **inspect_coverage(part)))
                    for replicates, count in part.groupby(
                            ["drug_parent", "dose_um", "cell_line"]).size().value_counts().sort_index().items():
                        replicate_counts.append(dict(**key, profiles_per_condition=replicates, n_conditions=count))
                    for (dose, cell), group in part.groupby(["dose_um", "cell_line"]):
                        conditions.append(dict(**key, dose_um=dose, cell_line=cell,
                                               **inspect_coverage(group)))
                test = split.loc[split.partition.eq("test")].copy()
                validation = split.loc[split.partition.eq("validation")].copy()
                origin = next(record for record in accepted["origins"]
                              if (record["source"], record["scenario"], record["fold"], record["method"])
                              == (source, scenario, fold, "Mean_cell"))
                baseline = baselines / source / scenario / str(fold) / "Mean_cell"
                baseline_receipt = read_json(verified(baseline / "manifest.json", origin["manifest_sha256"]))
                test["Mean_cell"] = check_scores(baseline / "test_per_sample.csv",
                                                baseline_receipt["outputs"]["test_per_sample.csv"], test)
                for method in METHODS[1:]:
                    rel = Path("runs") / source / scenario / f"fold-{fold}" / method
                    run = formal / rel
                    run_receipt = read_json(verified(run / "manifest.json",
                                                    accepted["scientific_artifacts"][str(rel / "manifest.json")]))
                    for partition, part in (("test", test), ("validation", validation)):
                        artifact = f"{partition}_per_sample.csv"
                        part[method] = check_scores(run / artifact, run_receipt["outputs"][artifact], part)
                tests.append(test)
                validation["similarity_bin"] = pd.cut(validation.nearest_fit_ecfp4, SIMILARITY_BOUNDS,
                                                       include_lowest=True, labels=SIMILARITY_LABELS)
                validation["gain"] = validation["latent-embedding"] - validation["no-dose"]
                for label in SIMILARITY_LABELS:
                    group = validation.loc[validation.similarity_bin.eq(label)]
                    by_parent = group.groupby("drug_parent")["gain"].mean()
                    development.append(dict(source=source, scenario=scenario, fold=fold,
                                            similarity_bin=label, n_profiles=len(group), n_parents=len(by_parent),
                                            n_finite_parents=int(by_parent.notna().sum()),
                                            gain_parent_equal=by_parent.mean(),
                                            similarity_low=group.nearest_fit_ecfp4.min(),
                                            similarity_high=group.nearest_fit_ecfp4.max()))
            scores = pd.concat(tests, ignore_index=True)
            if scores.profile_id.duplicated().any():
                raise ValueError("A published profile appears in multiple outer folds")
            for method in METHODS:
                published = pooled.loc[(pooled.source == source) & (pooled.scenario == scenario)
                                       & (pooled.seed == 42) & (pooled.method == method)]
                if (len(published) != 1 or published.n_samples.item() != len(scores)
                        or not np.isclose(published.pcc.item(), scores[method].mean(), rtol=0, atol=1e-12)):
                    raise ValueError(f"Cannot reconcile accepted pooled result: {source} {scenario} {method}")
            table, draws = summarize_scores(scores, current_seed[i * len(SCENARIOS) + j], args.draws)
            results.append(table.assign(source=source, scenario=scenario))
            bootstraps[source, scenario] = draws
            paired = paired_scores(scores)
            by_parent = paired.groupby("drug_parent").mean()
            by_parent["n_profiles"] = paired.groupby("drug_parent").size()
            per_parent.append(by_parent.reset_index().assign(source=source, scenario=scenario))
    scores = pd.concat(results, ignore_index=True)
    point = scores.set_index(["source", "scenario", "weight", "metric"]).value
    cross_domain = []
    for scenario in SCENARIOS:
        for (weight, metric), samples in bootstraps["tahoe_bulk", scenario].items():
            difference = samples - bootstraps["tenv", scenario][weight, metric]
            low, high = np.quantile(difference, [.025, .975])
            cross_domain.append(dict(scenario=scenario, weight=weight, metric=metric,
                                     tahoe_minus_tenv=point["tahoe_bulk", scenario, weight, metric]
                                     - point["tenv", scenario, weight, metric], ci_low=low, ci_high=high))
    validation = pd.DataFrame(development)
    folds = pd.DataFrame(conditions)
    dev_conditions = folds.loc[folds.partition.eq("validation")]
    support = []
    for (scenario, fold), group in dev_conditions.groupby(["scenario", "fold"]):
        domain = {name: group.loc[group.source.eq(name)] for name in SOURCES}
        intersection = set(zip(domain["tahoe_bulk"].dose_um, domain["tahoe_bulk"].cell_line)) & set(
            zip(domain["tenv"].dose_um, domain["tenv"].cell_line))
        bands = validation.loc[(validation.scenario == scenario) & (validation.fold == fold)]
        n_by_bin = bands.pivot(index="similarity_bin", columns="source", values="n_parents")
        bins_in_both = [label for label in SIMILARITY_LABELS if n_by_bin.loc[label].gt(0).all()]
        support.append(dict(scenario=scenario, fold=fold,
                            occupied_similarity_bins_in_both=";".join(bins_in_both),
                            shared_exact_doses_um=";".join(map(str, sorted(
            set(domain["tahoe_bulk"].dose_um) & set(domain["tenv"].dose_um)))),
            shared_exact_cell_lines=";".join(sorted(set(domain["tahoe_bulk"].cell_line) & set(
                domain["tenv"].cell_line))), shared_exact_dose_cell_groups=len(intersection)))
    return {
        "coverage": pd.DataFrame(coverage), "replicate_counts": pd.DataFrame(replicate_counts),
        "condition_coverage": folds, "development_similarity": validation,
        "development_joint_support": pd.DataFrame(support), "published_parent": pd.concat(per_parent),
        "published_scores": scores, "published_domain_gap": pd.DataFrame(cross_domain),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", required=True, type=Path, help="Materialized Step 03 path")
    parser.add_argument("--output", required=True, type=Path, help="Fresh report directory")
    parser.add_argument("--seed", type=int, default=261)
    parser.add_argument("--draws", type=int, default=4000)
    args = parser.parse_args()
    if args.draws < 1:
        parser.error("--draws must be positive")
    tables = assess(args)
    args.output.mkdir(parents=True, exist_ok=False)
    for name, table in tables.items():
        table.to_csv(args.output / f"{name}.csv", index=False, float_format="%.12g")
    table_digests = {path.name: file_digest_in_path(path)
                     for path in sorted(args.output.glob("*.csv"))}
    (args.output / "run.json").write_text(json.dumps(
        dict(schema="stage2-map-261-coverage-audit-v1", bootstrap_seed=args.seed,
             bootstrap_draws=args.draws, split_seed=226, published_fit_seed=42, tables_sha256=table_digests), indent=2) + "\n")


if __name__ == "__main__":
    main()
