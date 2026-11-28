"""Independently audit Step 05 published tables using standard-library CSV readers."""

import argparse
import csv
from itertools import product
import json
import math
from pathlib import Path
import statistics

from run_module import fit_coordinates as coordinates


OUTCOMES = ("Δ PCC", "Δ RMSE", "Δ Spearman", "Δ E-distance")
COMPARISONS = (
    ("XΔ-XΔΔ − XX-XXΔ", "XΔ-XΔΔ", "XX-XXΔ"),
    ("XΔ-XΔ − XΔ-XΔΔ", "XΔ-XΔ", "XΔ-XΔΔ"),
)
METRICS = ("JAC_top200", "Δ PCC", "Δ Spearman", "Δ RMSE",
           "Δ E-distance", "Δ MMD", "Δ R²")


def read_table(path, keys):
    rows = {}
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            key = tuple(int(row[name]) if name in {"fold", "seed"} else row[name] for name in keys)
            if key in rows:
                raise ValueError(f"duplicate {path.name} row {key}")
            rows[key] = row
    return rows


def value(row, name="value"):
    raw = row[name]
    if not raw:
        return None
    parsed = float(raw)
    if not math.isfinite(parsed):
        raise ValueError(f"nonfinite defined value in {name}")
    return parsed


def same(found, expected):
    if found is None or expected is None:
        return found is expected
    return math.isclose(found, expected, rel_tol=1e-12, abs_tol=1e-12)


def check_statistics(row, values, label, field="value", sd="sample_sd", count="n_defined"):
    defined = [number for number in values if number is not None]
    mean = statistics.mean(defined) if defined else None
    sample_sd = statistics.stdev(defined) if len(defined) > 1 else None
    if (int(row[count]) != len(defined) or int(row["n_coordinates"]) != len(values)
            or not same(value(row, field), mean) or not same(value(row, sd), sample_sd)):
        raise ValueError(f"invalid aggregate {label}")


def audit(report, smoke):
    required = coordinates(smoke)
    arms = tuple(dict.fromkeys(shape for shape, _, _ in required))
    runs = tuple(dict.fromkeys((fold, seed) for _, fold, seed in required))
    folds = tuple(dict.fromkeys(fold for fold, _ in runs))
    score_rows = read_table(report / "fold-scores.csv", ("shape", "fold", "seed", "metric"))
    summary = read_table(report / "delta-metrics.csv", ("shape", "metric"))
    if (set(score_rows) != {(*key, metric) for key, metric in product(required, METRICS)}
            or set(summary) != set(product(arms, METRICS))):
        raise ValueError("scores or summaries differ from the registered roster")

    scores = {}
    for key, row in score_rows.items():
        score = value(row)
        counts = tuple(int(row[name]) for name in ("n_conditions", "n_eligible", "n_undefined"))
        if any(count < 0 for count in counts) or counts[0] != sum(counts[1:]) or bool(counts[1]) != (score is not None):
            raise ValueError(f"invalid eligibility at {key}")
        scores[key] = score
    for shape, metric in product(arms, METRICS):
        keys = [(shape, fold, seed, metric) for fold, seed in runs]
        differences = []
        for key in keys:
            score, control = scores[key], scores["XX-XX", *key[1:]]
            difference = None if score is None or control is None else score - control
            if not same(value(score_rows[key], "score_difference"), difference):
                raise ValueError(f"invalid reference difference at {key}")
            differences.append(difference)
        aggregate = summary[shape, metric]
        for component in ("eligible", "undefined"):
            expected = sum(int(score_rows[key][f"n_{component}"]) for key in keys)
            if int(aggregate[f"n_{component}_conditions"]) != expected:
                raise ValueError(f"invalid aggregate eligibility for {(shape, metric)}")
        check_statistics(aggregate, [scores[key] for key in keys], (shape, metric))
        check_statistics(aggregate, differences, (shape, metric), "score_difference",
                         "score_difference_sd", "n_paired")

    contrasts = read_table(report / "primary-contrasts.csv", ("contrast", "metric", "fold", "seed"))
    by_fold = read_table(report / "primary-contrast-by-fold.csv", ("contrast", "metric", "fold"))
    overall = read_table(report / "primary-contrast-summary.csv", ("contrast", "metric"))
    groups = {(comparison[0], metric) for comparison, metric in product(COMPARISONS, OUTCOMES)}
    if (set(overall) != groups
            or set(contrasts) != {(*group, *run) for group, run in product(groups, runs)}
            or set(by_fold) != {(*group, fold) for group, fold in product(groups, folds)}):
        raise ValueError("contrast tables differ from the registered roster")
    for (label, left, right), metric in product(COMPARISONS, OUTCOMES):
        differences = []
        for fold, seed in runs:
            key = label, metric, fold, seed
            row = contrasts[key]
            left_score, right_score = scores[left, fold, seed, metric], scores[right, fold, seed, metric]
            difference = None if left_score is None or right_score is None else left_score - right_score
            if ((row["left_shape"], row["right_shape"]) != (left, right)
                    or not same(value(row, "left_value"), left_score)
                    or not same(value(row, "right_value"), right_score)
                    or not same(value(row), difference)):
                raise ValueError(f"invalid paired difference at {key}")
            differences.append(difference)
        check_statistics(overall[label, metric], differences, (label, metric))
        for fold in folds:
            within_fold = [diff for (run_fold, _), diff in zip(runs, differences, strict=True) if run_fold == fold]
            check_statistics(by_fold[label, metric, fold], within_fold, (label, metric, fold))
    return {"status": "PASS", "fold_scores": len(score_rows), "equal_coordinate_aggregates": len(summary),
            "paired_contrast_scores": len(contrasts), "foldwise_contrast_aggregates": len(by_fold),
            "overall_contrast_aggregates": len(overall)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    print(json.dumps(audit(args.report, args.smoke), sort_keys=True))


if __name__ == "__main__":
    main()
