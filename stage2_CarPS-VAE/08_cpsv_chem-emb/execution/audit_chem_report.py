"""Check Step 08 pairing and aggregates using the standard library."""

from itertools import product

from chem_report import METRICS
from chem_runs import REFERENCE, coordinates
from audit_measured_scores import check_statistics, read_table, same, value


def audit(report, smoke=False):
    roster = coordinates(smoke)
    arms = tuple(dict.fromkeys(arm for arm, _, _ in roster))
    folds = tuple(dict.fromkeys(fold for _, fold, _ in roster))
    score_rows = read_table(report / "fold-scores.csv", ("arm", "fold", "seed", "metric"))
    by_fold = read_table(report / "paired-by-fold.csv", ("arm", "metric", "fold"))
    summary = read_table(report / "paired-summary.csv", ("arm", "metric"))
    if (set(score_rows) != {(*key, metric) for key, metric in product(roster, METRICS)}
            or set(by_fold) != set(product(arms, METRICS, folds))
            or set(summary) != set(product(arms, METRICS))):
        raise ValueError("Step 08 score and summary tables must cover the registered roster")

    groups = {}
    for key, row in score_rows.items():
        arm, fold, seed, metric = key
        score = value(row)
        reference = value(score_rows[REFERENCE, fold, seed, metric])
        difference = None if score is None or reference is None else score - reference
        total, eligible, undefined = (int(row[name]) for name in ("n_conditions", "n_eligible", "n_undefined"))
        if (min(eligible, undefined) < 0 or total != eligible + undefined
                or bool(eligible) != (score is not None)
                or not same(value(row, "kpgt_value"), reference)
                or not same(value(row, "score_difference"), difference)):
            raise ValueError(f"invalid KPGT pair or scoring eligibility at {key}")
        for group in ((arm, metric), (arm, metric, fold)):
            groups.setdefault(group, []).append((score, reference, difference, eligible, undefined))

    for group, rows in groups.items():
        aggregate = (summary if len(group) == 2 else by_fold)[group]
        scores, references, differences, eligible, undefined = zip(*rows)
        check_statistics(aggregate, scores, group)
        check_statistics(aggregate, differences, group, "score_difference", "score_difference_sd", "n_paired")
        expected_reference = [reference for reference in references if reference is not None]
        if (not same(value(aggregate, "kpgt_value"),
                     sum(expected_reference) / len(expected_reference) if expected_reference else None)
                or int(aggregate["n_eligible_conditions"]) != sum(eligible)
                or int(aggregate["n_undefined_conditions"]) != sum(undefined)):
            raise ValueError(f"invalid reference mean or aggregate eligibility at {group}")
    return {"status": "PASS", "score_rows": len(score_rows), "foldwise_aggregates": len(by_fold),
            "overall_aggregates": len(summary)}
