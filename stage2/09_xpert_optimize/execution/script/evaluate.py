"""Seal validation promotions and independently replay the Round 1 test figures."""

import argparse
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from round_contract import DEFAULT_WORK, PROMOTION_SCHEMA, checked_payload, read_json, sha256, verify_committed_promotion, write_json
from round_inputs import CASE_COLUMNS, FULL_VIEW, MODEL_VIEWS, REFERENCES, SHARED_VIEW, STEP, RoundInputs, Significance
from round_metrics import CONTINUOUS, PR_FIELDS, average_precision, continuous, pooled_pr, paired_promotions, project_delta


RECORD = STEP / "execution/round1_record"
DEFAULT_ROUND = STEP / "output/round1"
DEFAULT_CORPUS = STEP / "output/corpus"
FIGURES = ("xpert_opti_r1_metrics.png", "xpert_opti_r1_PR.png")
TABLES = ("metrics.csv", "per_group.csv", "pr.csv")
METRIC_KEYS = ("auprc_full", "auprc_g", *CONTINUOUS, *(f"native_{key}" for key in CONTINUOUS))
BLOCK = 64


class Groups:
    def __init__(self, cases):
        self.positions, self.names = pd.factorize(cases.group_id, sort=False)
        self.counts = np.bincount(self.positions)

    def blocks(self):
        for start in range(0, len(self.positions), BLOCK):
            part = slice(start, start + BLOCK)
            yield part, self.positions[part]

    def empty(self, genes):
        return np.zeros((len(self.names), genes), dtype=np.float64)

    def means(self, values):
        return values / self.counts[:, None]


def observed_means(cohort, groups):
    result = groups.empty(len(cohort.genes))
    for part, positions in groups.blocks():
        expr = cohort.treated[part]
        control = cohort.controls[cohort.control_index[part]]
        if not np.isfinite(expr).all() or not np.isfinite(control).all():
            raise ValueError("non-finite observed expression or matched control")
        np.add.at(result, positions, expr.astype(np.float64) - control)
    return groups.means(result)


def predicted_means(delta, native, shared, groups, shared_columns):
    native_sums, shared_sums = groups.empty(len(native.genes)), groups.empty(len(shared.genes))
    for part, positions in groups.blocks():
        values = delta[part]
        if not np.isfinite(values).all():
            raise ValueError("non-finite held-out prediction")
        controls = native.controls[native.control_index[part]]
        shared_controls = shared.controls[shared.control_index[part]]
        np.add.at(native_sums, positions, values)
        np.add.at(shared_sums, positions, project_delta(values, controls, shared_columns, shared_controls))
    return groups.means(native_sums), groups.means(shared_sums)


def reference_delta(inputs, view):
    train = inputs.cohort(view, "train")
    grouped = observed_means(train, Groups(train.cases))
    return grouped.mean(axis=0)


def pr_fields(name, ranking, labels, all_positive, n_genes, group_eligible, row):
    positive = int(all_positive)
    result = average_precision(ranking, labels, positive) if group_eligible else None
    row[f"auprc_{name}"] = result.ap if result else None
    row[f"terminal_recall_{name}"] = float(result.recall[-1]) if result else None
    row[f"significant_genes_{name}"] = positive if group_eligible else None
    row[f"prevalence_{name}"] = positive / n_genes if group_eligible and n_genes else None
    row[f"ranked_genes_{name}"] = len(ranking)
    row[f"auprc_status_{name}"] = ("untestable_group" if not group_eligible
                                  else "defined" if result else "zero_significant")
    return result


def average_defined(rows, key):
    defined = [row[key] for row in rows if row[key] is not None]
    return (float(np.mean(defined)) if defined else None), len(defined)


def score_model(model, native, shared, groups, labels, obs_native, pred_native, obs_shared, pred_shared,
                full_count, minutes):
    label_indices = labels.genes.set_index("hgnc_id").index.get_indexer(native.genes.hgnc_id)
    if np.any(label_indices < 0):
        raise ValueError(f"{model}: prediction axis outside significance family")
    full_columns = label_indices < full_count
    shared_genes = len(shared.genes)
    rows, curves = [], {"full": [], "g": []}
    for i, group in enumerate(groups.names):
        metadata, significant, testable = labels.row(group)
        group_eligible = metadata.status == "eligible"
        row = {"model": model, "group_id": group, "wells": int(groups.counts[i]),
               "significance_status": metadata.status, "significance_reason": metadata.get("reason", ""),
               "untestable_family_genes": int((~testable).sum())}
        direction = labels.log2_fold_change[labels.index[group]]
        if np.any(significant & (~np.isfinite(direction) | (direction == 0))):
            raise ValueError(f"invalid significant direction for {group}")
        native_truth = significant[label_indices]
        signed_truth = np.where(native_truth, np.sign(direction[label_indices]), 0)
        for name, columns, positive, n_genes in (
            ("full", full_columns, significant[:full_count].sum(), full_count),
            ("g", slice(None), native_truth.sum(), len(native.genes)),
        ):
            result = pr_fields(name, pred_native[i, columns], signed_truth[columns], positive,
                               n_genes, group_eligible, row)
            if result:
                curves[name].append(result)
        row.update(continuous(pred_shared[i], obs_shared[i]))
        row.update({f"native_{name}": value for name, value in continuous(pred_native[i], obs_native[i]).items()})
        rows.append(row)

    summary = {"model": model, "status": "complete", "reason": "", "gene_count": len(native.genes),
               "cm_covered_genes": int(full_columns.sum()), "cm_coverage": float(full_columns.sum() / full_count),
               "shared_gene_count": shared_genes, "continuous_provisional": shared_genes < 200,
               "groups": len(rows), "significance_eligible_groups": int(labels.groups.status.eq("eligible").sum()),
               "training_minutes": minutes}
    for key in (*METRIC_KEYS, "terminal_recall_full", "terminal_recall_g", "prevalence_full", "prevalence_g"):
        summary[key], summary[f"{key}_defined_groups"] = average_defined(rows, key)
        summary[f"{key}_undefined_groups"] = len(rows) - summary[f"{key}_defined_groups"]
    pr = []
    for name, group_curves in curves.items():
        frame = pooled_pr(group_curves)
        if not frame.empty:
            frame.insert(0, "metric", name)
            frame.insert(0, "model", model)
            pr.append(frame)
    return summary, rows, pr


def unavailable(model, native, shared, groups, full_genes, reason, run):
    covered = int(native.genes.hgnc_id.isin(full_genes.hgnc_id).sum())
    result = {"model": model, "status": "incomplete", "reason": reason, "gene_count": len(native.genes),
              "cm_covered_genes": covered, "cm_coverage": covered / len(full_genes),
              "shared_gene_count": len(shared.genes), "continuous_provisional": len(shared.genes) < 200,
              "groups": len(groups.names), "significance_eligible_groups": None,
              "training_minutes": run["training_minutes"] if run else None}
    for key in (*METRIC_KEYS, "terminal_recall_full", "terminal_recall_g", "prevalence_full", "prevalence_g"):
        result[key], result[f"{key}_defined_groups"] = None, 0
        result[f"{key}_undefined_groups"] = len(groups.names)
    return result


def score(inputs, role, destination, promotion_sha=None):
    statuses = inputs.statuses(role)
    destination.mkdir()
    shared = inputs.cohort(SHARED_VIEW, role)
    full = inputs.cohort(FULL_VIEW, role)
    pd.testing.assert_frame_equal(full.cases[CASE_COLUMNS], shared.cases[CASE_COLUMNS])
    groups = Groups(full.cases)
    labels = Significance(inputs.root / "labels" / role, role, inputs.significance_genes, groups.names)
    expected_provenance = {
        "gene_axes_manifest_sha256": inputs.axis_sha,
        "predictor_manifest_sha256": inputs.views[FULL_VIEW].digest,
        "evaluator_manifest_sha256": inputs.evaluator_digests[FULL_VIEW],
        "landmark_predictor_manifest_sha256": inputs.views[MODEL_VIEWS["XPert_base"]].digest,
        "native_sha256": inputs.axis_manifest["source_hashes"]["tenv.h5ad"],
    }
    if any(labels.release.manifest[key] != digest for key, digest in expected_provenance.items()):
        raise ValueError("significance labels differ from the sealed Round 1 cohort")
    if role == "test" and labels.release.manifest["promotion_sha256"] != promotion_sha:
        raise ValueError("test significance labels do not bind the committed validation promotion")
    shared_observed, full_observed = observed_means(shared, groups), observed_means(full, groups)
    full_count = len(inputs.full_genes)
    reports, group_rows, pr = [], [], []
    evidence = {"view_evaluators": {SHARED_VIEW: shared.binding(), FULL_VIEW: full.binding()}, "fits": {}}

    def record(model, native, observed, predictions, shared_predictions, minutes):
        result, rows, curves = score_model(
            model, native, shared, groups, labels, observed, predictions, shared_observed, shared_predictions,
            full_count, minutes,
        )
        reports.append(result)
        group_rows.extend(rows)
        pr.extend(curves)

    for reference in REFERENCES:
        native_values = np.zeros(full_count) if reference == "no_change" else reference_delta(inputs, FULL_VIEW)
        shared_values = np.zeros(len(shared.genes)) if reference == "no_change" else reference_delta(inputs, SHARED_VIEW)
        record(reference, full, full_observed, np.broadcast_to(native_values, full_observed.shape),
               np.broadcast_to(shared_values, shared_observed.shape), 0.)
        reports[-1].update(fit_attempt_wall_minutes=0., fit_attempts=0, timed_fit_attempts=0)

    for model, view in MODEL_VIEWS.items():
        native = full if view == FULL_VIEW else inputs.cohort(view, role)
        pd.testing.assert_frame_equal(native.cases[CASE_COLUMNS], full.cases[CASE_COLUMNS])
        evidence["view_evaluators"][view] = native.binding()
        predictions, run, binding, reason = inputs.prediction(model, role, native, statuses[model])
        evidence["fits"][model] = binding
        if predictions is None:
            reports.append(unavailable(model, native, shared, groups, inputs.full_genes, reason, run))
        else:
            shared_columns = native.genes.set_index("hgnc_id").index.get_indexer(shared.genes.hgnc_id)
            if np.any(shared_columns < 0):
                raise ValueError(f"{model}: native axis misses shared genes")
            observed = full_observed if view == FULL_VIEW else observed_means(native, groups)
            predictions, shared_predictions = predicted_means(predictions, native, shared, groups, shared_columns)
            record(model, native, observed, predictions, shared_predictions, run["training_minutes"])
        attempts = binding["fit_attempts"]
        known = [item["seconds"] for item in attempts if item["seconds"] is not None]
        reports[-1].update(fit_attempt_wall_minutes=sum(known) / 60 if known else None,
                           fit_attempts=len(attempts), timed_fit_attempts=len(known))

    pd.DataFrame(reports).to_csv(destination / TABLES[0], index=False)
    pd.DataFrame(group_rows).to_csv(destination / TABLES[1], index=False)
    curve = pd.concat(pr, ignore_index=True) if pr else pd.DataFrame(columns=("model", "metric", *PR_FIELDS))
    curve.to_csv(destination / TABLES[2], index=False)
    manifest = {"schema": "step09-round1-score-v2", "direction_match": "DESeq2 log2 fold change",
                "pr_aggregation": "fixed-group micro at common absolute-delta thresholds",
                "ap_aggregation": "macro of exact per-group AP", "role": role, "fold": 0, "seed": 17,
                "plan_sha256": inputs.plan_sha, "gene_axes_manifest_sha256": inputs.axis_sha,
                "labels_sha256": labels.release.digest, "promotion_sha256": promotion_sha, **evidence,
                "files": {name: sha256(destination / name) for name in TABLES}}
    write_json(destination / "manifest.json", manifest)
    return reports, manifest


def require_same_release(published, replayed, names):
    for name in names:
        if sha256(published / name) != sha256(replayed / name):
            raise ValueError(f"independent replay differs: {published / name}")


def promotion_record(inputs, validation, results, manifest):
    return {"schema": PROMOTION_SCHEMA, "fold": 0, "seed": 17, "plan_sha256": inputs.plan_sha,
            "validation_manifest_sha256": sha256(validation / "manifest.json"),
            "validation_labels_sha256": manifest["labels_sha256"],
            "validation_metrics_sha256": manifest["files"]["metrics.csv"],
            "comparisons": paired_promotions(results)}


def seal(inputs, promotion_path):
    if (inputs.root / "test_scores").exists() or promotion_path.exists():
        raise FileExistsError("Round 1 already has a promotion seal or test scores")
    published = inputs.root / "validation_scores"
    with tempfile.TemporaryDirectory(dir=inputs.root, prefix=".validation-score-") as directory:
        temporary = Path(directory) / "validation_scores"
        results, manifest = score(inputs, "validation", temporary)
        promotion = promotion_record(inputs, temporary, results, manifest)
        if published.exists():
            require_same_release(published, temporary, (*TABLES, "manifest.json"))
        else:
            temporary.rename(published)
    promotion_path.parent.mkdir(parents=True, exist_ok=True)
    staged = promotion_path.with_suffix(".tmp")
    write_json(staged, promotion)
    staged.replace(promotion_path)
    return promotion


def validate_promotion(inputs, promotion_path):
    promotion = verify_committed_promotion(promotion_path)
    if promotion["plan_sha256"] != inputs.plan_sha:
        raise ValueError("committed promotion belongs to a different Round 1 plan")
    published = inputs.root / "validation_scores"
    if sha256(published / "manifest.json") != promotion["validation_manifest_sha256"]:
        raise ValueError("validation scoring manifest differs from the committed promotion")
    manifest = read_json(published / "manifest.json")
    for name in TABLES:
        checked_payload(published, manifest, name)
    with tempfile.TemporaryDirectory(dir=inputs.root, prefix=".validation-replay-") as directory:
        temporary = Path(directory) / "validation_scores"
        results, replayed_manifest = score(inputs, "validation", temporary)
        require_same_release(published, temporary, (*TABLES, "manifest.json"))
        if promotion_record(inputs, temporary, results, replayed_manifest) != promotion:
            raise ValueError("recomputed validation promotion differs from the committed seal")
    return sha256(promotion_path)


def test_release(inputs, promotion_path, destination):
    promotion_sha = validate_promotion(inputs, promotion_path)
    scores = destination / "test_scores"
    _, manifest = score(inputs, "test", scores, promotion_sha)
    from round_figures import render

    render(scores, destination, FIGURES)
    manifest["figures"] = {name: sha256(destination / name) for name in FIGURES}
    write_json(scores / "manifest.json", manifest)


def test(inputs, promotion_path):
    if (inputs.root / "test_scores").exists():
        raise FileExistsError("test scores are already published")
    with tempfile.TemporaryDirectory(dir=inputs.root, prefix=".test-score-") as directory:
        path = Path(directory)
        test_release(inputs, promotion_path, path)
        for name in FIGURES:
            (path / name).replace(inputs.root / name)
        (path / "test_scores").rename(inputs.root / "test_scores")


def replay(inputs, promotion_path):
    with tempfile.TemporaryDirectory(dir=inputs.root, prefix=".round1-replay-") as directory:
        path = Path(directory)
        test_release(inputs, promotion_path, path)
        require_same_release(inputs.root / "test_scores", path / "test_scores", (*TABLES, "manifest.json"))
        require_same_release(inputs.root, path, FIGURES)


def revision_release(inputs, promotion_path, destination):
    """Rescore frozen predictions without retroactively changing promotion."""
    promotion = verify_committed_promotion(promotion_path)
    if promotion["plan_sha256"] != inputs.plan_sha:
        raise ValueError("historical promotion belongs to a different plan")
    legacy = inputs.root / "legacy_unsigned"
    validation = legacy / "validation_scores" if legacy.exists() else inputs.root / "validation_scores"
    if sha256(validation / "manifest.json") != promotion["validation_manifest_sha256"]:
        raise ValueError("historical validation seal changed")
    old_manifest = read_json(validation / "manifest.json")
    for name in TABLES:
        checked_payload(validation, old_manifest, name)
    promotion_sha = sha256(promotion_path)
    score(inputs, "validation", destination / "validation_scores_directional")
    _, manifest = score(inputs, "test", destination / "test_scores", promotion_sha)
    from round_figures import render
    render(destination / "test_scores", destination, FIGURES)
    manifest["figures"] = {name: sha256(destination / name) for name in FIGURES}
    manifest["historical_promotion_sha256"] = promotion_sha
    write_json(destination / "test_scores/manifest.json", manifest)


def revise(inputs, promotion_path):
    legacy = inputs.root / "legacy_unsigned"
    if legacy.exists():
        raise FileExistsError("directional revision already published; use replay-revision")
    with tempfile.TemporaryDirectory(dir=inputs.root, prefix=".directional-revision-") as directory:
        path = Path(directory)
        revision_release(inputs, promotion_path, path)
        legacy.mkdir()
        for name in ("validation_scores", "test_scores", *FIGURES):
            (inputs.root / name).rename(legacy / name)
        for name in ("validation_scores_directional", "test_scores", *FIGURES):
            (path / name).rename(inputs.root / name)


def replay_revision(inputs, promotion_path):
    with tempfile.TemporaryDirectory(dir=inputs.root, prefix=".directional-replay-") as directory:
        path = Path(directory)
        revision_release(inputs, promotion_path, path)
        for name in ("validation_scores_directional", "test_scores"):
            require_same_release(inputs.root / name, path / name, (*TABLES, "manifest.json"))
        require_same_release(inputs.root, path, FIGURES)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("seal", "test", "replay", "revise", "replay-revision"))
    parser.add_argument("--round-root", type=Path, default=DEFAULT_ROUND)
    parser.add_argument("--work-root", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--plan", type=Path, default=RECORD / "plan.json")
    parser.add_argument("--promotion", "--output", type=Path, default=RECORD / "promotion.json")
    args = parser.parse_args()
    inputs = RoundInputs(args.round_root, args.work_root, args.corpus, args.plan)
    {"seal": seal, "test": test, "replay": replay,
     "revise": revise, "replay-revision": replay_revision}[args.action](inputs, args.promotion)
    print(f"Round 1 {args.action}: verified", flush=True)


if __name__ == "__main__":
    main()
