"""Score or exactly replay the matched TranSiGen axis evaluation report."""

import argparse
from hashlib import sha256 as hash_bytes
import json
import math
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parent
sys.path[:0] = [str(STEP / "output/api"), str(STAGE / "01_data_raw/output/api"),
                str(STAGE / "02_gene_axis/02.1_hvg-deg/output/api")]

from native_corpus import open_corpus
from prediction_axis import COMMON_SHA256, COMMON_SIZE, open_axes
from tsg_axis.disease_reference import disease_gene_sets, write_cm_directions
from tsg_axis.inputs import AxisInputs, score_fixed_prediction
from tsg_axis.metrics import (CONDITION, DISTRIBUTION, METRICS, condition_scores, distribution_scores,
                              fold_bandwidth, signed_disease_f1, summarize_conditions, summarize_folds)
from tsg_axis.protocol import MAX_EPOCHS, PATIENCE, PROTOCOL, REPORT_SCHEMA, RUN_SCHEMA
from training_data import TENV_AXIS_STEPS, TENV_BASELINE_AXIS, open_tenv_axis_training
from stage2_util.derivation.release import check_outputs, read_json as _read, sha256, write_json
from stage2_util.publication import publication_lock

PARTITIONS = {"train": 64, "valid": 24, "test": 24}
RUN_OUTPUTS = {"checkpoint.pt", "validation_prediction.npy", "learning_curve.csv", "training.json",
               "predictions.json", "prediction.npy", "predicted_delta.npy", "measured_delta.npy",
               "conditions.csv", "genes.csv"}
DISEASE_RELEASE = STAGE / "02_gene_axis/02.1_hvg-deg/output/release"


def _csv(path):
    return pd.read_csv(path, keep_default_na=False, float_precision="round_trip")


def _file_hashes(root):
    return {path.relative_to(root).as_posix(): sha256(path)
            for path in sorted(root.rglob("*")) if path.is_file()}


def _run_inputs(paths, smoke):
    runs = {}
    for path in paths:
        config_path, manifest_path = (path / name for name in ("registration.json", "manifest.json"))
        config, manifest = _read(config_path), _read(manifest_path)
        axis, fold = manifest["axis"], manifest["fold"]
        args = config["arguments"]
        if (config["schema"] != RUN_SCHEMA or manifest["schema"] != RUN_SCHEMA
                or config.get("protocol") != PROTOCOL or manifest.get("protocol") != PROTOCOL
                or config["patience"] != PATIENCE or config["max_epochs"] != (2 if smoke else MAX_EPOCHS)
                or axis not in TENV_AXIS_STEPS or type(fold) is not int or fold not in range(3)
                or (axis, fold, smoke) != (args["axis"], args["fold"], args["smoke"])
                or manifest["smoke"] is not smoke or manifest["seed"] != config["seed"]
                or config["seed"] != 42 or config["scenario"] != "cold-drug"
                or manifest["validation_replay"] != "exact" or manifest["test_replay"] != "exact"
                or manifest.get("validation_objective_replay") != "exact"
                or manifest["registration_sha256"] != sha256(config_path)
                or set(manifest["outputs"]) != RUN_OUTPUTS):
            raise ValueError(f"not an exactly replayed run for this report: {path}")
        check_outputs(path, manifest)
        if (axis, fold) in runs:
            raise ValueError(f"duplicate run coordinate: {axis} fold {fold}")
        runs[axis, fold] = (path, config, manifest)
    if smoke:
        if len(runs) != 1 or any(axis != TENV_BASELINE_AXIS for axis, _ in runs):
            raise ValueError("bounded reporting requires one HVG-DEG 5k smoke run")
    elif set(runs) != {(axis, fold) for axis in TENV_AXIS_STEPS for fold in range(3)}:
        raise ValueError("formal reporting requires 18 sealed fits on the six axes and three folds")
    return runs


def _verify_inputs(path, config, data, base, axis, fold, smoke):
    family, width = axis.rsplit("-", 1)
    axes = open_axes(family=family)
    released_axis = axes.axis(int(width.removesuffix("k")) * 1000)
    stage03 = data.release.manifest["inputs"]
    baseline = base.release.manifest["inputs"]
    if (config["axis_sha256"] != stage03["axis_sha256"]
            or config["axis_sha256"] != released_axis.sha256
            or config["axis_family_manifest_sha256"] != axes.manifest_sha256
            or config["training_sha256"] != data.release.manifest_sha256
            or config["raw_sha256"] != stage03["raw_sha256"]
            or config["raw_sha256"] != baseline["raw_sha256"]
            or config["folds"] != stage03["folds"] or config["folds"] != baseline["folds"]
            or not _csv(path / "genes.csv").equals(data.genes)):
        raise ValueError(f"axis, training, raw or frozen fold identity differs: {path}")
    for phase, limit in PARTITIONS.items():
        rows = data.rows(fold, phase)
        reference_rows = base.rows(fold, phase)
        if smoke:
            rows, reference_rows = rows[:limit], reference_rows[:limit]
        ids = data.condition_ids(rows)
        if (rows.tolist() != config["partitions"][phase] or list(ids) != config["conditions"][phase]
                or ids != base.condition_ids(reference_rows)):
            raise ValueError(f"registered {phase} condition IDs differ from the shared fold: {path}")
    test = np.asarray(config["partitions"]["test"], dtype=np.int64)
    reference_test = base.rows(fold, "test")
    return test, reference_test[:len(test)]


def _bandwidth(base, fold):
    data = base.release
    rows = base.rows(fold, "train")
    parts = []
    for start in range(0, len(data.treatments), 64):
        stop = min(start + 64, len(data.treatments))
        indexes = rows[(rows >= start) & (rows < stop)]
        if len(indexes):
            parts.append(data.read_treatments(start, stop, "delta")[indexes - start])
    measured = np.concatenate(parts)
    record = {"fold": fold, "bandwidth": fold_bandwidth(measured), "rows": len(measured),
              "width": len(data.genes), "pair_count": math.comb(len(measured), 2),
              "training_manifest_sha256": data.manifest_sha256,
              "training_delta_values_sha256": hash_bytes(measured.tobytes()).hexdigest(),
              "training_conditions_sha256": hash_bytes(base.condition_table(rows).to_csv(index=False).encode()).hexdigest()}
    return record


def _score_run(path, config, manifest, data, base, axis, fold, smoke, directions, disease, bandwidth):
    rows, reference_rows = _verify_inputs(path, config, data, base, axis, fold, smoke)
    conditions = _csv(path / "conditions.csv")
    if (not conditions.equals(data.condition_table(rows).reset_index(drop=True))
            or len(conditions) != manifest["test_conditions"]):
        raise ValueError(f"fixed prediction condition table differs: {path}")
    prediction = np.load(path / "prediction.npy", mmap_mode="r", allow_pickle=False)
    target = score_fixed_prediction(data, fold, prediction, conditions.condition_id, rows=rows)
    for name, expected in (("predicted_delta", target.predicted_delta),
                           ("measured_delta", target.measured_delta)):
        array = np.load(path / f"{name}.npy", mmap_mode="r", allow_pickle=False)
        if array.dtype != expected.dtype or not np.array_equal(array, expected):
            raise ValueError(f"fixed prediction {name} differs from paired source arrays: {path}")

    genes = data.genes.gene_index.to_numpy()
    paired, disease_scores = [], []
    for condition, cm_row, truth, predicted in zip(conditions.condition_id, reference_rows,
                                                   target.measured_delta, target.predicted_delta, strict=True):
        key = {"axis": axis, "fold": fold, "condition_id": condition}
        paired.append({**key, **condition_scores(truth, predicted, genes)})
        for name, members in disease.items():
            result = signed_disease_f1(predicted, genes, directions[cm_row], members)
            disease_scores.append({**key, "disease": name, "f1": result.f1,
                                   "tp": result.tp, "fp": result.fp, "fn": result.fn})
    paired = pd.DataFrame(paired)
    disease_scores = pd.DataFrame(disease_scores)
    dist = dict(zip(DISTRIBUTION, distribution_scores(target.measured_delta, target.predicted_delta,
                                                       bandwidth), strict=True))
    fold_metrics = []
    for metric in METRICS:
        score = summarize_conditions(paired[metric]) if metric in CONDITION else {
            "value": dist[metric], "n_conditions": len(paired),
            "n_eligible": len(paired), "n_undefined": 0,
        }
        fold_metrics.append({"axis": axis, "fold": fold, "metric": metric, **score})
    fold_disease = [{"axis": axis, "fold": fold, "disease": name,
                     **summarize_conditions(disease_scores.loc[disease_scores.disease.eq(name), "f1"])}
                    for name in disease]
    return paired, disease_scores, fold_metrics, fold_disease


def _build_report(root, runs, step03, raw, smoke):
    root.mkdir(exist_ok=True)
    reference = root / "reference"
    reference.mkdir()
    figures = root / "figure/data"
    figures.mkdir(parents=True)
    disease, evidence = disease_gene_sets(DISEASE_RELEASE)
    pd.DataFrame([(name, rank, gene) for name, genes in disease.items()
                  for rank, gene in enumerate(genes, start=1)],
                 columns=["disease", "disease_rank", "gene_index"]).to_csv(reference / "deg1000.csv", index=False)
    pairs, conditions_disease, metric_folds, disease_folds, provenance = [], [], [], [], []
    with open_tenv_axis_training(TENV_BASELINE_AXIS, step03) as training, open_corpus(raw) as corpus:
        if (corpus.manifest["common_axis"] != {"size": COMMON_SIZE, "sha256": COMMON_SHA256}
                or corpus.manifest_sha256 != training.manifest["inputs"]["raw_sha256"]
                or corpus.manifest["sources"]["tenv"]["sha256"] != training.manifest["inputs"]["raw_counts_sha256"]):
            raise ValueError("shared CM truth must use the native tenv counts of the released training data")
        base = AxisInputs(training)
        base.condition_table(np.arange(len(training.treatments))).to_csv(reference / "conditions.csv", index=False)
        direction_path = reference / "cm-measured-direction.npy"
        write_cm_directions(corpus.source("tenv"), training, direction_path)
        directions = np.load(direction_path, mmap_mode="r", allow_pickle=False)
        folds = sorted({fold for _, fold in runs})
        bands = {fold: _bandwidth(base, fold) for fold in folds}
        pd.DataFrame(bands.values()).to_csv(root / "bandwidths.csv", index=False)
        membership = []
        for fold in folds:
            for phase, size in PARTITIONS.items():
                rows = base.rows(fold, phase)
                table = base.condition_table(rows)
                included = np.arange(len(rows)) < size if smoke else True
                membership.append(table.assign(fold=fold, partition=phase, included_in_run=included))
        pd.concat(membership, ignore_index=True).to_csv(root / "fold-memberships.csv", index=False)
        for axis in TENV_AXIS_STEPS:
            for fold in folds:
                if (axis, fold) not in runs:
                    continue
                path, config, manifest = runs[axis, fold]
                with open_tenv_axis_training(axis, step03) as data:
                    scores = _score_run(path, config, manifest, AxisInputs(data), base, axis, fold,
                                        smoke, directions, disease, bands[fold]["bandwidth"])
                paired, by_disease, fold_metric, fold_disease = scores
                pairs.append(paired)
                conditions_disease.append(by_disease)
                metric_folds.extend(fold_metric)
                disease_folds.extend(fold_disease)
                provenance.append({"axis": axis, "fold": fold, "axis_sha256": config["axis_sha256"],
                                   "registration_sha256": manifest["registration_sha256"],
                                   "run_manifest_sha256": sha256(path / "manifest.json"),
                                   "training_manifest_sha256": config["training_sha256"]})

        common = {"raw_manifest_sha256": corpus.manifest_sha256,
                  "raw_tenv_counts_sha256": training.manifest["inputs"]["raw_counts_sha256"],
                  "common_axis_sha256": COMMON_SHA256, "common_axis_size": COMMON_SIZE,
                  "reference_training_manifest_sha256": training.manifest_sha256,
                  "disease_evidence_sha256": evidence,
                  "disease_axis_release_manifest_sha256": sha256(DISEASE_RELEASE / "manifest.json")}
    pd.concat(pairs, ignore_index=True).to_csv(root / "condition-metrics.csv", index=False)
    pd.concat(conditions_disease, ignore_index=True).to_csv(root / "condition-disease-f1.csv", index=False)
    metric_folds, disease_folds = pd.DataFrame(metric_folds), pd.DataFrame(disease_folds)
    metric_folds.to_csv(root / "fold-metrics.csv", index=False)
    disease_folds.to_csv(root / "fold-disease-f1.csv", index=False)
    summarize_folds(metric_folds, ["axis", "metric"]).to_csv(figures / "delta-metrics.csv", index=False)
    summarize_folds(disease_folds, ["axis", "disease"]).to_csv(figures / "disease-f1.csv", index=False)
    write_json(root / "manifest.json", {
        "schema": REPORT_SCHEMA, "protocol": PROTOCOL, "smoke": smoke, "reference": common, "runs": provenance,
        "metrics": METRICS, "kernel": "exp(-scaled_euclidean² / (2 * bandwidth²))",
        "bandwidth": "median distinct-row-pair scaled Euclidean distance on full HVG-DEG 5k training measured Delta, per fold",
        "aggregation": "equal defined conditions per fold, equal defined folds per axis",
        "interpretation": "retrospective comparison of published fixed axes conditional on measured matched plate DMSO",
        "outputs": _file_hashes(root),
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("score", "verify"))
    parser.add_argument("--run", action="append", type=Path, required=True, help="exactly replayed run, repeat per coordinate")
    parser.add_argument("--step03", type=Path, default=STAGE / "03_data_train")
    parser.add_argument("--raw", type=Path, default=STAGE / "01_data_raw/output")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true", help="score one bounded HVG-DEG 5k run")
    args = parser.parse_args()
    runs = _run_inputs(args.run, args.smoke)
    with publication_lock(args.output) as output:
        if args.command == "score" and output.exists():
            raise ValueError("score requires an unused output directory")
        with TemporaryDirectory(prefix=f".{output.name}-", dir=output.parent) as scratch:
            candidate = Path(scratch) / "report"
            _build_report(candidate, runs, args.step03, args.raw, args.smoke)
            if args.command == "score":
                candidate.rename(output)
            else:
                manifest = _read(candidate / "manifest.json")
                for name in ("manifest.json", *manifest["outputs"]):
                    if sha256(output / name) != sha256(candidate / name):
                        raise ValueError(f"evaluation replay differs: {name}")
    print(json.dumps({"result": "PASS", "manifest_sha256": sha256(output / "manifest.json")}, sort_keys=True))


if __name__ == "__main__":
    main()
