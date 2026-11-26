"""Gate each Step 06.4 KPGT reference before its reuse in the formal comparison."""

import argparse
import json
import os
from pathlib import Path

import pandas as pd
import torch

from chem_matrix import execute
from chem_report import METRICS, TENV_DEFAULT, AxisInputs, axis_evaluate, open_tenv_training, score_run
from chem_runs import (REFERENCE, ROUND_FOUR, STEP, AlphaAdapter, add_inputs, axis_run, bind, coordinates,
                       initialization, open_inputs)
from verify_alpha import (ACCEPTANCE, RELEASE, REFERENCE_ROSTER, verify_protocol, verify_reference,
                          verify_reference_release, verify_release)

SCORE_COLUMNS = ["value", "n_conditions", "n_eligible", "n_undefined"]
HISTORY_FIELDS = {"arguments", "nonce", "code_sha256", "schema", "protocol"}
ADDED_FIELDS = {"alpha_release_sha256", "representation"}
REUSE_SCHEMA = "carps-vae-alpha-chemical-reference-v1"
GATE = STEP / "log/363-gate.json"


def verify_bounded_gate():
    gate = axis_run.read_json(GATE)
    bounded, formal = coordinates(True), coordinates()
    if ((gate["status"], gate["ticket"], gate["independent_cuda_replays"]) != ("PASS", 363, len(bounded))
            or tuple(tuple(fit["coordinate"]) for fit in gate["fits"]) != bounded
            or tuple(map(tuple, gate["preflight"]["formal_roster"])) != formal
            or gate["metrics"]["independent_rescore"] != "exact"
            or gate["metrics"]["scores"] != len(bounded) * len(METRICS)
            or gate["figure_replay"]["independent_csv_only_renders"] != 2
            or gate["alpha"]["before_and_after_replay_sha256"] != axis_run.digest(ACCEPTANCE)
            or gate["alpha"]["release_sha256"] != axis_run.digest(RELEASE)):
        raise ValueError("formal Step 08 requires the accepted seven-arm bounded CUDA gate")


def verify_reference_scores(observed, recorded, registration_sha):
    reference = recorded.set_index("metric").reindex(METRICS)
    scored = observed.set_index("metric").reindex(METRICS)
    if (len(recorded) != len(METRICS) or len(observed) != len(METRICS)
            or not reference.run_sha256.eq(registration_sha).all()
            or not scored[SCORE_COLUMNS].equals(reference[SCORE_COLUMNS])):
        raise ValueError("rescored KPGT metrics, eligibility or registration do not exactly match Step 06.4")


def verify_live_alpha_dependencies(stage, release):
    """Check inherited loss modules against the code hashes in the frozen alpha release."""
    original = release["step_06_4"]["code_sha256"]
    paths = (
        "05_tsg_loss/05.1_loss_metrics/execution/objective.py",
        "05_tsg_loss/05.2_mse_l1/execution/additive_objective.py",
        "05_tsg_loss/05.2_mse_l1/execution/run_sweep.py",
    )
    for path in paths:
        if axis_run.digest(stage / path) != original[f"stage2/{path}"]:
            raise ValueError(f"inherited alpha runtime dependency differs from the pinned source: {path}")


def reference_plan(round_four, training, representations, source, bandwidth_training, *, accepted_plan=None):
    """Independently replay and score each valid reference, and return its verified fit and decision."""
    verify_bounded_gate()
    release = verify_release()
    verify_live_alpha_dependencies(STEP.parent, release)
    accepted = axis_run.read_json(ACCEPTANCE)
    release_sha = axis_run.digest(RELEASE)
    release_receipt = release["step_06_4"]["receipt_sha256"]
    if ((accepted["status"], accepted["release_sha256"], accepted["step_06_4_receipt_sha256"],
         accepted["validation_seed"], accepted["test_seed"])
            != ("PASS", release_sha, release_receipt, axis_run.VALIDATION_SEED, axis_run.TEST_SEED)
            or [(row["fold"], row["seed"]) for row in accepted["coordinates"]] != list(REFERENCE_ROSTER)):
        raise ValueError("published alpha parity receipt differs from the fixed release or reference roster")
    parity = {(row["fold"], row["seed"]): row for row in accepted["coordinates"]}
    defects = (OSError, ValueError, KeyError, RuntimeError)
    shared_defect = None
    try:
        ledger = verify_reference_release(round_four, release)
        scores = axis_evaluate._csv(round_four / "figure/data/fold-scores.csv")
        scores = scores.loc[scores["shape"].eq(release["default"]["architecture"])]
        bandwidths = axis_evaluate._csv(round_four / "log/bandwidths.csv").set_index("fold")
    except defects as error:
        shared_defect = type(error).__name__
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    torch.use_deterministic_algorithms(True)
    axis_run.kpgt.configure_backend("cuda")

    if accepted_plan is None:
        refits = {}
    else:
        records = accepted_plan["decisions"]
        if ([(row["coordinate"][1], row["coordinate"][2]) for row in records] != list(REFERENCE_ROSTER)
                or any(row["coordinate"][0] != REFERENCE or row["decision"] not in {"reuse", "refit"}
                       for row in records)):
            raise ValueError("invalid formal KPGT reference decisions")
        refits = {tuple(row["coordinate"]): row for row in records if row["decision"] == "refit"}
    decisions, reused, checked_bandwidths = [], {}, {}
    with open_inputs(training, representations) as (data, vectors), open_tenv_training(bandwidth_training) as raw:
        base = AxisInputs(raw)
        for fold, seed in REFERENCE_ROSTER:
            key = REFERENCE, fold, seed
            if key in refits:
                decisions.append(refits[key])
                continue
            row = {"coordinate": list(key)}
            stage = "sealed Step 06.4 ledgers"
            try:
                if shared_defect:
                    raise ValueError(shared_defect)
                stage = "sealed checkpoint lineage"
                path, config, replay, trained = verify_reference(ledger, release, fold, seed)
                stage = "current alpha and Step 08 protocol equivalence"
                verify_protocol(config, data, vectors, source)
                current = bind(key, data, vectors, source, False)
                if ({name: value for name, value in config.items() if name not in HISTORY_FIELDS}
                        != {name: value for name, value in current.items() if name not in HISTORY_FIELDS | ADDED_FIELDS}):
                    raise ValueError("registered KPGT model or inputs differ in Step 08")
                stage = "fresh exact alpha checkpoint replay"
                adapter = AlphaAdapter(source, current["model"]["transform"])
                model = adapter.build(seed, "cuda", n_genes=len(data.genes))
                initial = initialization(model, REFERENCE)
                registration_sha = replay["registration_sha256"]
                try:
                    values, test = axis_run.predict_test(path, adapter, model, data, config, registration_sha, trained)
                    result = axis_run.record_predictions(path, None, values, test, data, config, replay=True)
                finally:
                    del model
                predicted = {"fold": fold, "seed": seed, "registration_sha256": registration_sha,
                             "checkpoint_sha256": replay["outputs"]["checkpoint.pt"],
                             "validation_prediction_sha256": replay["outputs"]["validation_prediction.npy"],
                             "prediction_sha256": replay["outputs"]["prediction.npy"],
                             "validation_conditions": len(config["partitions"]["valid"]),
                             "test_replay": "exact", **result}
                if predicted != parity[fold, seed]:
                    raise ValueError("fresh replay differs from published Step 07 parity")
                stage = "independent exact Step 06.4 scoring"
                if fold not in checked_bandwidths:
                    bandwidth = axis_evaluate._bandwidth(base, fold)
                    if bandwidth != {"fold": fold, **bandwidths.loc[fold].to_dict()}:
                        raise ValueError("current training bandwidth differs from sealed Step 06.4")
                    checked_bandwidths[fold] = bandwidth["bandwidth"]
                _, metric_rows = score_run(path, config, replay, data, base, fold, False, checked_bandwidths[fold])
                verify_reference_scores(metric_rows, scores.loc[scores.fold.eq(fold) & scores.seed.eq(seed)],
                                        registration_sha)
                trace = axis_evaluate._csv(path / "learning_curve.csv")
                reused[key] = path, config, replay, trained, trace
            except defects as error:
                row.update(decision="refit", failed_gate=stage, error_type=type(error).__name__)
                print(json.dumps({**row, "details": str(error)}), flush=True)
            else:
                row.update(decision="reuse", registration_sha256=registration_sha,
                           checkpoint_sha256=predicted["checkpoint_sha256"],
                           prediction_sha256=predicted["prediction_sha256"],
                           validation_prediction_sha256=predicted["validation_prediction_sha256"],
                           initialization=initial)
                print(json.dumps({"coordinate": row["coordinate"], "decision": row["decision"]}), flush=True)
            decisions.append(row)
    return reused, {"schema": REUSE_SCHEMA, "alpha_release_sha256": release_sha,
                    "alpha_receipt_sha256": axis_run.digest(ACCEPTANCE), "step_06_4_receipt_sha256": release_receipt,
                    "decisions": decisions}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_inputs(parser)
    parser.add_argument("--round-four", type=Path, default=ROUND_FOUR)
    parser.add_argument("--bandwidth-training", type=Path, default=TENV_DEFAULT)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--logs", type=Path, required=True)
    parser.add_argument("--reference-only", action="store_true")
    args = parser.parse_args()
    reused, evidence = reference_plan(args.round_four, args.training, args.representations, args.source,
                                      args.bandwidth_training)
    if args.evidence.exists():
        if axis_run.read_json(args.evidence) != evidence:
            raise ValueError("new reference reuse decisions differ from the registered formal plan")
    else:
        with args.evidence.open("x", encoding="utf-8") as stream:
            json.dump(evidence, stream, sort_keys=True, indent=2)
            stream.write("\n")
    if args.reference_only:
        return
    fresh = tuple(key for key in coordinates() if key not in reused)
    result = execute(args.matrix, args.logs, args.training, args.representations, args.source, fresh=fresh)
    print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
