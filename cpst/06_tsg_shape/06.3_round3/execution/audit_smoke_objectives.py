"""Audit weighted objective components and seed behavior across the sealed CUDA smoke fits."""

import argparse
import csv
from itertools import groupby
import json
import math
from pathlib import Path

from composition_model import axis_run
from composition_roster import NO_X2_SHAPES
from matrix_composition import verify_matrix

RELATIVE_TOLERANCE = 1e-6
ABSOLUTE_TOLERANCE = 1e-4
EVALUATION_SEEDS = (22617, 100042)


def audit(matrix):
    paths, _ = verify_matrix(matrix, True)
    model_objectives = 0
    max_relative_error = 0.
    for (shape, fold), group in groupby(paths, key=lambda key: key[:2]):
        predictions, populations = set(), {}
        for coordinate in group:
            path = paths[coordinate]
            config = axis_run.read_json(path / "registration.json")
            if ((config["seed"], config["validation_seed"], config["test_seed"])
                    != (coordinate[2], *EVALUATION_SEEDS)):
                raise ValueError(f"invalid training or evaluation seed in {coordinate}")
            order = config["conditions"]
            if populations and populations != order:
                raise ValueError(f"training seed changed the registered conditions in {coordinate}")
            populations = order
            predictions.add(axis_run.digest(path / "prediction.npy"))

            with (path / "learning_curve.csv").open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            if [int(row["epoch"]) for row in rows] != [1, 2]:
                raise ValueError(f"expected two full epochs in {coordinate}")
            no_x2 = shape in NO_X2_SHAPES
            for row in rows:
                for phase in ("train", "validation"):
                    def component(name):
                        return float(row[f"{phase}_{name}"])

                    kls = [component(name) for name in ("control_kl", "treatment_kl", "prediction_kl")]
                    treatment = component("treatment_mse")
                    if any(term <= 0. for term in kls) or (treatment == 0.) != no_x2:
                        raise ValueError(f"objective is missing a registered term in {coordinate}, {phase}")
                    objective = (0.1 * math.fsum(kls) + component("control_mse") + treatment
                                 + (2 if no_x2 else 1) * component("delta_mse") + 1.75 * component("delta_l1"))
                    observed = component("loss")
                    if not math.isclose(observed, objective, rel_tol=RELATIVE_TOLERANCE,
                                        abs_tol=ABSOLUTE_TOLERANCE):
                        raise ValueError(f"weighted objective differs in {coordinate}, {phase}")
                    max_relative_error = max(max_relative_error, abs(observed - objective) / objective)
                    model_objectives += 1
        if len(predictions) != 2:
            raise ValueError(f"training seeds did not produce distinct predictions for {shape}, fold {fold}")
    return {
        "status": "PASS",
        "matrix_manifest_sha256": axis_run.digest(matrix / "manifest.json"),
        "fit_count": len(paths),
        "objective_rows": model_objectives,
        "objective_max_relative_error": max_relative_error,
        "objective_relative_tolerance": RELATIVE_TOLERANCE,
        "objective_absolute_tolerance": ABSOLUTE_TOLERANCE,
        "evaluation_seeds": {"validation": EVALUATION_SEEDS[0], "test": EVALUATION_SEEDS[1]},
        "conditions_invariant_and_predictions_distinct_across_training_seeds": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.matrix), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
