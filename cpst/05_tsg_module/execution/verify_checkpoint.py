"""Independently calculate each bounded fit's selected validation objective."""

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np
import torch

from run_module import axis_run, build_model
from matrix_measured import verify_matrix
from measured_roster import CONTROL_ONLY, MEASURED_DELTA


def prior(mean, logvar):
    mean, logvar = mean.double(), logvar.double()
    return (mean.square() + logvar.expm1() - logvar).sum().item() / 2


def alignment(pred_mu, pred_logvar, ref_mu, ref_logvar):
    pred_mu, ref_mu = pred_mu.double(), ref_mu.double()
    log_ratio = pred_logvar.double() - ref_logvar.double()
    scaled_distance = (pred_mu - ref_mu).square() / ref_logvar.double().exp()
    return (scaled_distance + log_ratio.expm1() - log_ratio).sum().item() / 2


def squared_distance(left, right):
    return (left.double() - right.double()).square().sum().item()


def verify_fit(run, data, config, args):
    device = torch.device(args["device"])
    native, model = build_model(args, config, data)
    state = torch.load(run / "checkpoint.pt", map_location="cpu", weights_only=True)
    model.load_state_dict(state["state_dict"])
    model.eval()
    shape = config["shape"]
    validation = np.asarray(config["partitions"]["valid"])
    inputs = axis_run.DoseView(data.model_input("valid"), native.transform)
    batch_size = config["model"]["batch_size"]
    contributions, predictions = [], []
    with axis_run.fixed_rng(config["validation_seed"], device), torch.no_grad():
        for start in range(0, len(validation), batch_size):
            rows = validation[start:start + batch_size]
            x1, x2, features = (
                torch.as_tensor(values, device=device)
                for values in (inputs.control[rows], inputs.expression[rows], inputs.fingerprints(rows))
            )
            x1_rec, mu1, logvar1, pred, pred_mu, pred_logvar, _ = model(x1, features)
            z2, mu2, logvar2 = model.encode_x2(x2)
            x2_rec = model.decode_x2(z2)
            control = squared_distance(x1, x1_rec)
            prediction = squared_distance(x2, pred)
            measured, measured_prior, prediction_alignment = 0.0, 0.0, 0.0
            if shape != CONTROL_ONLY:
                if shape in MEASURED_DELTA:
                    measured_delta = x2 - x1
                    delta_rec, mu2, logvar2 = model.reconstruct_measured_delta(measured_delta)
                    measured = squared_distance(measured_delta, delta_rec)
                else:
                    measured = squared_distance(x2, x2_rec)
                measured_prior = prior(mu2, logvar2)
                prediction_alignment = alignment(pred_mu, pred_logvar, mu2, logvar2)
            contributions.append((control, measured, prediction, prior(mu1, logvar1),
                                  measured_prior, prediction_alignment))
            predictions.append(pred.cpu().numpy())
    if not np.array_equal(np.concatenate(predictions), np.load(run / "validation_prediction.npy", allow_pickle=False)):
        raise ValueError(f"Independent checkpoint prediction differs: {run}")
    totals = [math.fsum(column) / len(validation) for column in zip(*contributions, strict=True)]
    observed = math.fsum(totals[:3]) + 0.1 * math.fsum(totals[3:])
    expected = state["validation_loss"]
    if not math.isclose(observed, expected, rel_tol=1e-6, abs_tol=1e-4):
        raise ValueError(f"Independent full validation loss differs: {run}: {observed} vs {expected}")
    return {"shape": shape, "selected_epoch": state["epoch"], "independent_loss": observed,
            "checkpoint_loss": expected, "components_per_condition": dict(zip(
                ("control_mse", "measured_delta_mse", "predicted_delta_mse", "control_prior",
                 "measured_prior", "predicted_to_measured_kl"), totals, strict=True))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    args = parser.parse_args()
    paths, _ = verify_matrix(args.matrix, True)
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    records = []
    for path in paths.values():
        config = axis_run.read_json(path / "registration.json")
        registered = config["arguments"]
        torch.use_deterministic_algorithms(True)
        axis_run.kpgt.configure_backend(registered["device"])
        with axis_run.open_inputs(registered["training"], registered["representations"]) as (data, _):
            data.bind_protocol(axis_run.shared_protocol())
            records.append(verify_fit(path, data, config, registered))
    print(json.dumps({"status": "PASS", "selected_objectives": records}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
