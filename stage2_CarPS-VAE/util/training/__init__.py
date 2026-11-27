"""Training identity for the registered matched-control TranSiGen studies."""


def experiment_configuration(arguments, bound, *, dose_conditioning="latent-numeric"):
    """Describe the model and frozen condition population in a run registration."""
    model = bound["model"]
    objective = bound.get("objective")
    if objective is None:
        objective = {"form": "native", "weights": {
            "mse_x1": 1, "mse_x2": 1, "mse_delta": 1,
            "kl_x1": model["beta"], "kl_x2": model["beta"], "kl_prediction": model["beta"],
        }}
    return {
        "architecture": bound.get("shape") or model.get("architecture") or model["method_identity"],
        "objective": objective,
        "representation": model["molecule_feature"],
        "dose": {"conditioning": model.get("method", dose_conditioning), "transform": model["transform"]},
        "axis": {"name": arguments["axis"], "sha256": bound["axis_sha256"]},
        "split": bound["scenario"],
        "fold": arguments["fold"],
        "seed": bound["seed"],
    }
