"""Assemble and validate the bounded Step 09 CUDA acceptance evidence."""

import argparse
import json
from pathlib import Path

from axis_matrix import verify_matrix
from axis_model import AxisAdapter, named_tensor_digest, paired_initialization
from axis_protocol import REFERENCE_AXIS, coordinates
from axis_runs import axis_run, configure_cuda
from axis_scores import REPORT_METRICS
from stage2_util.comparison import FIGURES
from audit_chem_resources import audit_resources


def process_resources(path):
    lines = path.read_text().splitlines()
    maximum = next(line for line in lines if "Maximum resident set size (kbytes)" in line)
    elapsed = next(line for line in lines if "Elapsed (wall clock) time" in line)
    if not any(line.strip() == "Exit status: 0" for line in lines):
        raise ValueError(f"unsuccessful process: {path.name}")
    return {"peak_rss_kib": int(maximum.rsplit(": ", 1)[1]),
            "elapsed": elapsed.rsplit(": ", 1)[1]}


def audit_initialized_projections(fresh, source):
    """Rebuild the bounded models to compare the full state of equal-width axes."""
    configure_cuda()
    signatures = {}
    for (axis, _, seed), (_, config, _, trained, _) in fresh.items():
        model_configuration = config["model"]
        width = model_configuration["n_genes"]
        adapter = AxisAdapter(source, model_configuration["transform"])
        model = adapter.build(seed, "cuda", n_genes=width)
        if paired_initialization(model) != trained["initialization"]:
            raise ValueError(f"fresh GPU construction differs from the bounded train initialization: {axis}")
        weights = model.state_dict()
        signature = named_tensor_digest(weights, weights.keys())
        if signatures.setdefault(width, signature) != signature:
            raise ValueError(f"equal-width axis initialization differs at {width} genes")
    return signatures


def gate_receipt(work, coordinator, source):
    matrix = work / "matrix"
    report = work / "report"
    fresh = verify_matrix(matrix, True)
    alpha = axis_run.read_json(work / "alpha.json")
    reference = axis_run.read_json(work / "reference.json")
    figures = axis_run.read_json(work / "figures.json")
    score = axis_run.read_json(report / "manifest.json")
    expected = [(fold, seed) for axis, fold, seed in coordinates() if axis == REFERENCE_AXIS]
    if alpha["status"] != "PASS" or [(row["fold"], row["seed"]) for row in alpha["coordinates"]] != expected:
        raise ValueError("Step 07 current exact replay does not cover all nine sealed alpha fits")
    decisions = reference["decisions"]
    if ([(row["coordinate"][1], row["coordinate"][2]) for row in decisions] != expected
            or any(row["coordinate"][0] != REFERENCE_AXIS or row["decision"] != "reuse"
                   or row["validation_replay"] != "exact" or row["test_replay"] != "exact" for row in decisions)):
        raise ValueError("Step 08 current exact replay does not cover all nine sealed UniMol fits")
    trained = fresh[REFERENCE_AXIS, 0, 17][3]
    if trained["initialization"] != decisions[0]["paired_initialization"]:
        raise ValueError("Step 09 bounded reference and sealed UniMol initialization differ")
    if (figures["status"] != "PASS" or figures["independent_csv_only_renders"] != 2
            or set(figures["png_sha256"]) != set(FIGURES)
            or score["scores"] != len(fresh) * len(REPORT_METRICS)):
        raise ValueError("gate did not score the full prediction axes and render the three canonical figures")
    for subdirectory, artifacts in (("figures", figures["png_sha256"]),
                                    ("figures/data", figures["csv_sha256"]), ("report", score["outputs"])):
        if axis_run.file_hashes(work / subdirectory, artifacts) != artifacts:
            raise ValueError(f"gate artifacts changed: {subdirectory}")
    coverage, resources = audit_resources(work / "matrix-audit", coordinator, coordinates(True), bracket=True)
    evidence = {
        "status": "PASS", "ticket": 367, "formal_coordinates": len(coordinates()),
        "bounded_coordinates": len(fresh), "fresh_alpha_replays": len(alpha["coordinates"]),
        "fresh_sealed_unimol_replays": len(decisions), "step08_acceptance_sha256": reference["step08_acceptance_sha256"],
        "matrix_sha256": axis_run.digest(matrix / "manifest.json"),
        "preflight_sha256": axis_run.digest(matrix / "preflight.json"),
        "reference_plan_sha256": axis_run.digest(work / "reference.json"),
        "alpha_replay_sha256": axis_run.digest(work / "alpha.json"),
        "report_sha256": axis_run.digest(report / "manifest.json"),
        "figures": figures, "fit_budgets": resources,
        "process_resources": {phase: process_resources(work / f"{phase}.time")
                              for phase in ("alpha", "gate", "reference", "score", "figures")},
        "equal_width_initialized_state_sha256": audit_initialized_projections(fresh, source),
        "bounded_fit_lineage": axis_run.read_json(matrix / "manifest.json")["fits"],
    }
    return evidence, coverage


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--coordinator", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    args = parser.parse_args()
    evidence, samples = gate_receipt(args.work, args.coordinator, args.source)
    args.output.mkdir()
    samples_path = args.output / "367-coordinator-resources.jsonl"
    samples_path.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in samples))
    evidence["coordinator_samples_sha256"] = axis_run.digest(samples_path)
    axis_run.write_json(args.output / "367-bounded.json", evidence)
    print(json.dumps({"status": evidence["status"], "bounded_coordinates": evidence["bounded_coordinates"]}))


if __name__ == "__main__":
    main()
