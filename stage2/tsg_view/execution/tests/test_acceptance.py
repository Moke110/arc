"""Acceptance receipt and Steps 01-06 invariance seams of ``stage2/tsg_view``.

Every test builds a synthetic study tree from scratch (report/figures/health included),
so the suite never touches the accepted multi-gigabyte products, a GPU, or the network.
"""
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXECUTION = HERE.parent
REPO = HERE.parents[3]
if str(EXECUTION) not in sys.path:
    sys.path.insert(0, str(EXECUTION))

import acceptance  # noqa: E402

SOURCES = ("tenv", "tahoe_bulk")
VIEWS = {
    "tenv": ("full", "merge", "solo", "solo-deepest", "full-pooled-label"),
    "tahoe_bulk": ("full", "merge", "solo", "full-pooled-label"),
}
ARMS = {
    "tenv": ("full", "merge", "solo", "full-pooled-label", "solo-deepest"),
    "tahoe_bulk": ("full", "merge", "solo"),
}
STEPS = ("01_bases", "02_models_tenv-5k", "03_tsg_drug_emb", "04_tsg_dose",
         "05_tsg_width_dose", "06_tsg_transfer")
STAMP = "2026-11-22T12:00:00Z"
ENVIRONMENT = {"python": "3.13.13", "torch": "2.9.0", "cuda": "13.0", "device": "cuda",
               "driver": "580.95", "numpy": "2.4.6", "matplotlib": "3.10.9"}

#: The verify artifact's transparency block: the dedup refactor moved two code digests after
#: the plan was sealed, which is allowed but must reach the receipt unchanged.
CODE_DRIFT = [
    {"file": "formal.py", "section": "code",
     "sealed_sha256": "1" * 64, "current_sha256": "2" * 64},
    {"file": "learner.py", "section": "learner_runtime",
     "sealed_sha256": "3" * 64, "current_sha256": "4" * 64},
]
SENSITIVITY_NOTE = ("the sealed anchor was measured at seed 42; this seed-43 arm is "
                    "reported for transparency and claims no anchor reproduction")
CODE_IDENTITY = {
    "status": "CODE_DRIFT", "code_drift": CODE_DRIFT, "driver_code_drift": CODE_DRIFT[:1],
    "sealed": {"formal.py": "1" * 64}, "current": {"formal.py": "2" * 64},
    "note": ("verify records every code digest that moved since the plan was sealed; the "
             "refactor-only drift is reported here and never hidden"),
}


def expected_origins():
    """The 45-fit schedule, spelled out independently of ``acceptance.plan_origins``."""
    rows = []
    for source in SOURCES:
        for view in ARMS[source]:
            seeds = (42,) if view == "solo-deepest" else (42, 43)
            for fold in (0, 1, 2):
                for seed in seeds:
                    rows.append({"source": source, "view": view, "fold": fold,
                                 "seed": seed, "method": "no-dose",
                                 "scenario": "cold-drug-dose"})
    return rows


ORIGINS = expected_origins()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, text=""):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(text, bytes):
        path.write_bytes(text)
    else:
        path.write_text(text)
    return path


def write_json(path, value):
    return write(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def record(path):
    path = Path(path)
    return {"file": path.name, "sha256": sha(path), "bytes": path.stat().st_size}


def csv(path, header, rows):
    text = ",".join(header) + "\n" + "".join(",".join(map(str, row)) + "\n" for row in rows)
    return write(path, text)


RENDERER = '''#!/usr/bin/env python3
"""Deterministic fixture renderer: PNG bytes come from the CSVs alone."""
import argparse, hashlib, json, pathlib

SUITES = {
    "tenv": ("learning_curves_tenv", "metrics_tenv", "metrics_per_dose_tenv"),
    "tahoe": ("learning_curves_tahoe_bulk", "metrics_tahoe_bulk",
              "metrics_per_dose_tahoe_bulk"),
    "cross": ("views_effect",),
}
SUITES["all"] = SUITES["tenv"] + SUITES["tahoe"] + SUITES["cross"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--suite", default="all")
    args = parser.parse_args()
    suite = "tahoe_bulk" if args.suite == "tahoe" else args.suite
    names = sorted(path.name for path in args.data.glob("*.csv"))
    digest = hashlib.sha256()
    inputs = {}
    for name in names:
        raw = (args.data / name).read_bytes()
        inputs[name] = hashlib.sha256(raw).hexdigest()
        digest.update(name.encode() + b"\\0" + raw)
    args.output.mkdir(parents=True, exist_ok=True)
    figures = {}
    for fig in SUITES[suite]:
        body = digest.hexdigest().encode() + b":" + fig.encode()
        target = args.output / (fig + ".png")
        target.write_bytes(b"PNG" + body)
        figures[fig + ".png"] = hashlib.sha256(target.read_bytes()).hexdigest()
    status = {
        "schema": "tsg-view-render-v1",
        "suite": suite,
        "dpi": 150,
        "matplotlib": "3.10.9",
        "metadata_software": None,
        "renderer_sha256": hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
        "figures": figures,
        "csv_inputs": inputs,
        "pending": {} if suite == "all" else {"tahoe_bulk": "pending"},
    }
    (args.output / "render_status.json").write_text(
        json.dumps(status, indent=2, sort_keys=True) + "\\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
'''


class Fixture:
    """A complete, digest-consistent ``tsg_view`` tree under ``root``."""

    def __init__(self, root, *, complete=True):
        self.root = Path(root)
        self.complete = complete

    # -- primitives ------------------------------------------------------- #
    def step(self, name):
        return self.root / "stage2" / "tsg_view" / name

    def build(self, *, git=True, baseline_mode="manifest", health=True, verify_all=True,
              report=True, figures=True, formal_verify=True):
        self._inputs()
        self._steps()
        self._smoke()
        self._views()
        self._freezes()
        self._truth()
        self._runs()
        if report:
            self._report()
        if figures:
            self._figures()
        # The commit first: health transcripts record the commit they were taken at, and a
        # sealed receipt refuses a transcript set that predates HEAD.
        if git:
            self._git()
        else:
            write(self.root / ".git-placeholder", "")
        if health:
            self._health(verify_all=verify_all)
        self._baseline(baseline_mode)
        return self.root

    def head(self):
        if not (self.root / ".git").exists():
            return None
        completed = subprocess.run(["git", "-C", str(self.root), "rev-parse", "HEAD"],
                                   check=True, capture_output=True, text=True)
        return completed.stdout.strip()

    # -- accepted inputs -------------------------------------------------- #
    def _inputs(self):
        write_json(self.root / "stage2/01_bases/output/view-5k/manifest.json",
                   {"schema": "stage2-bases-view-5k-release-v1", "size": 5000})
        write(self.root / "stage2/01_bases/output/gene-axis/gene_axis_5000.csv",
              ",".join(f"gene_{i}" for i in range(10)) + "\n")

    def _steps(self):
        for index, name in enumerate(STEPS):
            write(self.root / "stage2" / name / "README.md", f"# {name}\n")
            write(self.root / "stage2" / name / "AGENTS.md", f"routing {name}\n")
            write_json(self.root / "stage2" / name / "log" / "acceptance.json",
                       {"schema": f"{name}-acceptance-v1", "status": "passed",
                        "index": index})

    def _product_sha(self, source, view):
        return sha(self.step("output/views") / view / f"{source}.json")

    def _freeze_sha(self, source, view):
        return sha(self.step("execution/freezes") / view / source / "manifest.json")

    def _truth_shas(self, source, fold):
        return {"primary": sha(self.step("output/truth") / source / "primary"
                               / f"fold-{fold}" / "manifest.json"),
                "condition_mean": sha(self.step("output/truth") / source
                                      / "condition_mean" / f"fold-{fold}"
                                      / "manifest.json")}

    def _smoke(self):
        for index, source in enumerate(SOURCES, start=1):
            write_json(self.step(f"execution/checkpoints/merger-0{index}")
                       / f"smoke-full-{source}.json", {
                "schema": "tsg-view-smoke-v1", "status": "complete",
                "source": source, "view": "full", "fold": 0, "seed": 42,
                "device": "cuda", "epochs": 2, "rows": 64,
                "epochs_run": 2, "best_epoch": 2,
                "replay": "train -> replay -> test -> replay, one process per phase",
                "exact_replay": "validation and test predictions replayed bit-exactly in "
                                "independent fresh processes (phases 2 and 4)",
                "phases": [{"phase": "train", "returncode": 0},
                           {"phase": "replay", "returncode": 0},
                           {"phase": "test", "returncode": 0},
                           {"phase": "replay", "returncode": 0}]})

    # -- view products ---------------------------------------------------- #
    def _views(self):
        gene = record(self.root / "stage2/01_bases/output/gene-axis/gene_axis_5000.csv")
        release = record(self.root / "stage2/01_bases/output/view-5k/manifest.json")
        for source in SOURCES:
            for view in VIEWS[source]:
                folder = self.step("output/views") / view
                folder.mkdir(parents=True, exist_ok=True)
                profile_ids = write(folder / f"{source}.profile_ids.csv",
                                    "profile_id\n" + "".join(f"{source}-{i}\n" for i in range(4)))
                members = control = None
                if view in ("merge", "solo", "solo-deepest", "full-pooled-label"):
                    members = write(folder / f"{source}.members.csv",
                                    "profile_id,member\nx,y\n")
                    control = write(folder / f"{source}.control_sources.csv",
                                    "profile_id,control\nx,y\n")
                payload = alias = None
                if view == "full":
                    alias = {
                        "release": "stage2/01_bases/output/view-5k",
                        "manifest_sha256": release["sha256"],
                        "payload_sha256": {"tenv": "a" * 64, "tahoe_bulk": "b" * 64},
                        "outputs": {"tenv.h5ad": "c" * 64, "tahoe_bulk.h5ad": "d" * 64},
                        "control": {"rule": "mean of matched vehicles"},
                        "normalization": {"formula": "log1p"},
                        "recipe": "ADR 0077 unchanged",
                    }
                elif view == "full-pooled-label":
                    delta = write(folder / f"{source}.delta.npy", f"delta-{source}")
                    payload = {"kind": "delta", **record(delta), "rows": 4, "cols": 5000,
                               "dtype": "float32"}
                else:
                    payload_file = write(folder / f"{source}.h5ad", f"h5ad-{source}-{view}")
                    payload = {**record(payload_file), "rows": 4, "cols": 5000,
                               "dtype": "float32", "layers": ["control", "delta"]}
                manifest = {
                    "schema": "tsg-view-product-v1",
                    "product": f"{source}-{view}-5k-view",
                    "view": view, "source": source, "rows": 4, "conditions": 2,
                    "treated": 4,
                    "axis": {"file": "stage2/01_bases/output/gene-axis/gene_axis_5000.csv",
                             "sha256": gene["sha256"], "size": 5000},
                    "normalization": {"formula": "log1p(1e4 * count / total)"},
                    "recipe": f"{view} recipe", "control_rule": f"{view} control rule",
                    "seed": 20261122 if view == "solo" else None,
                    "selection": "largest retained-axis total"
                                 if view == "solo-deepest" else None,
                    "profile_ids": record(profile_ids),
                    "members": record(members) if members else None,
                    "control_sources": record(control) if control else None,
                    "retained_totals": {"min": 1, "median": 2.0, "max": 3, "sum": 10},
                    "source_manifests": {
                        "accepted_view": "stage2/01_bases/output/view-5k/manifest.json",
                        "accepted_view_sha256": release["sha256"],
                        "gene_axis": "stage2/01_bases/output/gene-axis/gene_axis_5000.csv",
                        "gene_axis_sha256": gene["sha256"],
                    },
                    "payload": payload, "alias": alias,
                }
                write_json(folder / f"{source}.json", manifest)

    # -- freezes ---------------------------------------------------------- #
    def _sealed(self, source):
        folder = (self.root
                  / "stage2/04_tsg_dose/execution/checkpoints/round1-axis-local-inputs"
                  / source)
        folder.mkdir(parents=True, exist_ok=True)
        population = csv(folder / "population.csv", ["profile_id"], [["p"]])
        splits = csv(folder / "splits.csv", ["scenario", "fold"], [["s", 0]])
        parents = csv(folder / "parents.csv", ["parent"], [["p"]])
        morgan = write(folder / "morgan.npy", b"morgan")
        return {"file": f"stage2/04_tsg_dose/execution/checkpoints/"
                        f"round1-axis-local-inputs/{source}",
                "population_sha256": sha(population), "splits_sha256": sha(splits),
                "parents_sha256": sha(parents), "morgan_sha256": sha(morgan)}

    def _freezes(self):
        for source in SOURCES:
            sealed = self._sealed(source)
            for view in VIEWS[source]:
                folder = self.step("execution/freezes") / view / source
                folder.mkdir(parents=True, exist_ok=True)
                files = {}
                for name in ("config.json", "population.csv", "all_obs.csv", "parents.csv",
                             "genes.csv", "splits.csv", "matched_controls.json",
                             "expression.npy", "control.npy", "morgan.npy",
                             "control_sources.csv"):
                    files[name] = write(folder / name, f"{view}-{source}-{name}")
                view_manifest = self.step("output/views") / view / f"{source}.json"
                mode = ("delta" if view in ("full-pooled-label", "solo-deepest")
                        else ("alias" if view == "full" else "payload"))
                write_json(folder / "manifest.json", {
                    "schema": "tsg-view-freeze-v1",
                    "delta_kind": "pooled" if view == "full-pooled-label"
                                  else "expression_minus_control",
                    "outputs": {name: sha(path) for name, path in files.items()},
                    "sealed": sealed,
                    "view_product": {"file": f"stage2/tsg_view/output/views/{view}/"
                                             f"{source}.json",
                                     "sha256": sha(view_manifest), "mode": mode},
                })

    # -- frozen truth ----------------------------------------------------- #
    def _truth(self):
        for source in SOURCES:
            folder = self.step("output/truth") / source
            noise = csv(folder / "noise_floor.csv", ["condition", "noise"], [["c", 0.3]])
            diagnostics = write_json(folder / "diagnostics.json", {
                "schema": "tsg-view-truth-diagnostics-v1", "source": source,
                "label": "heuristic diagnostic",
                "noise_floor": {"median": 0.34, "unit": "same as Delta RMSE"},
                "pcc_ceiling": {"value": 0.38, "label": "heuristic diagnostic"},
                "file": {"name": "noise_floor.csv", "sha256": sha(noise), "rows": 1,
                         "bytes": noise.stat().st_size},
            })
            primary, condition = {}, {}
            for fold in (0, 1, 2):
                primary[f"fold-{fold}"] = self._truth_fold(
                    folder / "primary" / f"fold-{fold}",
                    ("population.csv", "delta.npy", "control.npy", "morgan.npy",
                     "parents.csv"))
                condition[f"fold-{fold}"] = self._truth_fold(
                    folder / "condition_mean" / f"fold-{fold}",
                    ("population.csv", "delta.npy"))
            write_json(folder / "manifest.json", {
                "schema": "tsg-view-truth-v1", "source": source,
                "scenario": "cold-drug-dose", "study": "stage2/tsg_view",
                "view_independent": True, "folds": [0, 1, 2],
                "primary": primary, "condition_mean": condition,
                "diagnostics": {"file": "diagnostics.json", "sha256": sha(diagnostics),
                                "label": "heuristic diagnostic"},
                "noise_floor": {"file": "noise_floor.csv", "sha256": sha(noise),
                                "bytes": noise.stat().st_size},
                "counts": {"primary": {"fold-0": {"rows": 4}}},
            })

    def _truth_fold(self, folder, names):
        outputs = {}
        for name in names:
            outputs[name] = record(write(folder / name, f"{folder.name}-{name}"))
        manifest = write_json(folder / "manifest.json",
                              {"schema": "tsg-view-truth-fold-v1",
                               "outputs": {name: record(folder / name) for name in names}})
        return {"rows": 4, "conditions": 2,
                "manifest": record(manifest),
                "outputs": outputs}

    # -- runs ------------------------------------------------------------- #
    def _runs(self):
        runs = self.step("output/runs")
        plan = {"schema": "tsg-view-formal-plan-v1", "mode": "formal",
                "origins": ORIGINS}
        plan_file = write_json(runs / "formal_plan.json", plan)
        plan_sha = sha(plan_file)
        write_json(runs / "formal_plan_manifest.json", {
            "schema": "tsg-view-formal-plan-manifest-v1", "mode": "formal",
            "origins": len(ORIGINS), "plan_sha256": plan_sha})
        write_json(runs / "formal_cap.json", {
            "schema": "tsg-view-formal-cap-v1", "memory_max": "10G",
            "memory_swap_max": "8G", "unit": "arc-tsg-view-formal-1",
            "scope": {"Id": "arc-tsg-view-formal-1.scope",
                      "InvocationID": "1" * 32,
                      "MemoryMax": "10737418240", "MemorySwapMax": "8589934592"},
            "threads": 2})
        present = ORIGINS if self.complete else [row for row in ORIGINS
                                                 if row["source"] == "tenv"]
        verify_rows = []
        for origin in present:
            self._run(origin, runs, plan_sha)
            verify_rows.append({
                "origin": origin,
                "run_manifest_sha256": sha(self._run_dir(runs, origin) / "manifest.json"),
                "selection_sha256": sha(self._selection(runs, origin)),
            })
        write_json(runs / "run_records.json", {
            "schema": "tsg-view-formal-run-v1", "mode": "formal",
            "plan_sha256": plan_sha, "selected": len(ORIGINS), "runs": len(present),
            "complete": len(present)})
        if self.complete:
            write_json(runs / "formal_verify.json", {
                "schema": "tsg-view-formal-verify-v1", "status": "PASS",
                "mode": "formal", "plan_sha256": plan_sha, "origins": len(present),
                "code_identity": CODE_IDENTITY,
                "full_arm": {"anchor_seed": 42, "origins": [], "fold_mean": [],
                             "sensitivity": [{"seed": 43, "note": SENSITIVITY_NOTE}]},
                "runs": verify_rows})

    def _run_dir(self, runs, origin):
        return (runs / "runs" / origin["source"] / origin["view"]
                / f"fold-{origin['fold']}" / f"seed-{origin['seed']}")

    def _selection(self, runs, origin):
        return runs / "selection" / origin["source"] / f"fold-{origin['fold']}.json"

    def _run(self, origin, runs, plan_sha):
        folder = self._run_dir(runs, origin)
        folder.mkdir(parents=True, exist_ok=True)
        files = {}
        for name in ("config.json", "checkpoint.pt", "learning_curve.csv",
                     "validation_prediction.npy", "validation_metrics.json",
                     "validation_manifest.json", "terminal.json", "test_predictions.npy",
                     "test_per_sample.csv", "test_metrics.json"):
            files[name] = write(folder / name, f"{origin}-{name}")
        terminal = write_json(folder / "terminal.json", {
            "status": "complete", "device": "cuda", "best_epoch": 7, "epochs": 9,
            "seed": origin["seed"], "fold": origin["fold"], "source": origin["source"],
            "view": origin["view"], "validation_mse": 0.25, "peak_rss_kib": 1000,
            "replay": "exact fresh-process native sampled prediction",
            "replay_pid": 100, "training_pid": 101,
            "training_process_token": "a" * 64, "test_evaluation_process_token": "b" * 64,
            "replay_process_token": "c" * 64, "validation_replay_token": "d" * 64,
            "test_replay_token": "c" * 64, "test_samples": 4})
        outputs = {name: sha(path) for name, path in files.items()}
        outputs["terminal.json"] = sha(terminal)
        truth_shas = self._truth_shas(origin["source"], origin["fold"])
        manifest = write_json(folder / "manifest.json", {
            "schema": "tsg-view-run-v1", "mode": "formal", "origin": origin,
            "plan": plan_sha,
            "freeze": {"manifest_sha256": self._freeze_sha(origin["source"],
                                                           origin["view"]),
                       "sealed_splits_sha256": "f" * 64,
                       "view_product_sha256": self._product_sha(origin["source"],
                                                                origin["view"])},
            "truth": {"primary": {"manifest_sha256": truth_shas["primary"],
                                  "path": "p", "rows": 4},
                      "condition_mean": {"manifest_sha256": truth_shas["condition_mean"],
                                         "path": "c", "rows": None}},
            "outputs": outputs,
        })
        selection = write_json(self._selection(runs, origin), {
            "schema": "tsg-view-validation-selection-v1", "source": origin["source"],
            "fold": origin["fold"], "plan_sha256": plan_sha, "origins": [],
            "truth": {"primary": {"manifest_sha256": "1" * 64, "path": "p", "rows": 4}}})
        write_json(runs / "records" / origin["source"] / origin["view"]
                   / f"fold-{origin['fold']}" / f"seed-{origin['seed']}.json", {
            "schema": "tsg-view-run-record-v1", "mode": "formal", "origin": origin,
            "seed": origin["seed"], "fold": origin["fold"], "device": "cuda",
            "plan_sha256": plan_sha, "selection_sha256": sha(selection),
            "best_epoch": 7, "realized_epochs": 9, "validation_mse": 0.25,
            "runtime_seconds": 12.0, "training_seconds": 9.0, "peak_rss_kib": 1000,
            "peak_gpu_used_mib": 2000, "artifacts": outputs,
            "metrics": {"primary": {"pcc": 0.3, "rmse": 0.4, "n_samples": 4},
                        "condition_mean": {"pcc": 0.3, "rmse": 0.4, "n_samples": 2}},
            "truth": {"primary": {"manifest_sha256": self._truth_shas(
                origin["source"], origin["fold"])["primary"], "path": "p", "rows": 4}},
            "phases": []})
        self._resources(runs, origin, plan_sha)
        return manifest

    def _resources(self, runs, origin, plan_sha):
        base = (runs / "resources" / origin["source"] / origin["view"]
                / f"fold-{origin['fold']}" / f"seed-{origin['seed']}")
        phases = (("train", "train", "training"), ("replay", "validation", "validation"),
                  ("evaluate", "test", "test"), ("replay", "test", "test"))
        for index, (action, phase, tag) in enumerate(phases):
            folder = base / f"{phase}-{tag}-{index}"
            folder.mkdir(parents=True, exist_ok=True)
            start = 1795343861.0 + index
            write_json(folder / "phase.json", {
                "action": action, "phase": phase, "origin": origin,
                "memory_max": "10G", "memory_swap_max": "8G", "threads": 2,
                "plan_sha256": plan_sha, "started_utc_unix_seconds": start,
                "available_ram_kib": 20000000})
            samples = [{"utc_unix_seconds": start + 10 * step,
                        "mem_total_kib": 24612504, "mem_available_kib": 20000000 - step,
                        "swap_free_kib": 60000000,
                        "gpu_used_mib": 2000 + step, "gpu_free_mib": 14000,
                        "gpu_util_pct": 5} for step in range(3)]
            write(folder / "telemetry.jsonl",
                  "".join(json.dumps(sample, sort_keys=True) + "\n" for sample in samples))
            write(folder / "time.txt",
                  "\tCommand being timed: \"worker\"\n"
                  "\tMaximum resident set size (kbytes): 1000\n"
                  "\tExit status: 0\n")
            write_json(folder / "result.json", {
                "elapsed_seconds": 30.0, "peak_rss_kib": 1000, "plan_mode": "formal",
                "returncode": 0, "time_exit_status": 0,
                "telemetry": {"samples": 3, "peak_gpu_used_mib": 2002,
                              "min_mem_available_kib": 19999998}})
            write(folder / "subprocess.log", "ok\n")

    # -- report + figures ------------------------------------------------- #
    def _report(self):
        write(self.step("execution/report.py"), "#!/usr/bin/env python3\n")
        report = self.step("output/report")
        tables = {}
        for name in ("metrics_primary.csv", "metrics_condition_mean.csv",
                     "metrics_per_dose.csv", "diagnostics_beside_scores.csv",
                     "fold_spread.csv", "learning_curves.csv", "paired_vs_full.csv",
                     "paired_vs_full_solo_deepest.csv", "baseline_mean_cell.csv",
                     "run_status.csv", "seed_spread.csv"):
            tables[name] = sha(csv(report / name, ["metric", "value"], [["pcc", 0.3]]))
        figure_data = {}
        for name in ("learning_curves.csv", "metrics_bars.csv", "metrics_points.csv",
                     "views_effect.csv"):
            figure_data[f"figure/{name}"] = sha(csv(self.step("figure/data") / name,
                                                    ["panel", "value"], [["a", 0.3]]))
        write_json(report / "report_manifest.json", {
            "schema": "tsg-view-report-v1",
            "sources": list(SOURCES) if self.complete else ["tenv"],
            "complete_origins": len(ORIGINS) if self.complete else 27,
            "pending_origins": 0 if self.complete else 18,
            "tables": tables, "figure_data": figure_data,
            "estimands": {"fold": "equal-fold mean"}})

    def _figures(self):
        figure = self.step("figure")
        renderer = write(self.step("execution/render.py"), RENDERER)
        # The replay copies the declared closure, so the fixture tree must carry a sibling
        # module for every declared entry even though this stand-in renderer imports none.
        write(self.step("execution/common.py"),
              "SOURCES = (\"tenv\", \"tahoe_bulk\")\nVIEW_ORDER = {}\n"
              "def seeds_for(view):\n    return (42,)\n")
        write(self.step("figure/script/render.py"),
              "#!/usr/bin/env python3\nfrom pathlib import Path\n"
              "import sys\nsys.path.insert(0, str(Path(__file__).resolve().parents[2]"
              " / 'execution'))\nimport render\nraise SystemExit(render.main())\n")
        self._report_figure_data()
        env = {"PATH": "/usr/bin:/bin", "MPLCONFIGDIR": str(figure / ".mpl"),
               "PYTHONHASHSEED": "0", "HOME": str(figure / ".home")}
        subprocess.run([sys.executable, "-I", str(renderer), "--data", str(figure / "data"),
                        "--output", str(figure), "--suite", "all" if self.complete else "tenv"],
                       check=True, env=env)

    def _report_figure_data(self):
        # figure/data is the plotted subset duplicated from the report tables
        data = self.step("figure/data")
        for name in ("learning_curves.csv", "metrics_bars.csv", "metrics_points.csv",
                     "views_effect.csv"):
            if not (data / name).is_file():
                csv(data / name, ["panel", "value"], [["a", 0.3]])

    # -- health ----------------------------------------------------------- #
    def _health(self, *, verify_all=True):
        script = self.root / "script"
        write(script / "check_paths.sh", "#!/bin/sh\necho 'RESULT: PASS'\nexit 0\n")
        write(script / "check_tracked.sh", "#!/bin/sh\necho 'RESULT: PASS'\nexit 0\n")
        folder = self.step("log/acceptance/health")
        checks = {}
        for name in ("check_paths", "check_tracked"):
            log = write(folder / f"{name}.sh.log", f"== {name}: PASS ==\nexit: 0\n")
            checks[name + ".sh"] = {
                "returncode": 0, "script_sha256": sha(script / f"{name}.sh"),
                "log_sha256": sha(log), "runtime_seconds": 0.5, "command":
                ["bash", f"script/{name}.sh"]}
        if verify_all:
            log = write(folder / "verify_all.sh.log", "RESULT: PASS - verified 58/58\n")
            checks["verify_all.sh"] = {
                "returncode": 0,
                "script_sha256": sha(write(script / "verify_all.sh",
                                          "#!/bin/sh\necho pass\nexit 0\n")),
                "log_sha256": sha(log), "runtime_seconds": 273.0,
                "command": ["bash", "script/verify_all.sh"]}
        write_json(folder / "health.json", {
            "schema": "tsg-view-health-v1", "commit": self.head(), "checks": checks})

    # -- git + baseline --------------------------------------------------- #
    def _git(self):
        env = {"GIT_AUTHOR_NAME": "fixture", "GIT_AUTHOR_EMAIL": "f@example.test",
               "GIT_COMMITTER_NAME": "fixture", "GIT_COMMITTER_EMAIL": "f@example.test",
               "PATH": "/usr/bin:/bin", "HOME": str(self.root / ".home")}
        run = lambda *args: subprocess.run(["git", "-C", str(self.root), *args],
                                           check=True, env=env, capture_output=True)
        run("init", "-q", "-b", "dev")
        run("add", "-A")
        run("commit", "-q", "-m", "fixture")
        head = subprocess.run(["git", "-C", str(self.root), "rev-parse", "HEAD"],
                              check=True, capture_output=True, text=True).stdout.strip()
        write(self.root / ".scratch/tsg-view/research/baseline/refs.txt",
              head + "\n" + head + "\n")

    def _baseline(self, mode):
        lines = []
        for name in STEPS:
            folder = self.root / "stage2" / name
            for path in sorted(p for p in folder.rglob("*") if p.is_file()):
                lines.append(f"{sha(path)}  {path.relative_to(self.root).as_posix()}")
        baseline = (self.root / ".scratch/tsg-view/research/baseline/"
                    "steps01-06-per-path.sha256")
        if mode == "manifest":
            write(baseline, "\n".join(lines) + "\n")
        else:
            write(baseline, hashlib.sha256("\n".join(lines).encode()).hexdigest() + "  -\n")


def build(root, *, complete=True, **kwargs):
    return Fixture(root, complete=complete).build(**kwargs)


class AcceptanceTestCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.root = self.tmp / "arc"
        self.root.mkdir()
        build(self.root, **getattr(self, "fixture_kwargs", {}))

    def tearDown(self):
        self._tmp.cleanup()

    def receipt(self, **kwargs):
        options = {"phase": "final", "timestamp": STAMP, "environment": ENVIRONMENT,
                   "commands": [{"name": "receipt", "argv": [
                       "python", "stage2/tsg_view/execution/acceptance.py",
                       "receipt", "--phase", "final"]}],
                   "verify_binary": True, "preregistration_text": PREREGISTRATION}
        options.update(kwargs)
        return acceptance.build_receipt(self.root, **options)

    def git_head(self):
        """The fixture commit, read straight from git rather than from ``git_revs``."""
        completed = subprocess.run(["git", "-C", str(self.root), "rev-parse", "HEAD"],
                                   check=True, capture_output=True, text=True)
        return completed.stdout.strip()


PREREGISTRATION = """Pre-registered partial-contribution test.

A tenv Delta PCC gain from `merge` or `full-pooled-label` of at most 0.12 that stays
below the single-sample heuristic ceiling (about 0.38), with the tenv-Tahoe gap still
present, supports the partial-contribution reading. Any other outcome revises it.
"""


class InventoryTest(AcceptanceTestCase):
    def test_receipt_inventories_every_declared_product(self):
        receipt = self.receipt()
        self.assertEqual(receipt["schema"], "tsg-view-acceptance-v1")
        self.assertEqual(receipt["status"], "passed")
        self.assertEqual(receipt["formal_fits"], 45)
        self.assertEqual(len(receipt["origins"]), 45)
        self.assertEqual(receipt["seeds"],
                         {"main": [42, 43], "solo_deepest": [42], "solo_draw": 20261122})
        self.assertEqual(len(receipt["products"]["views"]), 9)
        self.assertEqual(len(receipt["products"]["freezes"]), 9)
        self.assertEqual(sorted(receipt["products"]["truth"]), sorted(SOURCES))
        self.assertEqual(len(receipt["products"]["report"]["tables"]), 11)
        self.assertEqual(len(receipt["products"]["figures"]), 7)
        self.assertTrue(receipt["caps"]["memory_max"] == "10G")
        self.assertEqual(receipt["csv_only_replay"]["byte_identical"], True)
        self.assertEqual(receipt["replay"]["all_exact"], True)
        # every inventoried path is repo-relative and digest-shaped
        for group in receipt["inventory"].values():
            for path, value in group.items():
                self.assertFalse(path.startswith("/"), path)
                self.assertFalse(path.startswith(".."), path)
                self.assertRegex(value["sha256"], r"^[0-9a-f]{64}$")
                self.assertGreater(value["bytes"], 0)

    def test_receipt_never_contains_absolute_tmp_paths(self):
        receipt = self.receipt()
        blob = json.dumps(receipt)
        self.assertNotIn(str(self.tmp), blob)
        self.assertNotIn("/tmp", blob)
        self.assertNotIn(str(self.root), blob)

    def test_receipt_is_deterministic_modulo_timestamp(self):
        first = self.receipt()
        second = self.receipt()
        self.assertEqual(acceptance.receipt_bytes(first), acceptance.receipt_bytes(second))
        later = self.receipt(timestamp="2026-11-23T00:00:00Z")
        first_without = dict(first)
        later_without = dict(later)
        first_without.pop("generated_at")
        later_without.pop("generated_at")
        self.assertEqual(first_without, later_without)
        self.assertNotEqual(acceptance.receipt_bytes(first),
                            acceptance.receipt_bytes(later))

    def test_preregistration_digest_is_recorded(self):
        receipt = self.receipt()
        self.assertEqual(receipt["pre_registration_sha256"],
                         hashlib.sha256(PREREGISTRATION.encode()).hexdigest())


class FailClosedTest(AcceptanceTestCase):
    def test_view_payload_digest_mismatch_fails_closed(self):
        payload = self.root / "stage2/tsg_view/output/views/merge/tenv.h5ad"
        write(payload, "corrupted")
        with self.assertRaisesRegex(acceptance.AcceptanceError, "digest mismatch"):
            self.receipt()

    def test_missing_artifacts_fail_closed_in_final(self):
        (self.root / "stage2/tsg_view/output/truth/tenv/noise_floor.csv").unlink()
        with self.assertRaisesRegex(acceptance.AcceptanceError, "missing"):
            self.receipt()

    def test_unknown_origin_manifest_fails_closed(self):
        manifest = self.root / ("stage2/tsg_view/output/runs/runs/tenv/full/fold-0/"
                                "seed-42/manifest.json")
        value = json.loads(manifest.read_text())
        value["plan"] = "0" * 64
        write_json(manifest, value)
        with self.assertRaisesRegex(acceptance.AcceptanceError, "plan"):
            self.receipt()

    def test_cap_violation_fails_closed(self):
        phase = next((self.root / "stage2/tsg_view/output/runs/resources").rglob("phase.json"))
        value = json.loads(phase.read_text())
        value["memory_max"] = "12G"
        write_json(phase, value)
        with self.assertRaisesRegex(acceptance.AcceptanceError, "cap"):
            self.receipt()

    def test_peak_rss_above_cap_fails_closed(self):
        phase = next((self.root / "stage2/tsg_view/output/runs/resources").rglob("time.txt"))
        write(phase, "\tMaximum resident set size (kbytes): 999999999999\n"
                     "\tExit status: 0\n")
        with self.assertRaisesRegex(acceptance.AcceptanceError, "cap"):
            self.receipt()

    def test_figure_digest_mismatch_fails_closed(self):
        png = self.root / "stage2/tsg_view/figure/metrics_tenv.png"
        write(png, b"not the pinned png")
        with self.assertRaisesRegex(acceptance.AcceptanceError, "figure"):
            self.receipt()

    def test_tampered_figure_data_fails_closed(self):
        data = self.root / "stage2/tsg_view/figure/data/metrics_bars.csv"
        write(data, data.read_text() + "extra,1\n")
        with self.assertRaisesRegex(acceptance.AcceptanceError, "figure"):
            self.receipt()

    def test_missing_replay_token_fails_closed(self):
        run = self.root / "stage2/tsg_view/output/runs/runs/tenv/full/fold-0/seed-42"
        terminal = run / "terminal.json"
        value = json.loads(terminal.read_text())
        value.pop("test_replay_token")
        write_json(terminal, value)
        manifest = json.loads((run / "manifest.json").read_text())
        manifest["outputs"]["terminal.json"] = sha(terminal)
        write_json(run / "manifest.json", manifest)
        record_file = (self.root / "stage2/tsg_view/output/runs/records/tenv/full/fold-0/"
                       "seed-42.json")
        record_value = json.loads(record_file.read_text())
        record_value["artifacts"]["terminal.json"] = sha(terminal)
        write_json(record_file, record_value)
        with self.assertRaisesRegex(acceptance.AcceptanceError, "replay"):
            self.receipt()

    def test_missing_health_transcript_fails_closed(self):
        (self.root / "stage2/tsg_view/log/acceptance/health/check_paths.sh.log").unlink()
        with self.assertRaisesRegex(acceptance.AcceptanceError, "health"):
            self.receipt()


class TelemetryTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def phase(self, *, samples=3, memory_max="10G", peak_rss_kib=1000, gap=10.0,
              returncode=0):
        folder = self.tmp / "phase"
        folder.mkdir(exist_ok=True)
        write_json(folder / "phase.json", {
            "action": "evaluate", "phase": "test",
            "memory_max": memory_max, "memory_swap_max": "8G",
            "origin": {"source": "tenv", "view": "full", "fold": 0, "seed": 42,
                       "method": "no-dose", "scenario": "cold-drug-dose"},
            "started_utc_unix_seconds": 100.0, "plan_sha256": "a" * 64, "threads": 2})
        rows = [{"utc_unix_seconds": 100.0 + index * gap,
                 "mem_total_kib": 24612504, "mem_available_kib": 20000000 - index,
                 "swap_free_kib": 60000000, "gpu_used_mib": 2000 + index,
                 "gpu_free_mib": 14000, "gpu_util_pct": 7} for index in range(samples)]
        write(folder / "telemetry.jsonl",
              "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows))
        write(folder / "time.txt",
              f"\tMaximum resident set size (kbytes): {peak_rss_kib}\n"
              f"\tExit status: {returncode}\n")
        write_json(folder / "result.json", {
            "elapsed_seconds": 30.0, "peak_rss_kib": peak_rss_kib, "returncode": returncode,
            "time_exit_status": returncode, "plan_mode": "formal",
            "telemetry": {"samples": samples}})
        return folder

    def test_telemetry_is_summarised_with_gaps_and_peaks(self):
        summary = acceptance.summarise_phase(self.phase())
        self.assertEqual(summary["samples"], 3)
        self.assertEqual(summary["max_gap_seconds"], 10.0)
        self.assertEqual(summary["peak_gpu_used_mib"], 2002)
        self.assertEqual(summary["min_available_ram_kib"], 19999998)
        self.assertEqual(summary["peak_rss_kib"], 1000)
        self.assertEqual(summary["memory_max_bytes"], 10 * 1024 ** 3)
        self.assertEqual(summary["swap_max_bytes"], 8 * 1024 ** 3)
        self.assertEqual(summary["returncode"], 0)

    def test_cap_assertion_rejects_a_wrong_cap(self):
        summary = acceptance.summarise_phase(self.phase(memory_max="12G"))
        with self.assertRaisesRegex(acceptance.AcceptanceError, "cap"):
            acceptance.assert_caps(summary)

    def test_cap_assertion_rejects_rss_over_the_cap(self):
        summary = acceptance.summarise_phase(self.phase(peak_rss_kib=999999999))
        with self.assertRaisesRegex(acceptance.AcceptanceError, "cap"):
            acceptance.assert_caps(summary)

    def test_cap_assertion_rejects_a_stalled_sampler(self):
        summary = acceptance.summarise_phase(self.phase(gap=60.0))
        with self.assertRaisesRegex(acceptance.AcceptanceError, "telemetry"):
            acceptance.assert_caps(summary)


class PhaseTest(AcceptanceTestCase):
    fixture_kwargs = {"complete": False, "verify_all": False}

    def test_partial_phase_reports_completeness(self):
        receipt = self.receipt(phase="partial")
        self.assertEqual(receipt["status"], "partial")
        self.assertEqual(receipt["complete_origins"], 27)
        self.assertEqual(len(receipt["pending_origins"]), 18)
        self.assertEqual(len(receipt["origins"]), 27)
        self.assertEqual(sorted({row["source"] for row in receipt["pending_origins"]}),
                         ["tahoe_bulk"])
        self.assertEqual(receipt["csv_only_replay"]["suite"], "tenv")
        self.assertEqual(receipt["health"]["verify_all"]["status"],
                         "documented-not-run")

    def test_final_phase_fails_on_an_incomplete_matrix(self):
        with self.assertRaisesRegex(acceptance.AcceptanceError, "incomplete|missing"):
            self.receipt()



class FinalGateTest(AcceptanceTestCase):
    def test_final_phase_requires_the_preregistration(self):
        with self.assertRaisesRegex(acceptance.AcceptanceError, "pre-registration"):
            self.receipt(preregistration_text=None)

    def test_final_phase_requires_verify_all_evidence(self):
        health = self.root / "stage2/tsg_view/log/acceptance/health/health.json"
        value = json.loads(health.read_text())
        value["checks"].pop("verify_all.sh")
        write_json(health, value)
        with self.assertRaisesRegex(acceptance.AcceptanceError, "verify_all"):
            self.receipt()


class FigureReplayTest(AcceptanceTestCase):
    def test_csv_only_replay_is_byte_identical(self):
        value = acceptance.figure_replay(self.root, self.step_figure())
        self.assertTrue(value["byte_identical"])
        self.assertEqual(len(value["figures"]), 7)
        self.assertEqual(value["suite"], "all")

    def test_nondeterministic_renderer_fails_closed(self):
        renderer = self.root / "stage2/tsg_view/execution/render.py"
        write(renderer, RENDERER.replace(
            "body = digest.hexdigest().encode() + b\":\" + fig.encode()",
            "import os\n        body = (digest.hexdigest() + str(os.getpid()))"
            ".encode() + b\":\" + fig.encode()"))
        with self.assertRaisesRegex(acceptance.AcceptanceError, "re-render"):
            acceptance.figure_replay(self.root, self.step_figure())

    def step_figure(self):
        return self.root / "stage2/tsg_view/figure"


class ClosureTest(unittest.TestCase):
    """The declared replay closure must cover what ``execution/render.py`` imports.

    ``figure_replay`` runs the renderer as a closed world under ``python -I``; an import of
    an ``execution`` module that is not declared in ``RENDERER_MODULES`` cannot be satisfied
    there.  This test re-derives the closure from the source, so the declaration cannot rot.
    """

    def test_declared_closure_covers_the_real_renderer_imports(self):
        closure = acceptance.renderer_import_closure(EXECUTION)
        declared = set(acceptance.RENDERER_MODULES)
        missing = sorted(set(closure) - declared)
        self.assertFalse(missing,
                         "execution/render.py needs modules missing from "
                         f"RENDERER_MODULES: {missing}")
        self.assertEqual(sorted(closure), sorted(declared))

    def test_closure_walk_is_transitive_and_skips_outside_imports(self):
        with tempfile.TemporaryDirectory() as temporary:
            execution = Path(temporary)
            write(execution / "render.py",
                  "import json\nfrom pathlib import Path\nimport common\n")
            write(execution / "common.py", "import helper\nimport matplotlib\n")
            write(execution / "helper.py", "import csv\n")
            write(execution / "unrelated.py", "pass\n")
            self.assertEqual(acceptance.renderer_import_closure(execution),
                             ("common.py", "helper.py", "render.py"))


class ReceiptClosureTest(AcceptanceTestCase):
    def test_receipt_pins_every_closure_module_digest(self):
        receipt = self.receipt()
        block = receipt["csv_only_replay"]["renderer_closure_sha256"]
        execution = self.root / "stage2/tsg_view/execution"
        self.assertEqual(sorted(block), sorted(acceptance.RENDERER_MODULES))
        for name in acceptance.RENDERER_MODULES:
            self.assertRegex(block[name], r"^[0-9a-f]{64}$")
            self.assertEqual(block[name], sha(execution / name), name)
        self.assertEqual(receipt["csv_only_replay"]["renderer_sha256"], block["render.py"])
        self.assertTrue(receipt["csv_only_replay"]["byte_identical"])


class CodeIdentityTest(AcceptanceTestCase):
    def test_code_drift_identity_is_carried_verbatim(self):
        block = self.receipt()["formal_verify"]
        self.assertEqual(block["status"], "PASS")
        self.assertEqual(block["origins"], len(ORIGINS))
        self.assertEqual(block["sha256"],
                         sha(self.root / "stage2/tsg_view/output/runs/formal_verify.json"))
        identity = block["code_identity"]
        self.assertEqual(identity["status"], "CODE_DRIFT")
        self.assertEqual(identity["moved"], CODE_DRIFT)
        self.assertEqual(identity["note"], CODE_IDENTITY["note"])
        self.assertEqual(identity["driver_code_drift"], ["formal.py"])
        full_arm = block["full_arm"]
        self.assertEqual(full_arm["anchor_seed"], 42)
        self.assertEqual(full_arm["sensitivity_seeds"], [43])
        self.assertIn("seed 42", full_arm["statement"])
        self.assertIn("seed-43", full_arm["statement"])

    def test_missing_formal_verify_fails_closed_in_final(self):
        (self.root / "stage2/tsg_view/output/runs/formal_verify.json").unlink()
        with self.assertRaisesRegex(acceptance.MissingArtifact, "formal verification"):
            self.receipt()

    def test_non_pass_formal_verify_fails_closed(self):
        path = self.root / "stage2/tsg_view/output/runs/formal_verify.json"
        value = json.loads(path.read_text())
        value["status"] = "DRIFT"
        write_json(path, value)
        with self.assertRaisesRegex(acceptance.AcceptanceError, "status 'DRIFT'"):
            self.receipt()

    def test_identity_status_must_agree_with_the_moved_digests(self):
        path = self.root / "stage2/tsg_view/output/runs/formal_verify.json"
        value = json.loads(path.read_text())
        value["code_identity"]["status"] = "MATCH"
        write_json(path, value)
        with self.assertRaisesRegex(acceptance.AcceptanceError, "code identity"):
            self.receipt()


class HealthFreshnessTest(AcceptanceTestCase):
    def test_health_block_records_the_commit_against_head(self):
        receipt = self.receipt()
        health = receipt["health"]
        self.assertEqual(health["head"], self.git_head())
        self.assertEqual(health["commit"], self.git_head())
        self.assertFalse(health["health_commit_stale"])

    def test_stale_health_commit_fails_closed_in_final(self):
        path = self.root / "stage2/tsg_view/log/acceptance/health/health.json"
        value = json.loads(path.read_text())
        value["commit"] = "0" * 40
        write_json(path, value)
        with self.assertRaisesRegex(acceptance.AcceptanceError,
                                    "health transcripts predate HEAD"):
            self.receipt()


class StaleHealthPartialTest(AcceptanceTestCase):
    fixture_kwargs = {"complete": False, "verify_all": False}

    def test_partial_phase_records_stale_health_instead_of_failing(self):
        path = self.root / "stage2/tsg_view/log/acceptance/health/health.json"
        value = json.loads(path.read_text())
        value["commit"] = "0" * 40
        write_json(path, value)
        receipt = self.receipt(phase="partial")
        self.assertTrue(receipt["health"]["health_commit_stale"])
        self.assertEqual(receipt["health"]["commit"], "0" * 40)
        self.assertEqual(receipt["health"]["head"], self.git_head())


class RealRendererReplayTest(unittest.TestCase):
    """The replay seam against the real renderer and the real ``figure/data`` payloads.

    ``FigureReplayTest`` runs a synthetic renderer that imports nothing, which is exactly
    how the ``render.py`` -> ``common`` regression slipped through.  This test renders with
    the real module and the real CSVs, then asks ``figure_replay`` to reproduce those bytes
    from a closure-only copy under ``python -I``; an import the copy forgets fails here.
    """

    def test_csv_only_replay_reproduces_the_real_renderer(self):
        figure = REPO / "stage2/tsg_view/figure"
        source = REPO / "stage2/tsg_view/execution/render.py"
        data = sorted((figure / "data").glob("*.csv"))
        if not data or not source.is_file():
            self.skipTest(f"real figure payloads are absent: {figure / 'data'}")
        with tempfile.TemporaryDirectory() as temporary:
            staged = Path(temporary) / "figure"
            (staged / "data").mkdir(parents=True)
            (staged / "home").mkdir()
            for path in data:
                shutil.copyfile(path, staged / "data" / path.name)
            environment = {"PATH": "/usr/bin:/bin",
                           "MPLCONFIGDIR": str(staged / "mpl"),
                           "PYTHONHASHSEED": "0",
                           "PYTHONDONTWRITEBYTECODE": "1",
                           "HOME": str(staged / "home")}
            completed = subprocess.run(
                [sys.executable, "-I", str(source), "--data", str(staged / "data"),
                 "--output", str(staged), "--suite", "all"],
                env=environment, capture_output=True, text=True, timeout=600)
            self.assertEqual(completed.returncode, 0, completed.stderr[-2000:])
            status = json.loads((staged / "render_status.json").read_text())
            replayed = acceptance.figure_replay(REPO, staged)
            self.assertTrue(replayed["byte_identical"])
            self.assertEqual(replayed["figures"], status["figures"])
            self.assertGreaterEqual(len(replayed["figures"]), 1)
            closure = replayed["renderer_closure_sha256"]
            self.assertEqual(sorted(closure), sorted(acceptance.RENDERER_MODULES))
            self.assertEqual(closure["render.py"], sha(source))


class InvarianceTest(AcceptanceTestCase):
    def baseline(self):
        return (self.root / ".scratch/tsg-view/research/baseline/"
                "steps01-06-per-path.sha256")

    def test_invariance_manifest_baseline_passes(self):
        value = acceptance.steps_invariance(self.root, self.baseline())
        self.assertEqual(value["verdict"], "PASS")
        self.assertEqual(sorted(value["steps"]), list(STEPS))
        for step in STEPS:
            self.assertEqual(value["steps"][step]["verdict"], "PASS")
            self.assertEqual(value["steps"][step]["changed"], [])
        self.assertEqual(value["baseline_file"]["mode"], "manifest")
        self.assertGreater(value["baseline_file"]["files"], 0)
        self.assertEqual(value["git_diff"]["status"], "PASS")
        self.assertTrue(value["git_diff"]["empty"])
        self.assertTrue(value["git_diff"]["worktree_clean"])

    def test_invariance_detects_a_committed_change_outside_the_baseline(self):
        # A source file the artifact manifest never listed is still caught by the
        # pinned-planning-ref diff.
        write(self.root / "stage2/01_bases/execution/late.py", "x = 1\n")
        env = {"GIT_AUTHOR_NAME": "f", "GIT_AUTHOR_EMAIL": "f@e",
               "GIT_COMMITTER_NAME": "f", "GIT_COMMITTER_EMAIL": "f@e",
               "PATH": "/usr/bin:/bin", "HOME": str(self.root / ".home")}
        subprocess.run(["git", "-C", str(self.root), "add", "-A"], check=True,
                       capture_output=True)
        subprocess.run(["git", "-C", str(self.root), "commit", "-q", "-m", "late"],
                       check=True, capture_output=True, env=env)
        value = acceptance.steps_invariance(self.root, self.baseline())
        self.assertEqual(value["verdict"], "CHANGED")
        block = value["steps"]["01_bases"]
        self.assertEqual(block["verdict"], "CHANGED")
        self.assertEqual(block["changed"], [])
        self.assertEqual(block["committed_changes"],
                         ["stage2/01_bases/execution/late.py"])
        self.assertEqual(value["git_diff"]["status"], "CHANGED")

    def test_invariance_detects_an_uncommitted_change(self):
        write(self.root / "stage2/04_tsg_dose/execution/edit.py", "y = 2\n")
        value = acceptance.steps_invariance(self.root, self.baseline())
        self.assertEqual(value["verdict"], "CHANGED")
        self.assertEqual(value["steps"]["04_tsg_dose"]["verdict"], "CHANGED")
        self.assertEqual(value["steps"]["04_tsg_dose"]["worktree_changes"],
                         ["stage2/04_tsg_dose/execution/edit.py"])
        self.assertFalse(value["git_diff"]["worktree_clean"])

    def test_invariance_detects_a_changed_tracked_file(self):
        write(self.root / "stage2/03_tsg_drug_emb/README.md", "# changed\n")
        value = acceptance.steps_invariance(self.root, self.baseline())
        self.assertEqual(value["verdict"], "CHANGED")
        self.assertEqual(value["steps"]["03_tsg_drug_emb"]["verdict"], "CHANGED")
        self.assertEqual(value["steps"]["01_bases"]["verdict"], "PASS")
        changed = value["steps"]["03_tsg_drug_emb"]["changed"]
        self.assertEqual([entry["path"] for entry in changed],
                         ["stage2/03_tsg_drug_emb/README.md"])
        self.assertNotEqual(changed[0]["baseline"], changed[0]["current"])

    def test_invariance_opaque_baseline_falls_back_to_the_pinned_ref(self):
        # The real capture is an opaque single digest; the ref in refs.txt is the
        # reproducible before/after anchor.
        lines = self.baseline().read_text().splitlines()
        self.baseline().write_text(
            hashlib.sha256("\n".join(lines).encode()).hexdigest() + "  -\n")
        value = acceptance.steps_invariance(self.root, self.baseline())
        self.assertEqual(value["verdict"], "PASS")
        self.assertEqual(value["baseline_file"]["mode"], "opaque")
        self.assertFalse(value["baseline_file"]["matches_recipe"])
        self.assertEqual(value["baseline_ref"]["rev"],
                         (self.root / ".scratch/tsg-view/research/baseline/refs.txt")
                         .read_text().split()[0])
        # a committed change moves HEAD past the pinned ref -> CHANGED
        write(self.root / "stage2/05_tsg_width_dose/README.md", "# changed\n")
        subprocess.run(["git", "-C", str(self.root), "add", "-A"], check=True,
                       capture_output=True)
        subprocess.run(["git", "-C", str(self.root), "commit", "-q", "-m", "change"],
                       check=True, capture_output=True,
                       env={"GIT_AUTHOR_NAME": "f", "GIT_AUTHOR_EMAIL": "f@e",
                            "GIT_COMMITTER_NAME": "f", "GIT_COMMITTER_EMAIL": "f@e",
                            "PATH": "/usr/bin:/bin", "HOME": str(self.root / ".home")})
        value = acceptance.steps_invariance(self.root, self.baseline())
        self.assertEqual(value["verdict"], "CHANGED")
        self.assertEqual(value["steps"]["05_tsg_width_dose"]["verdict"], "CHANGED")

    def test_receipt_pins_the_invariance_proof(self):
        receipt = self.receipt()
        proof = receipt["steps01_06"]
        self.assertEqual(proof["verdict"], "PASS")
        self.assertRegex(proof["baseline_file"]["sha256"], r"^[0-9a-f]{64}$")
        self.assertEqual(proof["baseline_ref"]["rev"], (self.root / ".scratch/tsg-view/"
                                                        "research/baseline/refs.txt")
                         .read_text().split()[0])
        self.assertGreater(proof["sealed"]["files"], 0)

    def test_receipt_pins_the_sealed_schema_and_git_evidence(self):
        receipt = self.receipt()
        self.assertEqual(receipt["schema"], "tsg-view-acceptance-v1")
        self.assertEqual(receipt["phase"], "final")
        self.assertEqual(receipt["status"], "passed")
        self.assertEqual(receipt["formal_fits"], 45)
        self.assertEqual(receipt["complete_origins"], 45)
        self.assertEqual(receipt["pending_origins"], [])
        self.assertTrue(receipt["csv_only_replay"]["byte_identical"])
        self.assertTrue(receipt["replay"]["all_exact"])
        self.assertEqual(receipt["caps"]["memory_max"], "10G")
        self.assertEqual(receipt["caps"]["memory_swap_max"], "8G")
        self.assertRegex(receipt["pre_registration_sha256"], r"^[0-9a-f]{64}$")
        self.assertRegex(receipt["plan_sha256"], r"^[0-9a-f]{64}$")
        self.assertRegex(receipt["input_manifest_sha256"], r"^[0-9a-f]{64}$")
        self.assertEqual(len(receipt["resources"]), 45)
        first = next(iter(receipt["resources"].values()))
        self.assertEqual(first["memory_max_bytes"], 10737418240)
        self.assertEqual(first["swap_max_bytes"], 8589934592)
        self.assertGreaterEqual(first["telemetry_samples"], 1)
        proof = receipt["steps01_06"]
        self.assertEqual(proof["git_diff"]["status"], "PASS")
        for step in STEPS:
            block = proof["steps"][step]
            self.assertRegex(block["digest"], r"^[0-9a-f]{64}$")
            self.assertRegex(block["baseline_digest"], r"^[0-9a-f]{64}$")
            self.assertEqual(block["baseline_files"], block["files"])
        paths = json.dumps(receipt)
        self.assertNotIn('"/tmp', paths)
        self.assertNotIn("/home/chang", paths)


class HealthTest(AcceptanceTestCase):
    def test_health_runner_records_transcripts_and_runtimes(self):
        output = self.tmp / "health-run"
        value = acceptance.run_health(self.root, output, scripts=("check_paths",
                                                                  "check_tracked"))
        self.assertEqual(value["schema"], "tsg-view-health-v1")
        for name in ("check_paths.sh", "check_tracked.sh"):
            self.assertEqual(value["checks"][name]["returncode"], 0)
            self.assertTrue((output / (name + ".log")).is_file())
            self.assertGreaterEqual(value["checks"][name]["runtime_seconds"], 0.0)
            self.assertRegex(value["checks"][name]["log_sha256"], r"^[0-9a-f]{64}$")

    def test_failing_health_script_fails_closed(self):
        write(self.root / "script/check_paths.sh", "#!/bin/sh\nexit 1\n")
        with self.assertRaisesRegex(acceptance.AcceptanceError, "check_paths"):
            acceptance.run_health(self.root, self.tmp / "health-fail",
                                  scripts=("check_paths",))


class PlanTest(unittest.TestCase):
    def test_plan_origins_is_the_45_fit_schedule(self):
        self.assertEqual(acceptance.plan_origins(), ORIGINS)

    def test_origin_key_is_stable(self):
        self.assertEqual(acceptance.origin_key(ORIGINS[0]),
                         "tenv/full/fold-0/seed-42")


if __name__ == "__main__":
    unittest.main()
