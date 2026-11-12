"""One-command candidate -> acceptance -> publish for the #165 tenv matrix.

The declared matrix (ADR 0063) is 2 outer cold-drug folds x 2 training seeds of
the PRnet entrant plus the deterministic ``Mean`` reference.  Runs are built and
trained **serialized**, one heavy job at a time through ``execution/capped.sh``
under ``/usr/bin/time -v``, in a candidate tree outside the formal ``output/``.
Only after every run has reached ``complete`` and the shared #164 driver has
scored the whole set are the stable products published, the tracked plot-value
CSVs copied and every figure rendered and replay-verified.

    python execution/run_matrix.py --plan-only
    python execution/run_matrix.py --all
    python execution/run_matrix.py --all --resume

``--plan-only`` prints the declared identities and changes nothing.  A failed or
partial run aborts before any stable path is touched (fail closed).
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

STEP_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = STEP_ROOT / "execution" / "script"
EXPECTED_GENES = 18164
ENTRIES = ("prnet", "mean")
FOLDS, SEEDS = 2, 2
STABLE_PRODUCTS = ("inputs", "runs", "evaluation", "figure_data")
FIGURE_FILES = {"metrics": ("metrics_values.csv", "metrics.png"),
                "metrics_indirect": ("metrics_indirect_values.csv",
                                     "metrics_indirect.png"),
                "learning_curve": ("learning_curve_values.csv", "learning_curve.png")}
DEFAULT_BASES = Path("/home/chang/arc/stage2/01_bases/output")


class MatrixError(RuntimeError):
    """A matrix step failed; nothing is published."""


def sha256(path):
    import hashlib
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def repo_root():
    for parent in [STEP_ROOT, *STEP_ROOT.parents]:
        if (parent / ".git").exists():
            return parent
    raise MatrixError("cannot find the repository root")


def resolve_bases(explicit=None):
    """Prefer the checkout's payloads, else the canonical checkout's."""
    if explicit:
        return Path(explicit)
    local = repo_root() / "stage2/01_bases/output"
    if (local / "training-corpora/tenv.h5ad").is_file():
        return local
    return DEFAULT_BASES


def declared_identities(fold_product, folds=FOLDS, seeds=SEEDS):
    """Declared matrix cells, selected by row prefix through the fold/seed index."""
    sys.path.insert(0, str(SCRIPT_DIR))
    import mean_reference
    identities = []
    for fold_row in range(folds):
        for seed_row in range(seeds):
            identities.append(mean_reference.select_identities(
                fold_product, fold_row, seed_row, folds, seeds))
    labels = [f"fold_{cell['fold']}/seed_{cell['seed']}" for cell in identities]
    if len(set(labels)) != len(labels):
        raise MatrixError("declared matrix repeats a (fold, seed) identity")
    return identities


def role_limits(cohorts, fold):
    """Exact well limits of one fold, so prepare_inputs selects every well."""
    import pandas as pd
    scoped = cohorts[cohorts.fold.astype(int).eq(int(fold))]
    limits = {}
    for role in ("train", "validation", "test"):
        limits[role] = int(scoped.role.eq(role).sum())
        if not limits[role]:
            raise MatrixError(f"fold {fold} has no {role} wells")
    return limits


TIME_WALL = re.compile(r"Elapsed \(wall clock\) time.*: ([0-9:.]+(?: \([0-9:.]+\))?)")
TIME_RSS = re.compile(r"Maximum resident set size \(kbytes\): (\d+)")
TIME_USER = re.compile(r"User time \(seconds\): ([0-9.]+)")
TIME_SYSTEM = re.compile(r"System time \(seconds\): ([0-9.]+)")


def parse_time(path):
    """Parse the resource fields of one ``/usr/bin/time -v`` record."""
    text = Path(path).read_text(errors="replace")
    value = {"record": str(Path(path).resolve()), "config": "GNU time -v"}
    match = TIME_RSS.search(text)
    if not match:
        raise MatrixError(f"no maximum resident set size in {path}")
    value["peak_rss_kib"] = int(match.group(1))
    wall = TIME_WALL.search(text)
    if not wall:
        raise MatrixError(f"no wall clock time in {path}")
    token = wall.group(1).split(" ")[0]
    parts = [float(part) for part in token.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0.0)
    value["elapsed_seconds"] = round(parts[0] * 3600 + parts[1] * 60 + parts[2], 3)
    value["user_seconds"] = float(TIME_USER.search(text).group(1))
    value["system_seconds"] = float(TIME_SYSTEM.search(text).group(1))
    return value


def heavy(command, *, stdout_path, time_path, env=None, cwd=None):
    """Run one heavy command under the step's 10 GiB cap and GNU time."""
    stdout_path, time_path = Path(stdout_path), Path(time_path)
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    time_path.parent.mkdir(parents=True, exist_ok=True)
    wrapped = [str(STEP_ROOT / "execution/capped.sh"), "/usr/bin/time", "-v",
               *[str(item) for item in command]]
    environment = dict(os.environ)
    environment.update({"OMP_NUM_THREADS": "2", "MKL_NUM_THREADS": "2",
                        "OPENBLAS_NUM_THREADS": "2", "NUMEXPR_NUM_THREADS": "2",
                        "PYTHONPATH": str(SCRIPT_DIR)})
    environment.update(env or {})
    with open(stdout_path, "wb") as out, open(time_path, "wb") as err:
        completed = subprocess.run(wrapped, stdout=out, stderr=err,
                                   cwd=str(cwd or STEP_ROOT), env=environment)
    if completed.returncode != 0:
        raise MatrixError(
            f"command failed with exit {completed.returncode}: {' '.join(wrapped[:3])} "
            f"(see {stdout_path})")
    return parse_time(time_path)


def check_terminal(run_dir):
    """Fail closed: only an explicit ``complete`` terminal status counts."""
    terminal = read_json(Path(run_dir) / "terminal.json")
    status = terminal.get("status")
    if status != "complete":
        raise MatrixError(f"run {run_dir} terminal status is {status!r}, not complete")
    return terminal


def build_release(cell, *, bases, candidate, python):
    """Build one fold's complete predictor-input release (serialized)."""
    import pandas as pd
    cohorts = pd.read_csv(bases / "fold_n_seed/cohorts.csv", keep_default_na=False,
                          dtype={"sample_id": str, "role": str})
    limits = role_limits(cohorts, cell["fold"])
    release = candidate / "inputs" / f"fold_{cell['fold']}" / f"seed_{cell['seed']}"
    if release.is_dir() and read_json(release / "manifest.json")["identity"] == {
            "fold_row": cell["fold_row"], "fold": cell["fold"],
            "seed_row": cell["seed_row"], "seed": cell["seed"]}:
        return release, limits
    command = [python, "execution/script/prepare_inputs.py",
               "--corpus", str(bases / "training-corpora/tenv.h5ad"),
               "--fold-product", str(bases / "fold_n_seed"),
               "--output", str(release),
               "--fold-row", str(cell["fold_row"]), "--seed-row", str(cell["seed_row"]),
               "--folds", str(FOLDS), "--seeds", str(SEEDS),
               "--train-wells", str(limits["train"]),
               "--validation-wells", str(limits["validation"]),
               "--test-wells", str(limits["test"])]
    heavy(command, stdout_path=release / "prepare.out",
          time_path=release / "prepare.time")
    return release, limits


def train_prnet(cell, *, bases, candidate, config, python, log_root):
    """Train and predict one PRnet matrix cell (serialized, 10 GiB cap)."""
    release = candidate / "inputs" / f"fold_{cell['fold']}" / f"seed_{cell['seed']}"
    run_dir = candidate / "runs/prnet" / f"fold_{cell['fold']}" / f"seed_{cell['seed']}"
    if (run_dir / "terminal.json").is_file():
        try:
            check_terminal(run_dir)
            return run_dir
        except MatrixError:
            pass
    identity = {"model": "PRnet", **{key: cell[key] for key in
                                     ("fold_row", "fold", "seed_row", "seed")}}
    request_path = candidate / "requests" / f"fold_{cell['fold']}_seed_{cell['seed']}.json"
    write_json(request_path, {
        "identity": identity,
        "release": str(release),
        "n_genes": EXPECTED_GENES,
        "target_kind": "delta",
        "output_directory": str(run_dir),
    })
    environment = {"PYTHONHASHSEED": "0"}
    resources = heavy([python, "execution/script/prnet_fit_predict.py",
                       "--request", str(request_path), "--config", str(config)],
                      stdout_path=log_root / f"prnet-fold_{cell['fold']}-seed_{cell['seed']}.out",
                      time_path=log_root / f"prnet-fold_{cell['fold']}-seed_{cell['seed']}.time",
                      env=environment)
    terminal = check_terminal(run_dir)
    write_json(run_dir / "resources.json",
               {"resources": resources, "terminal": terminal})
    return run_dir


def build_mean(cell, *, bases, candidate, python, log_root):
    """Write the deterministic ``Mean`` run for one matrix cell (serialized)."""
    run_dir = candidate / "runs/mean" / f"fold_{cell['fold']}" / f"seed_{cell['seed']}"
    if (run_dir / "terminal.json").is_file():
        try:
            check_terminal(run_dir)
            return run_dir
        except MatrixError:
            pass
    resources = heavy([python, "execution/script/mean_reference.py",
                       "--corpus", str(bases / "training-corpora/tenv.h5ad"),
                       "--fold-product", str(bases / "fold_n_seed"),
                       "--output", str(run_dir),
                       "--fold-row", str(cell["fold_row"]),
                       "--seed-row", str(cell["seed_row"]),
                       "--folds", str(FOLDS), "--seeds", str(SEEDS)],
                      stdout_path=log_root / f"mean-fold_{cell['fold']}-seed_{cell['seed']}.out",
                      time_path=log_root / f"mean-fold_{cell['fold']}-seed_{cell['seed']}.time")
    check_terminal(run_dir)
    write_json(run_dir / "resources.json", {"resources": resources})
    return run_dir


def score(*, bases, candidate, python):
    """Score the whole run set with the shared #164 driver; fail on any omission."""
    output = candidate / "evaluation"
    figure_data = candidate / "figure_data"
    cells = declared_identities(bases / "fold_n_seed")
    command = [python, "execution/script/evaluation_driver.py",
               "--runs-root", str(candidate / "runs"),
               "--entries", ",".join(ENTRIES),
               "--fold-rows", ",".join(str(row) for row in sorted({cell["fold_row"] for cell in cells})),
               "--seeds", ",".join(str(seed) for seed in sorted({cell["seed"] for cell in cells})),
               "--fold-product", str(bases / "fold_n_seed"),
               "--truth", str(bases / "tenv-truth"),
               "--corpus", str(bases / "training-corpora/tenv.h5ad"),
               "--output", str(output), "--figure-data", str(figure_data)]
    completed = subprocess.run(command, cwd=str(STEP_ROOT), capture_output=True, text=True)
    if completed.returncode != 0:
        raise MatrixError(f"evaluation driver failed: {completed.stderr.strip()}")
    record = read_json(output / "evaluation.json")
    runs = record.get("runs") or []
    expected = len(ENTRIES) * FOLDS * SEEDS
    if len(runs) != expected:
        raise MatrixError(f"scored {len(runs)} runs, expected {expected}")
    not_complete = [run["label"] for run in runs
                    if run.get("terminal_status") != "complete"]
    if not_complete:
        raise MatrixError(f"unscored or incomplete runs: {not_complete}")
    tables = record.get("tables") or {}
    for name in ("per_well", "per_group", "per_run", "per_run_curve", "confusion"):
        if name not in tables or not Path(tables[name]["path"]).is_file():
            raise MatrixError(f"evaluation table is missing: {name}")
    if not (figure_data / "metrics_values.csv").is_file():
        raise MatrixError("metrics plot-value table is missing")
    if not (figure_data / "metrics_indirect_values.csv").is_file():
        raise MatrixError("reliability plot-value table is missing")
    if not (figure_data / "learning_curve_values.csv").is_file():
        raise MatrixError("learning-curve plot-value table is missing")
    return record


def verify_figure_values(figure_data):
    """Check the recorded tables are the shared ADR 0064 schema before publishing."""
    sys.path.insert(0, str(repo_root() / "stage2/01_bases/output/api"))
    from bases import plots
    metrics = plots.read_plot_values(Path(figure_data) / "metrics_values.csv")
    if metrics["schema"] != plots.METRICS_SCHEMA or len(metrics["panels"]) != 4:
        raise MatrixError("metrics plot values are not the four direct ADR 0064 panels")
    if (metrics["layout"]["rows"], metrics["layout"]["columns"]) != (1, 4):
        raise MatrixError("direct metrics plot values are not one row by four columns")
    if [panel["panel"] for panel in metrics["panels"]] != list(plots.PANEL_ORDER):
        raise MatrixError("metrics panel order differs from ADR 0064")
    if [panel["kind"] for panel in metrics["panels"]] != ["bar"] * 4:
        raise MatrixError("metrics panel kinds differ from ADR 0064")
    indirect = plots.read_plot_values(Path(figure_data) / "metrics_indirect_values.csv")
    if indirect["schema"] != plots.INDIRECT_SCHEMA or len(indirect["panels"]) != 8:
        raise MatrixError("reliability plot values are not the eight ADR 0064 panels")
    if (indirect["layout"]["rows"], indirect["layout"]["columns"]) != (2, 4):
        raise MatrixError("reliability plot values are not two rows by four columns")
    if [panel["panel"] for panel in indirect["panels"]] != list(plots.INDIRECT_PANEL_ORDER):
        raise MatrixError("reliability panel order differs from ADR 0064")
    expected_kinds = ["curve"] * 4 + ["bar", "bar", "heatmap", "bar"]
    if [panel["kind"] for panel in indirect["panels"]] != expected_kinds:
        raise MatrixError("reliability panel kinds differ from ADR 0064")
    curve = plots.read_plot_values(Path(figure_data) / "learning_curve_values.csv")
    if curve["schema"] != plots.LEARNING_CURVE_SCHEMA or len(curve["panels"]) != 1:
        raise MatrixError("learning curve must have one PRnet subplot")
    return metrics, indirect, curve


def publish_tree(source, destination):
    """Replace one stable product path only from a verified candidate tree."""
    source, destination = Path(source), Path(destination)
    if not source.is_dir():
        raise MatrixError(f"candidate product is missing: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = destination.with_name(destination.name + ".publishing")
    superseded = destination.with_name(destination.name + ".superseded")
    for path in (staging, superseded):
        if path.exists():
            shutil.rmtree(path)
    shutil.copytree(source, staging)
    if destination.exists():
        os.replace(destination, superseded)
    os.replace(staging, destination)
    if superseded.exists():
        shutil.rmtree(superseded)
    return destination


def render(figure_data, figure_dir, python):
    """Render every figure from the tracked CSVs alone; return the PNG paths."""
    command = [python, "figure/script/render_figures.py",
               "--figure-data", str(figure_data), "--out-dir", str(figure_dir)]
    completed = subprocess.run(command, cwd=str(STEP_ROOT), capture_output=True, text=True)
    if completed.returncode != 0:
        raise MatrixError(f"figure renderer failed: {completed.stderr.strip()}")
    paths = {name: Path(figure_dir) / png for name, (_, png) in FIGURE_FILES.items()}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise MatrixError(f"renderer did not produce: {missing}")
    return paths


def verify_replay(figure_data, figure_dir, python):
    """Prove every PNG rebuilds with the step's ``output/`` moved away."""
    output = STEP_ROOT / "output"
    hidden = STEP_ROOT / "output.replay-hidden"
    moved = False
    if output.exists():
        os.replace(output, hidden)
        moved = True
    try:
        render(figure_data, figure_dir, python)
    finally:
        if moved:
            os.replace(hidden, output)
    return {name: sha256(Path(figure_dir) / png) for name, (_, png) in FIGURE_FILES.items()}


def publish(candidate, *, figure_data, figure_dir, python):
    """Publish stable products, tracked plot values and every figure."""
    step_output = STEP_ROOT / "output"
    published = {name: str(publish_tree(candidate / name, step_output / name))
                 for name in STABLE_PRODUCTS}
    figure_data = Path(figure_data)
    figure_data.mkdir(parents=True, exist_ok=True)
    for name, (csv_name, _) in FIGURE_FILES.items():
        source = candidate / "figure_data" / csv_name
        if not source.is_file():
            raise MatrixError(f"candidate plot values are missing: {source}")
        shutil.copyfile(source, figure_data / csv_name)
    pngs = render(figure_data, figure_dir, python)
    digests = verify_replay(figure_data, figure_dir, python)
    return published, {name: str(path) for name, path in pngs.items()}, digests


def device_of(run_dir):
    """The recorded torch device of one PRnet run, wherever the identity nests it."""
    identity = read_json(Path(run_dir) / "identity.json")
    for block in (identity, identity.get("identity") or {},
                  identity.get("environment") or {}):
        if isinstance(block, dict) and block.get("device"):
            return str(block["device"])
    return "unknown"


def acceptance_record(*, cells, candidate, bases, evaluation, figures, digests,
                      environments, started_utc, finished_utc):
    """Assemble the compact acceptance record with hashes, statuses and peaks."""
    runs = []
    for cell in cells:
        for entry in ENTRIES:
            run_dir = candidate / "runs" / entry / f"fold_{cell['fold']}" / f"seed_{cell['seed']}"
            terminal = check_terminal(run_dir)
            resources = read_json(run_dir / "resources.json")["resources"]
            identity = read_json(run_dir / "identity.json")
            checkpoint_path = run_dir / "checkpoint.json"
            checkpoint = read_json(checkpoint_path) if checkpoint_path.is_file() else {}
            runs.append({
                "entry": entry, "label": identity["identity"].get("model"),
                "fold": cell["fold"], "seed": cell["seed"],
                "terminal_status": terminal["status"],
                "terminal_sha256": sha256(run_dir / "terminal.json"),
                "identity_sha256": sha256(run_dir / "identity.json"),
                "prediction_sha256": identity.get("prediction_sha256"),
                "case_ids_sha256": identity.get("case_ids_sha256"),
                "best_epoch": checkpoint.get("best_epoch"),
                "best_validation_loss": checkpoint.get("best_validation_loss"),
                "stop_reason": checkpoint.get("stop_reason"),
                "history_epochs": len(checkpoint.get("history") or []),
                "elapsed_seconds": resources["elapsed_seconds"],
                "peak_rss_kib": resources["peak_rss_kib"],
                "resources_record": resources["record"],
                "device": device_of(run_dir) if entry == "prnet" else "cpu",
            })
    released = []
    for cell in cells:
        release = candidate / "inputs" / f"fold_{cell['fold']}" / f"seed_{cell['seed']}"
        manifest = read_json(release / "manifest.json")
        released.append({"fold": cell["fold"], "seed": cell["seed"],
                         "identity": manifest["identity"],
                         "roles": manifest["roles"],
                         "manifest_sha256": sha256(release / "manifest.json"),
                         "files_sha256": manifest["sha256"],
                         "source": manifest["source"]})
    figure_hashes = {name: sha256(Path(STEP_ROOT / "figure") / png)
                     for name, (_, png) in FIGURE_FILES.items()}
    metric_rows = {}
    for name in ("per_well", "per_group", "per_run"):
        path = candidate / "evaluation" / f"{name}.csv"
        metric_rows[name] = {"path": f"output/evaluation/{name}.csv",
                             "sha256": sha256(path), "rows": sum(1 for _ in open(path)) - 1}
    return {
        "schema": "stage2-tenv-acceptance-v1",
        "ticket": 165,
        "created_utc": finished_utc,
        "started_utc": started_utc,
        "branch": subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"],
                                 cwd=STEP_ROOT, capture_output=True, text=True).stdout.strip(),
        "head_sha": subprocess.run(["git", "rev-parse", "HEAD"], cwd=STEP_ROOT,
                                   capture_output=True, text=True).stdout.strip(),
        "declared": {"entries": list(ENTRIES), "fold_rows": [0, 1], "seeds": [17, 42],
                     "matrix": "2 folds x 2 seeds + deterministic Mean reference",
                     "gene_count": EXPECTED_GENES},
        "environments": environments,
        "bases": {"root": str(bases),
                  "corpus_sha256": sha256(bases / "training-corpora/tenv.h5ad"),
                  "fold_product_manifest_sha256": sha256(bases / "fold_n_seed/manifest.json"),
                  "truth_manifest_sha256": sha256(bases / "tenv-truth/manifest.json")},
        "releases": released,
        "runs": runs,
        "scoring": {"record": "output/evaluation/evaluation.json",
                    "sha256": sha256(candidate / "evaluation/evaluation.json"),
                    "per_run_metrics": read_json(candidate / "evaluation/evaluation.json").get("runs")},
        "tables": metric_rows,
        "figures": {"data": {name: {"path": f"figure/data/{csv}",
                                    "sha256": sha256(STEP_ROOT / "figure/data" / csv)}
                             for name, (csv, _) in FIGURE_FILES.items()},
                    "png": {name: {"path": f"figure/{png}", "sha256": figure_hashes[name]}
                            for name, (_, png) in FIGURE_FILES.items()},
                    "replay": "every PNG rebuilt from figure/data alone with output/ moved away",
                    "replay_sha256": digests},
        "tolerances": {
            "stochastic": ("metric values differ across seeds; the four runs are a spread "
                           "indicator, not a confidence interval, and no run-to-run numeric "
                           "tolerance is asserted"),
            "figures": ("figure/data CSVs must equal the driver output byte for byte and each "
                        "PNG must be byte-identical to a fresh CSV-only replay"),
            "hashes": "sha256 over exact bytes",
        },
        "replay_commands": [
            "python figure/script/render_figures.py --figure-data figure/data --out-dir figure",
            "python execution/run_matrix.py --all",
        ],
        "deviations": [
            "no GPU on this boot (torch.cuda.is_available() false): CPU fallback per ADR 0063",
            "training/prediction runs use 2 CPU threads and the 10 GiB cap",
        ],
    }


def environments_for(*, python, prepare_python, train_python, config):
    def probe(executable, code):
        completed = subprocess.run([str(executable), "-c", code], capture_output=True, text=True)
        return completed.stdout.strip() or completed.stderr.strip()

    return {
        "orchestration": {"executable": str(python),
                          "version": probe(python, "import sys;print(sys.version.split()[0])")},
        "input_preparation": {"executable": str(prepare_python),
                              "version": probe(prepare_python,
                                               "import sys;print(sys.version.split()[0])")},
        "training": {"executable": str(train_python),
                     "version": probe(train_python,
                                      "import sys,torch,numpy;print(sys.version.split()[0], "
                                      "torch.__version__, numpy.__version__, "
                                      "torch.cuda.is_available())")},
        "train_config": {"path": str(config), "sha256": sha256(config)},
    }


def run_all(*, bases, candidate, python, prepare_python, train_python, config,
            figure_dir, resume):
    """Serialized matrix -> score -> verify -> publish.  Fail closed."""
    started = datetime.now(timezone.utc).isoformat()
    if candidate.exists() and not resume:
        shutil.rmtree(candidate)
    candidate.mkdir(parents=True, exist_ok=True)
    log_root = STEP_ROOT / "log" / "resources"
    cells = declared_identities(bases / "fold_n_seed")
    for cell in cells:
        release, limits = build_release(cell, bases=bases, candidate=candidate,
                                        python=prepare_python)
        print(f"inputs fold {cell['fold']} seed {cell['seed']}: {limits}", flush=True)
        train_prnet(cell, bases=bases, candidate=candidate, config=config,
                    python=train_python, log_root=log_root)
        print(f"prnet fold {cell['fold']} seed {cell['seed']}: complete", flush=True)
        build_mean(cell, bases=bases, candidate=candidate, python=python, log_root=log_root)
        print(f"mean  fold {cell['fold']} seed {cell['seed']}: complete", flush=True)
    evaluation = score(bases=bases, candidate=candidate, python=python)
    metrics, indirect, curve = verify_figure_values(candidate / "figure_data")
    environments = environments_for(python=python, prepare_python=prepare_python,
                                    train_python=train_python, config=config)
    # Publish only from the fully verified candidate (fail closed above).
    published, figure_paths, digests = publish(
        candidate, figure_data=STEP_ROOT / "figure/data",
        figure_dir=figure_dir, python=python)
    finished = datetime.now(timezone.utc).isoformat()
    record = acceptance_record(cells=cells, candidate=candidate, bases=bases,
                               evaluation=evaluation, figures=figure_paths,
                               digests=digests, environments=environments,
                               started_utc=started, finished_utc=finished)
    record["publish"] = {"products": published,
                         "figures": figure_paths,
                         "candidate": str(candidate)}
    write_json(STEP_ROOT / "log" / "165-acceptance.json", record)
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bases", default=None,
                        help="accepted 01_bases/output payload root")
    parser.add_argument("--candidate", default=str(
        STEP_ROOT / "execution/checkpoints/matrix-165"))
    parser.add_argument("--python", default="/home/chang/miniconda3/bin/python")
    parser.add_argument("--prepare-python",
                        default="/home/chang/miniconda3/envs/arc/bin/python")
    parser.add_argument("--train-python",
                        default="/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python")
    parser.add_argument("--train-config", default=str(
        STEP_ROOT / "execution/models/prnet/config.json"))
    parser.add_argument("--figure-dir", default=str(STEP_ROOT / "figure"))
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--resume", action="store_true",
                        help="reuse completed runs in the candidate tree")
    arguments = parser.parse_args(argv)
    bases = resolve_bases(arguments.bases)
    plan = {"bases": str(bases), "entries": list(ENTRIES),
            "cells": declared_identities(bases / "fold_n_seed")}
    if arguments.plan_only or not arguments.all:
        print(json.dumps(plan, indent=2, sort_keys=True))
        return 0
    try:
        record = run_all(bases=bases, candidate=Path(arguments.candidate),
                         python=arguments.python,
                         prepare_python=arguments.prepare_python,
                         train_python=arguments.train_python,
                         config=Path(arguments.train_config),
                         figure_dir=Path(arguments.figure_dir),
                         resume=arguments.resume)
    except MatrixError as error:
        print(f"matrix failed (nothing published): {error}", file=sys.stderr)
        return 2
    summary = {"head_sha": record["head_sha"],
               "runs": [{"entry": run["entry"], "fold": run["fold"],
                         "seed": run["seed"], "status": run["terminal_status"],
                         "elapsed_seconds": run["elapsed_seconds"],
                         "peak_rss_kib": run["peak_rss_kib"]} for run in record["runs"]],
               "figures": record["figures"]["png"]}
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
