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

``--accept-only`` is the #170 no-training acceptance/replay path (ADR 0064).  It
never builds inputs, trains or re-infers: it verifies the eight already-published
runs and their identity digests, rescores them twice with the shared driver
(threads pinned to 2 so the BLAS reductions are deterministic), compares the two
trees byte for byte, republishes ``output/evaluation`` and ``output/figure_data``
only after that comparison, re-renders the three figures from the tracked
``figure/data`` and writes ``log/166-acceptance.json``.

    python execution/run_matrix.py --accept-only \
        --python /home/chang/miniconda3/envs/arc/bin/python
"""
import argparse
import csv
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
DEFAULT_OUTPUT = Path("/home/chang/arc/stage2/02_models_tenv/output")
ACCEPT_SCHEMA = "stage2-tenv-acceptance-v2"
ACCEPT_RECORD = STEP_ROOT / "log" / "166-acceptance.json"
ACCEPT_CHECKPOINT = STEP_ROOT / "execution/checkpoints/accept-170"
REPLAY_TABLES = ("per_well", "per_group", "per_run", "per_run_curve", "confusion")
# The accepted ADR 0064 figures and their complete plotted values.  A replay that
# changes any of these bytes is refused rather than silently rewriting them.
ACCEPTED_DATA = {
    "metrics": "c16b8c84e161214fa8e24567e4a8e0887a4715e11a05edae5fd37b13b9bdab88",
    "metrics_indirect":
        "6fdb0bc6520db686fc7a8e76101f75248061ec0c277939af9a235810c24dc980",
    "learning_curve":
        "4d054c03352a052b3aee8e2391c4ce8d01a8b75f8a941cac9e70963335eaff1c",
}
ACCEPTED_FIGURES = {
    "metrics": "3a58326edd5e31c9968a6f2f36e99c4b4d82d0ff8e131e9259e57046275f8ad1",
    "metrics_indirect":
        "e80d29054f1a49351fa832077b6d55e7b5d7f71f9cf84f49140ac325347bc3ba",
    "learning_curve":
        "f26d7b79049e0ecf653c6ee959ed1fa0e6b4793c1015422418681c35ca00f33f",
}


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


def resolve_output(explicit=None):
    """The stable step output: this checkout when it holds runs, else canonical."""
    if explicit:
        return Path(explicit).resolve()
    local = STEP_ROOT / "output"
    if (local / "runs").is_dir():
        return local
    return DEFAULT_OUTPUT


def _relative(path):
    try:
        return str(Path(path).resolve().relative_to(STEP_ROOT))
    except ValueError:
        return str(Path(path).resolve())


def verify_accepted_digests(paths, expected):
    """Refuse to replace accepted bytes with a different payload."""
    mismatches = []
    for name, path in paths.items():
        actual = sha256(path)
        if actual != expected[name]:
            mismatches.append(f"{name}: {actual} != accepted {expected[name]}")
    if mismatches:
        raise MatrixError("re-emitted payload does not match the accepted sha256; "
                          "refusing to rewrite: " + "; ".join(mismatches))


def verify_published_runs(output_root, cells):
    """Bind every declared run to its complete terminal record and digests."""
    rows = []
    for cell in cells:
        for entry in ENTRIES:
            run_dir = (Path(output_root) / "runs" / entry
                       / f"fold_{cell['fold']}" / f"seed_{cell['seed']}")
            if not run_dir.is_dir():
                raise MatrixError(f"published run is missing: {run_dir}")
            terminal = check_terminal(run_dir)
            identity = read_json(run_dir / "identity.json")
            for field, filename in (("prediction_sha256", "prediction.npy"),
                                    ("case_ids_sha256", "case_ids.json")):
                recorded = identity.get(field)
                if not recorded:
                    raise MatrixError(f"identity.json in {run_dir} is missing {field}")
                actual = sha256(run_dir / filename)
                if actual != recorded:
                    raise MatrixError(f"{field} digest differs in {run_dir}: "
                                      f"{actual} != {recorded}")
            rows.append({
                "entry": entry,
                "label": (identity.get("identity") or {}).get("model", entry),
                "fold": int(cell["fold"]), "seed": int(cell["seed"]),
                "terminal_status": terminal["status"],
                "terminal_sha256": sha256(run_dir / "terminal.json"),
                "identity_sha256": sha256(run_dir / "identity.json"),
                "prediction_sha256": identity.get("prediction_sha256"),
                "case_ids_sha256": identity.get("case_ids_sha256"),
            })
    return rows


def driver_command(*, bases, runs_root, output, figure_data, python):
    """The shared #168 evaluator over the already-published predictions."""
    cells = declared_identities(Path(bases) / "fold_n_seed")
    return [str(python), "execution/script/evaluation_driver.py",
            "--runs-root", str(runs_root),
            "--entries", ",".join(ENTRIES),
            "--fold-rows", ",".join(str(row) for row in
                                    sorted({cell["fold_row"] for cell in cells})),
            "--seeds", ",".join(str(seed) for seed in
                                sorted({cell["seed"] for cell in cells})),
            "--fold-product", str(Path(bases) / "fold_n_seed"),
            "--truth", str(Path(bases) / "tenv-truth"),
            "--corpus", str(Path(bases) / "training-corpora/tenv.h5ad"),
            "--output", str(output), "--figure-data", str(figure_data)]


def render_command(*, figure_data, out_dir, python):
    return [str(python), "figure/script/render_figures.py",
            "--figure-data", str(figure_data), "--out-dir", str(out_dir)]


def accept_commands(*, bases, output_root, python, candidate=None, figure_dir=None):
    """The only commands the acceptance path may run: rescore twice and render."""
    candidate = Path(candidate or ACCEPT_CHECKPOINT / "candidate")
    output_root = Path(output_root)
    figure_dir = Path(figure_dir or STEP_ROOT / "figure")
    return [
        driver_command(bases=bases, runs_root=output_root / "runs",
                       output=candidate / "evaluation",
                       figure_data=candidate / "figure_data", python=python),
        driver_command(bases=bases, runs_root=output_root / "runs",
                       output=output_root / "evaluation",
                       figure_data=output_root / "figure_data", python=python),
        render_command(figure_data=STEP_ROOT / "figure/data", out_dir=figure_dir,
                       python=python),
    ]


def compare_product_trees(first, second, names):
    """Prove a replay tree reproduces the published table bytes exactly."""
    mismatches = []
    for name in names:
        filename = name if Path(name).suffix else f"{name}.csv"
        left, right = Path(first) / filename, Path(second) / filename
        if not left.is_file() or not right.is_file():
            mismatches.append(f"{filename}: missing")
        elif sha256(left) != sha256(right):
            mismatches.append(f"{filename}: {sha256(left)} != {sha256(right)}")
    if mismatches:
        raise MatrixError("replay mismatch between the candidate and the published "
                          "products: " + "; ".join(mismatches))


def _normalised_record(record, root):
    """An evaluation record with its absolute product root removed."""
    text = json.dumps(record, sort_keys=True).replace(str(Path(root).resolve()), "<root>")
    return json.loads(text)


def compare_evaluation_records(first, first_root, second, second_root):
    if _normalised_record(first, first_root) != _normalised_record(second, second_root):
        raise MatrixError("replay mismatch between the candidate and the published "
                          "evaluation.json (values other than the product paths differ)")


def _restore_product(output_root, backup, names):
    for name in names:
        source = Path(backup) / name
        destination = Path(output_root) / name
        if source.exists():
            if destination.exists():
                shutil.rmtree(destination)
            shutil.copytree(source, destination)


def ensure_accepted_file(source, destination, expected, label):
    """Publish a tracked value only when it already equals the accepted bytes."""
    if sha256(source) != expected:
        raise MatrixError(f"re-emitted {label} does not match the accepted sha256; "
                          "refusing to rewrite the tracked artifact")
    destination = Path(destination)
    if destination.is_file():
        if sha256(destination) != expected:
            raise MatrixError(f"tracked accepted {label} was modified; refusing to "
                              "rewrite it")
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    if sha256(destination) != expected:
        raise MatrixError(f"published tracked {label} does not match the accepted sha256")
    return destination


def environment_record(python, commands):
    """Record the machine, interpreter, packages and GPU state of the replay."""
    probe = (
        "import json,platform,sys;"
        "import numpy,pandas,matplotlib;"
        "print(json.dumps({'python':sys.version.split()[0],"
        "'executable':sys.executable,'platform':platform.platform(),"
        "'release':platform.release(),'version':platform.version(),"
        "'machine':platform.machine(),'numpy':numpy.__version__,"
        "'pandas':pandas.__version__,'matplotlib':matplotlib.__version__,"
        "'wsl':('microsoft' in platform.release().lower() "
        "or 'microsoft' in platform.version().lower())}))")
    completed = subprocess.run([str(python), "-c", probe],
                               capture_output=True, text=True)
    if completed.returncode != 0:
        raise MatrixError(f"environment probe failed: {completed.stderr.strip()}")
    packages = json.loads(completed.stdout.strip().splitlines()[-1])
    try:
        gpu = subprocess.run(["nvidia-smi"], capture_output=True, text=True)
        gpu_record = {"command": "nvidia-smi", "available": gpu.returncode == 0,
                      "returncode": gpu.returncode,
                      "stdout": gpu.stdout.strip(), "stderr": gpu.stderr.strip()}
    except FileNotFoundError as error:
        gpu_record = {"command": "nvidia-smi", "available": False,
                      "returncode": None, "stdout": "", "stderr": str(error)}
    git_tool = lambda *args: subprocess.run(  # noqa: E731
        ["git", *args], cwd=STEP_ROOT, capture_output=True, text=True).stdout.strip()
    return {
        "git": {"branch": git_tool("rev-parse", "--abbrev-ref", "HEAD"),
                "head_sha": git_tool("rev-parse", "HEAD"),
                "status_porcelain": git_tool("status", "--porcelain")},
        "machine": packages,
        "gpu": gpu_record,
        "commands": [" ".join(str(item) for item in command) for command in commands],
    }


def acceptance_record_v2(*, cells, published, output_root, bases, candidate,
                         resources, environment, started_utc, finished_utc,
                         commands):
    """Assemble the self-describing ADR 0064 acceptance record."""
    output_root = Path(output_root)
    resources = {name: {**value, "record": _relative(value["record"])}
                 for name, value in resources.items()}
    evaluation = read_json(output_root / "evaluation/evaluation.json")
    rows = list(csv.DictReader(open(output_root / "evaluation/per_run.csv")))
    by_key = {(row["entry"], int(row["fold"]), int(row["seed"])): row for row in rows}
    runs = []
    for run in published:
        row = by_key[(run["entry"], run["fold"], run["seed"])]
        runs.append({
            **run,
            "label": row["label"],
            "wells": int(row["wells"]), "groups": int(row["groups"]),
            "null_groups": int(row["null_groups"]),
            "changed_groups": int(row["changed_groups"]),
            "M": float(row["M"]), "d_star": float(row["d_star"]),
            "d_star_k": int(row["d_star_k"]),
            "d_star_prime": float(row["d_star_prime"]),
            "Precision@d*": float(row["Precision@d*"]),
            "Recall@d*": float(row["Recall@d*"]),
            "F1@d*": float(row["F1@d*"]),
            "Spearman@d*": float(row["Spearman@d*"]),
            "Reliability@d*": float(row["Reliability@d*"]),
            "null_Reliability@d*": float(row["null_Reliability@d*"]),
            "changed_Reliability@d*": float(row["changed_Reliability@d*"]),
            "measured_effects": int(row["measured_effects"]),
        })
    tables = {name: {"path": f"output/evaluation/{name}.csv",
                     "sha256": sha256(output_root / "evaluation" / f"{name}.csv"),
                     "rows": len(list(csv.DictReader(
                         open(output_root / "evaluation" / f"{name}.csv"))))}
              for name in (*REPLAY_TABLES,)}
    evaluation_record = {
        "path": "output/evaluation/evaluation.json",
        "sha256": sha256(output_root / "evaluation/evaluation.json"),
        "schema": evaluation.get("schema"),
        "threshold_grid": evaluation.get("threshold_grid"),
        "direct_metrics": evaluation.get("direct_metrics"),
        "reliability_metrics": evaluation.get("reliability_metrics"),
    }
    figure_data = {name: {"path": f"figure/data/{csv_name}",
                          "sha256": sha256(STEP_ROOT / "figure/data" / csv_name)}
                   for name, (csv_name, _) in FIGURE_FILES.items()}
    figure_png = {name: {"path": f"figure/{png}",
                         "sha256": sha256(STEP_ROOT / "figure" / png)}
                  for name, (_, png) in FIGURE_FILES.items()}
    candidate = Path(candidate)
    candidate_tables = {
        name: sha256(candidate / "evaluation" / f"{name}.csv") for name in REPLAY_TABLES}
    published_tables = {
        name: sha256(output_root / "evaluation" / f"{name}.csv") for name in REPLAY_TABLES}
    return {
        "schema": ACCEPT_SCHEMA,
        "ticket": 170,
        "supersedes": "log/165-acceptance.json",
        "started_utc": started_utc,
        "created_utc": finished_utc,
        "branch": environment["git"]["branch"],
        "head_sha": environment["git"]["head_sha"],
        "environment": environment,
        "declared": {
            "entries": list(evaluation.get("entries") or ENTRIES),
            "fold_rows": list(evaluation.get("fold_rows") or [0, 1]),
            "seeds": list(evaluation.get("seeds") or [17, 42]),
            "matrix": "2 folds x 2 seeds + deterministic Mean reference",
            "gene_count": EXPECTED_GENES,
            "threshold_grid": evaluation.get("threshold_grid"),
        },
        "bases": {
            "root": str(bases),
            "corpus_sha256": sha256(Path(bases) / "training-corpora/tenv.h5ad"),
            "fold_product_manifest_sha256":
                sha256(Path(bases) / "fold_n_seed/manifest.json"),
            "truth_manifest_sha256": sha256(Path(bases) / "tenv-truth/manifest.json"),
        },
        "runs": runs,
        "evaluation": evaluation_record,
        "tables": tables,
        "figures": {
            "data": figure_data,
            "png": figure_png,
            "replay": ("each PNG was re-rendered from the tracked figure/data alone "
                       "and compared to the accepted sha256 before being installed"),
            "accepted_sha256": {"data": ACCEPTED_DATA, "png": ACCEPTED_FIGURES},
        },
        "resources": resources,
        "replay": {
            "method": ("the shared driver was executed twice with OMP, MKL, OPENBLAS "
                       "and NUMEXPR pinned to 2 threads: once into a candidate tree "
                       "and once into the stable output; every evaluation table, "
                       "evaluation record and figure-data CSV was compared byte for "
                       "byte, and the three figures were re-rendered from the tracked "
                       "figure/data alone"),
            "candidate": _relative(candidate),
            "candidate_tables": candidate_tables,
            "published_tables": published_tables,
            "tables_byte_identical": True,
            "evaluation_record_equal_except_paths": True,
            "figure_data_byte_identical_to_accepted": True,
            "figures_byte_identical_to_accepted": True,
            "per_run_values_reproduced": True,
        },
        "tolerances": {
            "rescore": ("deterministic: the driver is run twice and every evaluation "
                        "table and plot-value CSV must match byte for byte; no "
                        "stochastic term exists in the rescoring path"),
            "per_run_values": "exact equality; no numeric tolerance is applied",
            "figures": ("each PNG must equal its accepted sha256 and is a function of "
                        "the tracked figure/data alone"),
            "hashes": "sha256 over exact bytes",
        },
        "replay_commands": [
            "python execution/run_matrix.py --accept-only "
            "--python /home/chang/miniconda3/envs/arc/bin/python",
            "python figure/script/render_figures.py "
            "--figure-data figure/data --out-dir figure",
            *[" ".join(str(item) for item in command) for command in commands],
        ],
        "deviations": [
            ("the declared 2-fold matrix scores 918 test groups (fold 0: 459 with 90 "
             "null results; fold 1: 459 with 79 null results), not the corpus-wide "
             "1,377 parent-dose groups (1,104 non-empty, 273 null). Fold 2 has no runs "
             "and re-inference is out of scope: accepted, documented deviation"),
            ("PRnet's peak Reliability@d is at the strictest grid point k=100 (d'=1) "
             "in all four runs; at that point the mean window is 0.03 genes, so the "
             "peak is a nearly empty-table solution whose value is the empty-window "
             "null-group share (fold 0 90/459 = 0.196078; fold 1 76/459 = 0.165577, "
             "with 76 of the 79 null groups empty at k=100)"),
            ("delivered table size W(d*) is a property of each run and its d*: PRnet's "
             "delivered table is near-empty, so downstream phenotype work must treat "
             "table size as part of the intervention rather than as noise"),
        ],
        "notes": [
            ("peak selection is optimistic: the same group data selects d* and reports "
             "the value; the full curve and the four per-run d* values are recorded"),
            ("all four PRnet runs record M=5.926281 because the same cold-drug test "
             "wells define the predicted magnitude scale; this is a genuine data "
             "property, not a copied value"),
            ("the two Mean entries per fold are identical seed copies of the "
             "deterministic reference, not independent fits"),
        ],
    }


def mark_superseded(path, pointer):
    """Add the forward pointer without rewriting the historical record's layout."""
    path = Path(path)
    if not path.is_file():
        return
    text = path.read_text()
    record = json.loads(text)
    if record.get("superseded_by") == pointer:
        return
    stripped = text.rstrip()
    if not stripped.endswith("}"):
        raise MatrixError(f"cannot add a superseded_by pointer to {path}")
    nested = json.dumps(pointer, indent=2, sort_keys=True).replace("\n", "\n  ")
    path.write_text(f"{stripped[:-1].rstrip()},\n  \"superseded_by\": {nested}\n}}\n")
    if json.loads(path.read_text()).get("superseded_by") != pointer:
        raise MatrixError(f"superseded_by pointer did not round-trip in {path}")


def run_accept(*, bases, output_root, python, figure_dir, checkpoint_root=None,
               record_path=None):
    """Rescore and replay the accepted matrix without training or re-inference."""
    started = datetime.now(timezone.utc).isoformat()
    output_root = Path(output_root).resolve()
    figure_dir = Path(figure_dir).resolve()
    cells = declared_identities(Path(bases) / "fold_n_seed")
    published = verify_published_runs(output_root, cells)
    checkpoint = Path(checkpoint_root or ACCEPT_CHECKPOINT)
    candidate = checkpoint / "candidate"
    replay_figures = checkpoint / "replay-figures"
    for path in (candidate, replay_figures):
        if path.exists():
            shutil.rmtree(path)
    candidate.mkdir(parents=True, exist_ok=True)
    replay_figures.mkdir(parents=True, exist_ok=True)
    commands = accept_commands(bases=bases, output_root=output_root, python=python,
                               candidate=candidate, figure_dir=figure_dir)
    resources = {}
    resources["rescore_verification"] = heavy(
        commands[0], stdout_path=checkpoint / "rescore-verification.out",
        time_path=checkpoint / "rescore-verification.time")
    verify_figure_values(candidate / "figure_data")
    verify_accepted_digests(
        {name: candidate / "figure_data" / csv_name
         for name, (csv_name, _) in FIGURE_FILES.items()}, ACCEPTED_DATA)
    backup = checkpoint / "previous"
    if backup.exists():
        shutil.rmtree(backup)
    backup.mkdir(parents=True)
    for name in ("evaluation", "figure_data"):
        if (output_root / name).is_dir():
            shutil.copytree(output_root / name, backup / name)
    try:
        resources["rescore_publish"] = heavy(
            commands[1], stdout_path=checkpoint / "rescore-publish.out",
            time_path=checkpoint / "rescore-publish.time")
        compare_product_trees(candidate / "evaluation", output_root / "evaluation",
                              REPLAY_TABLES)
        compare_evaluation_records(
            read_json(candidate / "evaluation/evaluation.json"), candidate,
            read_json(output_root / "evaluation/evaluation.json"), output_root)
        verify_accepted_digests(
            {name: output_root / "figure_data" / csv_name
             for name, (csv_name, _) in FIGURE_FILES.items()}, ACCEPTED_DATA)
    except BaseException:
        _restore_product(output_root, backup, ("evaluation", "figure_data"))
        raise
    for name, (csv_name, _) in FIGURE_FILES.items():
        ensure_accepted_file(candidate / "figure_data" / csv_name,
                             STEP_ROOT / "figure/data" / csv_name,
                             ACCEPTED_DATA[name], f"plot-value table {name}")
    replay_command = render_command(figure_data=STEP_ROOT / "figure/data",
                                    out_dir=replay_figures, python=python)
    resources["figure_replay"] = heavy(
        replay_command, stdout_path=checkpoint / "figure-replay.out",
        time_path=checkpoint / "figure-replay.time")
    verify_accepted_digests(
        {name: replay_figures / png_name
         for name, (_, png_name) in FIGURE_FILES.items()}, ACCEPTED_FIGURES)
    for name, (_, png_name) in FIGURE_FILES.items():
        shutil.copyfile(replay_figures / png_name, Path(figure_dir) / png_name)
    finished = datetime.now(timezone.utc).isoformat()
    record = acceptance_record_v2(
        cells=cells, published=published, output_root=output_root, bases=bases,
        candidate=candidate, resources=resources,
        environment=environment_record(python, [commands[0], commands[1],
                                                replay_command]),
        started_utc=started, finished_utc=finished, commands=commands)
    record["publish"] = {"output_root": _relative(output_root),
                         "evaluation": _relative(output_root / "evaluation"),
                         "figure_data": _relative(output_root / "figure_data"),
                         "figure_dir": _relative(figure_dir)}
    target = Path(record_path or ACCEPT_RECORD)
    write_json(target, record)
    mark_superseded(STEP_ROOT / "log/165-acceptance.json",
                    {"record": "log/166-acceptance.json",
                     "schema": ACCEPT_SCHEMA, "ticket": 170})
    return record


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
    parser.add_argument("--accept-only", action="store_true",
                        help="rescore/replay the published runs; never train or re-infer")
    parser.add_argument("--output-root", default=None,
                        help="stable step output; defaults to this checkout, else canonical")
    parser.add_argument("--record", default=str(ACCEPT_RECORD),
                        help="acceptance record path for --accept-only")
    arguments = parser.parse_args(argv)
    if arguments.accept_only and (arguments.all or arguments.plan_only
                                  or arguments.resume):
        print("--accept-only cannot be combined with --all, --plan-only or --resume",
              file=sys.stderr)
        return 2
    bases = resolve_bases(arguments.bases)
    output_root = resolve_output(arguments.output_root)
    plan = {"bases": str(bases), "output": str(output_root), "entries": list(ENTRIES),
            "cells": declared_identities(bases / "fold_n_seed")}
    if arguments.accept_only:
        try:
            record = run_accept(bases=bases, output_root=output_root,
                                python=arguments.python,
                                figure_dir=Path(arguments.figure_dir),
                                record_path=Path(arguments.record))
        except MatrixError as error:
            print(f"acceptance failed (accepted products untouched): {error}",
                  file=sys.stderr)
            return 2
        summary = {"head_sha": record["head_sha"], "schema": record["schema"],
                   "record": str(Path(arguments.record).resolve()),
                   "runs": [{"entry": run["entry"], "fold": run["fold"],
                             "seed": run["seed"], "M": run["M"],
                             "d_star_k": run["d_star_k"],
                             "Reliability@d*": run["Reliability@d*"]}
                            for run in record["runs"]],
                   "figures": {name: value["sha256"]
                               for name, value in record["figures"]["png"].items()}}
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0
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
