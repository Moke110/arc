"""
Per-env model runner harness (ticket 06).

Drives a model adapter across an environment boundary: serializes a fold's
`(train_cases, test_cases)` to disk, runs the model's entry script inside its own
pinned conda env (`arc-stage2-v1-<model>`) as a subprocess, and reads back the
predicted Δ matrix aligned to the common gene universe. The workflow never
imports model internals; all model access is through the `predict` seam.

A model run that fails (non-zero exit, missing or mis-shaped output) raises
loudly and produces no score. The harness leaves the shared `arc` env untouched.
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from adapter import Cases  # noqa: E402


def _save_cases(path, cases):
    arrays = {
        "baseline": cases.baseline,
        "dose_log": cases.dose_log,
    }
    if cases.delta is not None:
        arrays["delta"] = cases.delta
    np.savez_compressed(path, **arrays)
    meta = {
        "drug": list(cases.drug),
        "dataset": list(cases.dataset),
        "smiles": list(cases.smiles) if cases.smiles is not None else [None] * len(cases),
        "group_key": [list(k) for k in cases.group_key],
        "n": len(cases),
        "n_genes": cases.n_genes,
    }
    path.with_suffix(".json").write_text(json.dumps(meta))


def _load_cases(path):
    z = np.load(path)
    meta = json.loads(path.with_suffix(".json").read_text())
    return Cases(
        baseline=z["baseline"], drug=meta["drug"], dataset=meta["dataset"],
        dose_log=z["dose_log"], group_key=[tuple(k) for k in meta["group_key"]],
        smiles=meta["smiles"], delta=z["delta"] if "delta" in z else None)


def run_model(env, entry_script, train, test, timeout=3600, extra_args=None):
    """Run `entry_script` in conda `env`; return the predicted Δ (n_test, n_genes)."""
    entry = Path(entry_script).resolve()
    if not entry.exists():
        raise FileNotFoundError(f"model entry script not found: {entry}")
    with tempfile.TemporaryDirectory(prefix="cm_bench_") as td:
        td = Path(td)
        _save_cases(td / "train.npz", train)
        _save_cases(td / "test.npz", test)
        out_path = td / "pred.npz"
        cmd = ["conda", "run", "--no-capture-output", "-n", env,
               "python", str(entry), str(td / "train.npz"), str(td / "test.npz"),
               str(out_path)]
        if extra_args:
            cmd += list(extra_args)
        print(f"[runner] ({env}) {' '.join(cmd[3:])}", flush=True)
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if proc.returncode != 0:
            raise RuntimeError(
                f"model '{entry.name}' in env '{env}' exited {proc.returncode}\n"
                f"--- stdout ---\n{proc.stdout[-4000:]}\n"
                f"--- stderr ---\n{proc.stderr[-4000:]}")
        if not out_path.exists():
            raise RuntimeError(
                f"model '{entry.name}' produced no prediction file "
                f"(stdout tail: {proc.stdout[-500:]})")
        pred = np.load(out_path)["pred"]
    expected = (len(test), test.n_genes)
    if pred.shape != expected:
        raise RuntimeError(
            f"model '{entry.name}' returned shape {pred.shape}, expected {expected}")
    return pred


class SubprocessAdapter:
    """Adapter that delegates `predict` to an entry script in a pinned env."""

    supports_dose = True

    def __init__(self, name, env, entry_script, extra_args=None):
        self.name = name
        self.env = env
        self.entry_script = entry_script
        self.extra_args = extra_args

    def predict(self, train, test):
        return run_model(self.env, self.entry_script, train, test,
                         extra_args=self.extra_args)
