#!/usr/bin/env python3
"""
Gate for the model runner harness (ticket 06).

Proves the per-env subprocess bridge on CPU with the stub model, and that a
failing model surfaces loudly instead of silently producing a score.

Run:
    python 03_stub/gate.py [--env arc-stage2-v1-stub]
"""

import argparse
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_lib"))
from adapter import Cases  # noqa: E402
from model_runner import run_model  # noqa: E402

ENTRY = ROOT / "03_stub" / "entry.py"


def _cases(n, g, seed):
    rng = np.random.default_rng(seed)
    base = rng.standard_normal((n, g))
    delta = rng.standard_normal((n, g))
    return Cases(baseline=base, drug=["A"] * n, dataset=["dtoxs_primary"] * n,
                 dose_log=np.zeros(n), group_key=[("dtoxs_primary", "A", 0)] * n,
                 smiles=["CCO"] * n, delta=delta)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--env", default="arc-stage2-v1-stub")
    args = ap.parse_args()

    train, test = _cases(6, 20, 1), _cases(3, 20, 2)
    pred = run_model(args.env, ENTRY, train, test)
    assert pred.shape == (3, 20), f"unexpected stub shape {pred.shape}"
    assert np.allclose(pred, train.delta.mean(axis=0)), "stub prediction wrong"
    print("[gate] stub round-trip OK (shape + values)")

    # a failing model must raise, not return a score
    with tempfile.TemporaryDirectory() as td:
        bad = Path(td) / "bad_entry.py"
        bad.write_text("import sys\nsys.exit(7)\n")
        try:
            run_model(args.env, bad, train, test)
        except RuntimeError as exc:
            assert "exited 7" in str(exc), str(exc)
            print("[gate] failing model surfaced loudly OK")
        else:
            raise AssertionError("failing model did not raise")

    # a mis-shaped output must raise
    with tempfile.TemporaryDirectory() as td:
        bad = Path(td) / "bad_shape.py"
        bad.write_text(
            "import sys, numpy as np\n"
            "np.savez_compressed(sys.argv[3], pred=np.zeros((2, 2)))\n")
        try:
            run_model(args.env, bad, train, test)
        except RuntimeError as exc:
            assert "expected" in str(exc), str(exc)
            print("[gate] mis-shaped output surfaced loudly OK")
        else:
            raise AssertionError("mis-shaped output did not raise")

    print("03_stub/gate: OK")


if __name__ == "__main__":
    main()
