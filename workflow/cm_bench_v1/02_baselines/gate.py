#!/usr/bin/env python3
"""
In-sample sanity gate for the baselines (ticket 05).

`ridge-fp` and `nearest-drug` must beat `no-change` when trained and evaluated
on the *same* wells (sanity only — not a benchmark claim). This catches broken
fingerprint/dose plumbing before the cold-drug run.

    python 02_baselines/gate.py
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_lib"))
sys.path.insert(0, str(ROOT / "02_baselines"))
from adapter import Cases  # noqa: E402
from baselines import (MeanEffectBaseline, NearestDrugBaseline,  # noqa: E402
                       NoChangeAdapter, RidgeFingerprintBaseline)
from bench_io import load_bench, load_true_sets  # noqa: E402
from metric import direction_aware_auprc  # noqa: E402


def _macro(adapter, meta, baseline, delta, groups, labels):
    sub = meta[meta["group_key"].isin(set(groups.loc[groups["reliable"],
                                                      "group_key"]))]
    mask = meta["group_key"].isin(set(sub["group_key"])).to_numpy()
    train = Cases.from_frame(meta, baseline, delta)
    test = Cases.from_frame(meta[mask], baseline[mask])
    pred = adapter.predict(train, test)
    lab_of = {k: labels[i] for i, k in enumerate(groups["group_key"])}
    vals = []
    for key, idx in sub.groupby("group_key").indices.items():
        p = pred[np.asarray(idx)].mean(axis=0)
        a = direction_aware_auprc(p, lab_of[key])
        if np.isfinite(a):
            vals.append(a)
    return float(np.mean(vals))


def main():
    bench = load_bench()
    groups, labels = load_true_sets()
    meta, baseline, delta = bench["meta"], bench["baseline"], bench["delta"]
    scores = {
        "no-change": _macro(NoChangeAdapter(), meta, baseline, delta, groups, labels),
        "mean-effect": _macro(MeanEffectBaseline(), meta, baseline, delta, groups, labels),
        "ridge-fp": _macro(RidgeFingerprintBaseline(), meta, baseline, delta, groups, labels),
        "nearest-drug": _macro(NearestDrugBaseline(), meta, baseline, delta, groups, labels),
    }
    for k, v in scores.items():
        print(f"[gate] in-sample macro AUPRC {k}: {v:.4f}")
    assert scores["ridge-fp"] > scores["no-change"], "ridge-fp does not beat no-change"
    assert scores["nearest-drug"] > scores["no-change"], \
        "nearest-drug does not beat no-change"
    print("02_baselines/gate: OK")


if __name__ == "__main__":
    main()
