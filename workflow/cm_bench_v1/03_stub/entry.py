#!/usr/bin/env python3
"""
Stub model entry script (runs inside the model's own conda env).

Exercises the runner harness end-to-end without any real model: it predicts the
train-mean Δ for every test sample. It imports only numpy so it can run in a
minimal throwaway environment.

    python entry.py <train.npz> <test.npz> <out.npz>
"""

import json
import sys

import numpy as np


def main():
    train_path, test_path, out_path = sys.argv[1:4]
    train = np.load(train_path)
    test = np.load(test_path)
    n_test = test["baseline"].shape[0]
    n_genes = test["baseline"].shape[1]
    if "delta" in train and train["delta"].shape[0] > 0:
        pred = np.tile(train["delta"].mean(axis=0), (n_test, 1))
    else:
        pred = np.zeros((n_test, n_genes))
    # sanity: the meta sidecar must be readable (as a real model would need it)
    _ = json.loads(open(test_path[:-4] + ".json").read())
    np.savez_compressed(out_path, pred=pred)


if __name__ == "__main__":
    main()
