#!/usr/bin/env python3
"""全面扫描 training/data/LINCS2020_set/processed_data.h5 的数值健康度。

用法:
  python3 .cursor/audit_lincs2020_h5.py
  python3 .cursor/audit_lincs2020_h5.py --h5 /path/to/processed_data.h5
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--h5",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "training/data/LINCS2020_set/processed_data.h5",
    )
    args = parser.parse_args()
    path = args.h5

    with h5py.File(path, "r") as f:
        keys = list(f.keys())
        print("datasets:", keys)
        for k in keys:
            d = f[k]
            print(f"  {k}: shape={d.shape} dtype={d.dtype}")

        x1 = np.asarray(f["x1"])
        x2 = np.asarray(f["x2"])

    n = x1.shape[0]
    g = x1.shape[1]

    print("\n=== 全局 ===")
    print("samples", n, "genes", g)
    for name, a in ("x1", x1), ("x2", x2):
        print(f"\n--- {name} ---")
        print("  nan count:", np.isnan(a).sum(), f"({100 * np.isnan(a).mean():.6f}%)")
        print("  inf count:", np.isinf(a).sum())
        finite = np.isfinite(a)
        print("  finite min/max:", float(np.min(a[finite])), float(np.max(a[finite])))

        row_std = np.std(a, axis=1)
        zv = row_std == 0
        print("  rows with std==0 (all finite):", int(zv.sum()))
        if zv.any():
            idx = np.where(zv)[0]
            print("  zero-variance row indices (first 50):", idx[:50].tolist())

        row_has_nan = np.isnan(a).any(axis=1)
        print("  rows with any nan:", int(row_has_nan.sum()))

    std1 = np.std(x1, axis=1)
    std2 = np.std(x2, axis=1)
    bad_x1 = std1 == 0
    bad_x2 = std2 == 0
    bad_either = bad_x1 | bad_x2

    print("\n=== 与 split 重叠（引用常数 profile 的样本数）===")
    split_dir = path.parent / "split"
    if split_dir.is_dir():
        for fname in ("train_id.h5", "valid_id.h5", "test_id.h5"):
            sp = split_dir / fname
            if not sp.is_file():
                continue
            with h5py.File(sp, "r") as sf:
                li = np.asarray(sf["LINCS_index"])
            hit_bad = np.isin(li, np.where(bad_either)[0]).sum()
            hit_x1 = np.isin(li, np.where(bad_x1)[0]).sum()
            print(f"  {fname}: n={len(li)}  hits(zero-std x1 or x2)={hit_bad}  hits(zero-std x1 only)={hit_x1}")

    import json

    payload = {
        "zero_std_x1": np.where(bad_x1)[0].astype(int).tolist(),
        "zero_std_x2": np.where(bad_x2)[0].astype(int).tolist(),
        "zero_std_either": np.where(bad_either)[0].astype(int).tolist(),
    }
    out_json = Path(__file__).resolve().parent / "LINCS2020_bad_row_indices.json"
    text = json.dumps(payload, indent=2)
    try:
        out_json.write_text(text, encoding="utf-8")
        print(f"\n已写入: {out_json}")
    except OSError as exc:
        print(f"\n无法写入 {out_json} ({exc!r})，JSON 如下：\n{text}")


if __name__ == "__main__":
    main()
