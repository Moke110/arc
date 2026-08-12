#!/usr/bin/env python3
"""为 processed_data.h5 中方差为 0 的 x1/x2 行添加极小数值抖动，避免 Pearson 等指标未定义。

默认读取项目内 LINCS2020_set/processed_data.h5，写出到 --output（可指向 .cursor 或原路径）。

用法:
  python3 .cursor/fix_lincs2020_constant_profiles.py --dry-run
  python3 .cursor/fix_lincs2020_constant_profiles.py -o .cursor/processed_data_fixed.h5
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import h5py
import numpy as np


def _to_bytes_array(value: np.ndarray) -> np.ndarray:
    if value.size == 0:
        return value
    if np.issubdtype(value.dtype, np.str_):
        return value.astype(np.bytes_)
    if value.dtype.kind == "O":
        flat = value.reshape(-1)
        if all(isinstance(x, str) for x in flat):
            return np.array(flat, dtype=str).reshape(value.shape).astype(np.bytes_)
    return value


def save_h5(path: Path, data: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        for key, item in data.items():
            value = np.asarray(item)
            value = _to_bytes_array(value)
            f[key] = value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "training/data/LINCS2020_set/processed_data.h5",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "processed_data_fixed.h5",
        help="输出路径；与 input 相同时需配合 --in-place",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--eps",
        type=float,
        default=1e-7,
        help="抖动系数：对需修复行加上 eps * N(0,1)（逐基因）",
    )
    parser.add_argument("--dry-run", action="store_true", help="只打印将修改的行，不写文件")
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="将 input 备份为 .bak 后写入 output（若 output==input）",
    )
    args = parser.parse_args()

    inp = args.input
    out = args.output
    rng = np.random.default_rng(args.seed)

    with h5py.File(inp, "r") as f:
        data = {k: np.asarray(f[k]) for k in f}

    x1 = data["x1"].astype(np.float32, copy=True)
    x2 = data["x2"].astype(np.float32, copy=True)

    std1 = np.std(x1, axis=1)
    std2 = np.std(x2, axis=1)
    fix1 = std1 == 0
    fix2 = std2 == 0

    print("rows to fix (x1 zero std):", int(fix1.sum()), "indices:", np.where(fix1)[0].tolist())
    print("rows to fix (x2 zero std):", int(fix2.sum()), "indices:", np.where(fix2)[0].tolist())

    noise_scale = float(args.eps)
    for i in np.where(fix1)[0]:
        x1[i] = x1[i] + (noise_scale * rng.standard_normal(x1.shape[1], dtype=np.float64)).astype(np.float32)
    for i in np.where(fix2)[0]:
        x2[i] = x2[i] + (noise_scale * rng.standard_normal(x2.shape[1], dtype=np.float64)).astype(np.float32)

    assert not np.any(np.std(x1, axis=1) == 0), "x1 still has zero std rows"
    assert not np.any(np.std(x2, axis=1) == 0), "x2 still has zero std rows"

    data["x1"] = x1
    data["x2"] = x2

    if args.dry_run:
        print("dry-run: 未写入文件")
        return

    if args.in_place and out.resolve() == inp.resolve():
        bak = inp.with_suffix(inp.suffix + ".bak")
        print(f"备份: {inp} -> {bak}")
        shutil.copy2(inp, bak)

    print(f"写入: {out}")
    save_h5(out, data)


if __name__ == "__main__":
    main()
