#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import pickle
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from rdkit import Chem

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MODEL_PATH = REPO_ROOT / "models" / "pretrained" / "base" / "base.pth"
DEFAULT_CACHE_ROOT = REPO_ROOT / "feature_extract" / ".cache"
DEFAULT_DATASET_NAME = "kpgt_tmp_dataset"


def _smiles_ok(s: str) -> bool:
    return Chem.MolFromSmiles(s) is not None


def _env() -> dict:
    env = os.environ.copy()
    env.setdefault("DGLBACKEND", "pytorch")
    return env


def _run(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=str(cwd), check=True, env=_env())


def prepare_downstream_csv(smiles_csv: Path, dataset_dir: Path, max_molecules: Optional[int]) -> Path:
    df = pd.read_csv(smiles_csv)
    if "canonical_smiles" not in df.columns:
        raise ValueError(f"CSV 必须包含 canonical_smiles 列，当前列: {list(df.columns)}")
    raw = df["canonical_smiles"].astype(str).str.strip()
    if max_molecules is not None:
        raw = raw.iloc[:max_molecules]

    mask = raw.map(_smiles_ok)
    invalid_count = int((~mask).sum())
    if invalid_count:
        print(f"[warn] 丢弃无法被 RDKit 解析的 SMILES: {invalid_count} 条", file=sys.stderr)

    valid_smiles = raw[mask].reset_index(drop=True)
    out_df = pd.DataFrame(
        {
            "smiles": valid_smiles,
            "_kpgt_dummy": np.zeros(len(valid_smiles), dtype=np.float32),
        }
    )
    dataset_dir.mkdir(parents=True, exist_ok=True)
    out_csv = dataset_dir / f"{DEFAULT_DATASET_NAME}.csv"
    out_df.to_csv(out_csv, index=False)
    print(f"已写入下游 CSV: {out_csv} （{len(out_df)} 行）", flush=True)
    return out_csv


def run_preprocess(cache_root: Path, n_jobs: int) -> None:
    scripts_dir = REPO_ROOT / "scripts"
    _run(
        [
            sys.executable,
            "preprocess_downstream_dataset.py",
            "--data_path",
            str(cache_root),
            "--dataset",
            DEFAULT_DATASET_NAME,
            "--n_jobs",
            str(n_jobs),
        ],
        cwd=scripts_dir,
    )


def run_extract_features(cache_root: Path, model_path: Path, config: str) -> None:
    scripts_dir = REPO_ROOT / "scripts"
    _run(
        [
            sys.executable,
            "extract_features.py",
            "--config",
            config,
            "--model_path",
            str(model_path),
            "--data_path",
            str(cache_root),
            "--dataset",
            DEFAULT_DATASET_NAME,
        ],
        cwd=scripts_dir,
    )


def merge_to_pickle(npz_path: Path, dataset_csv: Path, output_pickle: Path) -> None:
    pack = np.load(npz_path)
    fps = pack["fps"]
    df = pd.read_csv(dataset_csv)
    if len(df) != len(fps):
        raise RuntimeError(f"行数不一致: csv={len(df)} vs fps={len(fps)}")

    smiles_list = df["smiles"].astype(str).tolist()
    mapping = {s: np.asarray(fps[i], dtype=np.float32) for i, s in enumerate(smiles_list)}
    output_pickle.parent.mkdir(parents=True, exist_ok=True)
    with open(output_pickle, "wb") as f:
        pickle.dump(mapping, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"已保存: {output_pickle} （{len(mapping)} 条, 维度 {fps.shape[1]}）", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="将 smiles CSV 转为同目录 KPGT.pickle")
    parser.add_argument("--smiles_csv", type=Path, required=True, help="输入 CSV 路径（需含 canonical_smiles 列）")
    parser.add_argument("--model_path", type=Path, default=DEFAULT_MODEL_PATH, help="预训练权重 base.pth 路径")
    parser.add_argument("--config", type=str, default="base", help="model_config 配置名")
    parser.add_argument("--n_jobs", type=int, default=16, help="预处理并行进程数")
    parser.add_argument("--max-molecules", type=int, default=None, help="仅处理前 N 条（调试）")
    parser.add_argument("--cache_root", type=Path, default=DEFAULT_CACHE_ROOT, help="中间缓存根目录")
    parser.add_argument("--keep_cache", action="store_true", help="保留 feature_extract/.cache 中间文件")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    smiles_csv = args.smiles_csv.resolve()
    if not smiles_csv.is_file():
        print(f"未找到输入 CSV: {smiles_csv}", file=sys.stderr)
        sys.exit(1)

    model_path = args.model_path.resolve()
    if not model_path.is_file():
        print(
            f"未找到预训练权重: {model_path}\n"
            "请从 https://figshare.com/s/d488f30c23946cf6898f 下载并解压到 models/pretrained/base/",
            file=sys.stderr,
        )
        sys.exit(1)

    output_pickle = smiles_csv.parent / "KPGT.pickle"
    cache_root = args.cache_root.resolve()
    dataset_dir = cache_root / DEFAULT_DATASET_NAME
    dataset_csv = dataset_dir / f"{DEFAULT_DATASET_NAME}.csv"
    npz_path = dataset_dir / f"kpgt_{args.config}.npz"

    # 每次运行前清空旧缓存，避免不同任务间残留文件互相污染。
    if dataset_dir.exists():
        shutil.rmtree(dataset_dir)

    try:
        prepare_downstream_csv(smiles_csv, dataset_dir, args.max_molecules)
        run_preprocess(cache_root, args.n_jobs)
        run_extract_features(cache_root, model_path, args.config)
        merge_to_pickle(npz_path, dataset_csv, output_pickle)
    finally:
        if not args.keep_cache and dataset_dir.exists() and output_pickle.exists():
            shutil.rmtree(dataset_dir)
            if cache_root.exists() and not any(cache_root.iterdir()):
                cache_root.rmdir()
            print(f"已清理缓存目录: {dataset_dir}", flush=True)


if __name__ == "__main__":
    main()
