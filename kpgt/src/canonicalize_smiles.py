#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem


def to_canonical_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def pick_smiles_column(df: pd.DataFrame) -> str:
    for name in ("smiles", "canonical_smiles"):
        if name in df.columns:
            return name
    raise ValueError(f"CSV 缺少 SMILES 列，当前列: {list(df.columns)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="将 CSV 中的 SMILES 转为 canonical SMILES")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("feature_extract/INO/smiles.csv"),
        help="输入 CSV 路径",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("feature_extract/INO/smiles_canonical.csv"),
        help="输出 CSV 路径",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = args.input.resolve()
    output_path = args.output.resolve()

    if not input_path.is_file():
        raise FileNotFoundError(f"未找到输入文件: {input_path}")

    df = pd.read_csv(input_path)
    smiles_col = pick_smiles_column(df)

    src = df[smiles_col].fillna("").astype(str).str.strip()
    canonical = src.map(lambda s: to_canonical_smiles(s) if s else "")

    out_df = df.copy()
    out_df["canonical_smiles"] = canonical
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_path, index=False)

    total = len(out_df)
    invalid = int(((src != "") & (canonical == "")).sum())
    print(f"已写入: {output_path}")
    print(f"总条数: {total}, 无法解析: {invalid}")


if __name__ == "__main__":
    main()
