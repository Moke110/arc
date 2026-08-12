from __future__ import annotations

import pickle
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.dataset import TranSiGenDataset
from src.model import TranSiGen
from src.utils import save_to_HDF, setup_seed


def build_test_input_files(test_dir: Path) -> dict[str, Path]:
    setup_seed(42)
    n_samples = 8
    n_genes = 978
    n_features = 2304
    n_molecules = 4

    mids = np.arange(n_molecules)
    mid2smi = {int(i): f"SMI_{i}" for i in mids}
    smi2emb = {f"SMI_{i}": np.random.randn(n_features).astype(np.float32) for i in mids}

    sample_mol = np.array([0, 1, 2, 3, 0, 1, 2, 3], dtype=np.int64)
    cids = np.array(["CELL_A", "CELL_A", "CELL_B", "CELL_B", "CELL_A", "CELL_A", "CELL_B", "CELL_B"])
    sig = np.array([f"{cp}_{cid}" for cp, cid in zip(sample_mol, cids)])

    x1 = np.random.randn(n_samples, n_genes).astype(np.float32)
    perturb = np.random.randn(n_samples, n_genes).astype(np.float32) * 0.05
    x2 = x1 + perturb

    input_dir = test_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)

    input_h5 = input_dir / "processed_data.h5"
    save_to_HDF(
        str(input_h5),
        {
            "x1": x1,
            "x2": x2,
            "sig": sig,
        },
    )

    mid2smi_file = input_dir / "mid2smi.pickle"
    with mid2smi_file.open("wb") as f:
        pickle.dump(mid2smi, f)
    kpgt_file = input_dir / "example_KPGT.pickle"
    with kpgt_file.open("wb") as f:
        pickle.dump(smi2emb, f)

    split_file = input_dir / "dataset_index.csv"
    pd.DataFrame(
        {
            "LINCS_index": np.arange(n_samples),
            "mid": sample_mol,
            "cid": cids,
        }
    ).to_csv(split_file, index=False)

    return {
        "input_h5": input_h5,
        "mid2smi": mid2smi_file,
        "kpgt": kpgt_file,
        "split": split_file,
    }


def run_end_to_end_test(test_dir: Path) -> None:
    input_dir = test_dir / "input"
    output_dir = test_dir / "output"
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    inputs = build_test_input_files(test_dir)
    index_df = pd.read_csv(inputs["split"])
    dataset = TranSiGenDataset(
        LINCS_index=index_df["LINCS_index"].to_numpy(),
        mol_feature_type="KPGT",
        mid=index_df["mid"].to_numpy(),
        cid=index_df["cid"].to_numpy(),
        data_root=input_dir,
        processed_data_file=inputs["input_h5"].name,
        mid2smi_file=inputs["mid2smi"].name,
        molecule_emb_file=inputs["kpgt"].name,
        molecule_emb_root=input_dir,
    )
    dataloader = torch.utils.data.DataLoader(dataset=dataset, batch_size=4, shuffle=False, drop_last=False, num_workers=0)

    model = TranSiGen(
        n_genes=978,
        n_latent=64,
        n_en_hidden=[256, 128],
        n_de_hidden=[128, 256],
        features_dim=2304,
        features_embed_dim=[256],
        beta=0.1,
        dropout=0.1,
        device=torch.device("cpu"),
    )
    model.eval()

    x2_pred_all = []
    mid_all = []
    cid_all = []
    with torch.no_grad():
        for x1_batch, _, feat_batch, mid_batch, cid_batch, _ in dataloader:
            _, _, _, x2_pred, _, _, _ = model(x1_batch, feat_batch)
            x2_pred_all.append(x2_pred.numpy())
            mid_all.extend(mid_batch.numpy().tolist())
            cid_all.extend(list(cid_batch))

    x2_pred_all = np.concatenate(x2_pred_all, axis=0)
    output_h5 = output_dir / "predicted_profiles.h5"
    save_to_HDF(
        str(output_h5),
        {
            "x2_pred": x2_pred_all.astype(np.float32),
            "mid": np.asarray(mid_all, dtype=np.int64),
            "cid": np.asarray(cid_all),
        },
    )

    output_csv = output_dir / "prediction_summary.csv"
    pd.DataFrame(
        {
            "mid": mid_all,
            "cid": cid_all,
            "pred_mean": x2_pred_all.mean(axis=1),
            "pred_std": x2_pred_all.std(axis=1),
        }
    ).to_csv(output_csv, index=False)

    (test_dir / "README.txt").write_text(
        "# 1_update_model 测试目录说明\n\n"
        "本目录用于验证升级后的 TranSiGen 核心模块（dataset/model/vae）是否可以完整执行一次\n"
        "“读取数据 -> 模型预测 -> 输出结果”的流程。\n\n"
        "## 目录结构\n"
        "- input/: 测试输入数据\n"
        "- output/: 模型输出结果\n"
        "- run_test.py: 端到端测试脚本\n\n"
        "## 输入文件（input）\n"
        "- processed_data.h5\n"
        "  - 作用: 测试表达数据主文件，包含 x1、x2、sig。\n"
        "  - 说明: x1 为基线表达，x2 为扰动后表达，sig 为样本唯一标识。\n"
        "- dataset_index.csv\n"
        "  - 作用: 数据索引与映射文件，驱动 Dataset 逐样本读取。\n"
        "  - 字段: LINCS_index, mid, cid。\n"
        "- mid2smi.pickle\n"
        "  - 作用: 分子整数 ID（mid）到 SMILES 字符串的映射字典。\n"
        "- example_KPGT.pickle\n"
        "  - 作用: SMILES到分子向量（2304维）的映射字典。\n\n"
        "## 输出文件（output）\n"
        "- predicted_profiles.h5\n"
        "  - 作用: 保存模型预测的完整结果。\n"
        "  - 内容: x2_pred（预测表达谱）、mid、cid。\n"
        "- prediction_summary.csv\n"
        "  - 作用: 预测结果的样本级摘要，便于快速检查。\n"
        "  - 字段: mid, cid, pred_mean, pred_std。\n\n"
        "## 运行方式\n"
        "- 在项目根目录执行:\n"
        "  - .venv/bin/python test/1_update_model/run_test.py\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    base_dir = Path(__file__).resolve().parent
    run_end_to_end_test(base_dir)
