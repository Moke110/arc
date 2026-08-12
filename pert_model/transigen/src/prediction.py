from __future__ import annotations

import argparse
from pathlib import Path
import pickle
import sys
import warnings

import numpy as np
import pandas as pd
import torch
from cmapPy.pandasGEXpress.parse import parse

import src.dataset as src_dataset
import src.model as src_model
import src.utils as src_utils
from src.dataset import TranSiGenDatasetScreening
from src.utils import save_to_HDF, seed_worker, setup_seed

warnings.filterwarnings("ignore")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Arguments for TranSiGen phenotype prediction")
    parser.add_argument("--model_path", type=str, default="training/models/best_model.pt")
    parser.add_argument("--predict_data_root", type=str, default="prediction/example")
    parser.add_argument("--compound_info", type=str, default="screening_compound.csv")
    parser.add_argument("--compound_emb", type=str, default="KPGT_emb.pickle")
    parser.add_argument("--cell", type=str, default="YAPC")
    parser.add_argument("--seed", type=int, default=364039)
    parser.add_argument("--dev", type=str, default="cuda:0")
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--num_workers", type=int, default=2)
    return parser.parse_args()


def _resolve_device(dev_arg: str) -> torch.device:
    if not torch.cuda.is_available():
        return torch.device("cpu")
    requested = dev_arg.strip()
    if requested == "cuda":
        return torch.device("cuda:0")
    if not requested.startswith("cuda"):
        return torch.device("cpu")
    try:
        idx = int(requested.split(":", 1)[1])
    except (IndexError, ValueError):
        return torch.device("cuda:0")
    if idx < 0 or idx >= torch.cuda.device_count():
        return torch.device("cuda:0")
    return torch.device(requested)


def _load_legacy_model(model_path: Path, device: torch.device):
    # Torch 2.6+ defaults to weights_only=True; legacy full-object checkpoints need class imports.
    sys.modules.setdefault("dataset", src_dataset)
    sys.modules.setdefault("model", src_model)
    sys.modules.setdefault("utils", src_utils)
    model = torch.load(model_path, map_location="cpu", weights_only=False)
    model.dev = device
    model.to(device)
    return model


def _load_infer_weight(geneinfo_path: Path, infer_weight_path: Path) -> pd.DataFrame:
    df_gene = pd.read_csv(geneinfo_path)
    df_landmark_gene = df_gene[(df_gene["pr_is_bing"] == 1) & (df_gene["pr_is_lm"] == 1)]
    df_best_infer_gene = df_gene[(df_gene["pr_is_bing"] == 1) & (df_gene["pr_is_lm"] == 0)]
    landmark_ids = df_landmark_gene["pr_id"].tolist()
    best_infer_ids = df_best_infer_gene["pr_id"].tolist()

    infer_weight = parse(str(infer_weight_path), cid=["OFFSET"] + landmark_ids, rid=best_infer_ids)
    infer_weight_df_tmp = infer_weight.data_df
    infer_weight_df = infer_weight_df_tmp[["OFFSET"] + landmark_ids]
    infer_weight_df = infer_weight_df.loc[best_infer_ids]
    return infer_weight_df


def prediction_profiles(args: argparse.Namespace) -> None:
    project_root = Path(__file__).resolve().parent.parent
    model_path = project_root / args.model_path
    predict_data_root = project_root / args.predict_data_root
    data_path = predict_data_root / args.compound_info
    molecule_feature_embed_path = predict_data_root / args.compound_emb
    prediction_info_root = project_root / "prediction" / "info"
    geneinfo_path = prediction_info_root / "geneinfo_processed.csv"
    infer_weight_path = prediction_info_root / "infer_weight.gctx"
    modz_x1_path = prediction_info_root / "modz_x1.pickle"
    output_path = predict_data_root / f"prediction_profile_{args.cell}.h5"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    for required_path in [model_path, data_path, molecule_feature_embed_path, geneinfo_path, infer_weight_path, modz_x1_path]:
        if not required_path.is_file():
            raise FileNotFoundError(f"Required file not found: {required_path}")

    infer_weight_df = _load_infer_weight(geneinfo_path, infer_weight_path)
    with modz_x1_path.open("rb") as f:
        dict_modz_x1_all_cid = pickle.load(f)

    dev = _resolve_device(args.dev)
    model = _load_legacy_model(model_path, dev)

    selected_cid = args.cell
    random_seed = args.seed
    if selected_cid not in dict_modz_x1_all_cid:
        raise ValueError(f"Cell id {selected_cid!r} not found in modz_x1 dictionary.")

    with molecule_feature_embed_path.open("rb") as f:
        smi2emb = pickle.load(f)
    df_screening = pd.read_csv(data_path)

    emb_array = []
    mid_array = []
    for _, row in df_screening.iterrows():
        smi = row["canonical_smiles"]
        if smi not in smi2emb:
            continue
        emb_array.append(smi2emb[smi])
        mid_array.append(row["cp_id"])
    if len(emb_array) == 0:
        raise ValueError("No screening compounds matched molecule embedding dictionary.")

    emb_array = np.asarray(emb_array, dtype=np.float32)
    mid_array = np.asarray(mid_array)
    cid_array = np.asarray([selected_cid] * emb_array.shape[0])
    x1_array = dict_modz_x1_all_cid[selected_cid]
    x1_array = np.repeat(x1_array, emb_array.shape[0], axis=0).astype(np.float32)

    test = TranSiGenDatasetScreening(x1=x1_array, mol_feature=emb_array, mid=mid_array, cid=cid_array)
    test_loader = torch.utils.data.DataLoader(
        dataset=test,
        batch_size=args.batch_size,
        shuffle=False,
        drop_last=False,
        num_workers=args.num_workers,
        worker_init_fn=seed_worker,
    )

    setup_seed(random_seed)
    x2_pred_array, mid_array, cid_array = model.predict_profile_for_x1(test_loader)
    ddict_data: dict[str, np.ndarray] = {
        "x1": x1_array,
        "x2_pred": x2_pred_array,
        "mid": mid_array,
        "cid": cid_array,
    }

    sig_ls = [f"{int(ddict_data['mid'][idx])}_{ddict_data['cid'][idx]}" for idx in range(ddict_data["cid"].shape[0])]
    ddict_data["sig"] = np.asarray(sig_ls)

    for data_type in ["x1", "x2_pred"]:
        x = ddict_data[data_type]
        inferred = np.dot(x, np.asarray(infer_weight_df.T[1:])) + np.asarray(infer_weight_df.T.loc["OFFSET"])
        x_tmp = np.concatenate((x, inferred), axis=1)
        ddict_data[f"{data_type}_inferred"] = x_tmp

    for k in ddict_data:
        print(k, ddict_data[k].shape)
    save_to_HDF(str(output_path), ddict_data)
    print(f"Saved prediction profiles to: {output_path}")


if __name__ == "__main__":
    prediction_profiles(parse_args())
