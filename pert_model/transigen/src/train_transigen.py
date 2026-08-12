from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import pickle
import sys
import warnings

import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt

from src.dataset import TranSiGenDataset
from src.model import TranSiGen
from src.utils import load_from_HDF, setup_seed, seed_worker

warnings.filterwarnings("ignore")


def _str2bool(v):
    if isinstance(v, bool):
        return v
    value = str(v).strip().lower()
    if value in {"1", "true", "t", "yes", "y"}:
        return True
    if value in {"0", "false", "f", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {v}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Arguments for training TranSiGen")
    parser.add_argument("--training_data_root", type=str, default="training/data")
    parser.add_argument("--dataset_name", type=str, default="example", help="dataset name without suffix '_set'")
    parser.add_argument("--molecule_emb_root", type=str, default="mol_emb")
    parser.add_argument("--pretrained_vae1_path", type=str, default="training/pretrained/vae1.pt")
    parser.add_argument("--pretrained_vae2_path", type=str, default="training/pretrained/vae2.pt")
    parser.add_argument("--use_pretrained_vae", type=_str2bool, default=True)
    parser.add_argument("--output_root", type=str, default="training/models")

    parser.add_argument("--dev", type=str, default="cuda:2")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--molecule_feature",
        type=str,
        default="KPGT",
        help="Molecule feature type or optional direct .pickle path; defaults to auto-load from dataset directory.",
    )
    parser.add_argument("--n_epochs", type=int, default=300)
    parser.add_argument("--n_latent", type=int, default=100)
    parser.add_argument("--molecule_feature_embed_dim", nargs="+", type=int, default=[400])
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--learning_rate", type=float, default=1e-3)
    parser.add_argument("--beta", type=float, default=0.1)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--weight_decay", type=float, default=1e-5)
    parser.add_argument(
        "--patient_cons",
        type=int,
        default=10,
        help="Early stop when validation loss increases for this many consecutive epochs",
    )
    parser.add_argument(
        "--patience_turn",
        type=int,
        default=30,
        help="Early stop when validation loss has increased this many times in total",
    )

    parser.add_argument("--train_flag", type=_str2bool, default=True)
    parser.add_argument("--eval_metric", type=_str2bool, default=False)
    parser.add_argument("--predict_profile", type=_str2bool, default=False)
    return parser.parse_args()


def _load_pretrained_state_dict(model_path: Path) -> dict:
    # Torch 2.6 defaults to weights_only=True; these legacy files need class-based load.
    import src.vae_x1 as src_vae_x1
    import src.vae_x2 as src_vae_x2

    sys.modules.setdefault("vae_x1", src_vae_x1)
    sys.modules.setdefault("vae_x2", src_vae_x2)
    model_obj = torch.load(model_path, map_location="cpu", weights_only=False)
    return model_obj.state_dict() if hasattr(model_obj, "state_dict") else model_obj


def _resolve_device(dev_arg: str) -> torch.device:
    if not torch.cuda.is_available():
        return torch.device("cpu")

    requested = dev_arg.strip()
    if not requested.startswith("cuda"):
        return torch.device("cpu")

    cuda_count = torch.cuda.device_count()
    if requested == "cuda":
        return torch.device("cuda:0")

    try:
        requested_index = int(requested.split(":", 1)[1])
    except (IndexError, ValueError):
        return torch.device("cuda:0")

    if requested_index < 0 or requested_index >= cuda_count:
        return torch.device("cuda:0")

    return torch.device(requested)


def _rank_desc(values: np.ndarray) -> np.ndarray:
    order = np.argsort(-values)
    ranks = np.empty_like(order, dtype=np.int64)
    ranks[order] = np.arange(1, len(values) + 1)
    return ranks


def _plot_pred_pearson_rmse(x2_true: np.ndarray, x2_pred: np.ndarray, save_path: Path) -> None:
    pearson = float(np.corrcoef(x2_true, x2_pred)[0, 1])
    rmse = float(np.sqrt(np.mean((x2_pred - x2_true) ** 2)))
    fit_k, fit_b = np.polyfit(x2_true, x2_pred, deg=1)

    plt.figure(figsize=(7, 6), dpi=180)
    plt.scatter(x2_true, x2_pred, c="black", s=14, alpha=0.8)
    line_x = np.linspace(float(np.min(x2_true)), float(np.max(x2_true)), 200)
    plt.plot(line_x, fit_k * line_x + fit_b, color="tab:orange", linewidth=2, label="Linear fit")
    plt.plot(line_x, line_x, color="tab:green", linestyle="--", linewidth=2, label="$X'_2 = X_2$")
    plt.xlabel("expression in $X_2$")
    plt.ylabel("expression in $X_2'$")
    plt.title(f"pred_pearson_rmse (Pearson={pearson:.4f}, RMSE={rmse:.4f})")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()


def _plot_precision_at_100(x2_true: np.ndarray, x2_pred: np.ndarray, save_path: Path) -> None:
    true_rank = _rank_desc(x2_true)
    pred_rank = _rank_desc(x2_pred)
    n_genes = len(true_rank)

    point_colors = np.full(n_genes, "gray", dtype=object)
    top_true = true_rank <= 100
    top_pred = pred_rank <= 100
    bottom_true = true_rank > 878
    bottom_pred = pred_rank > 878
    point_colors[top_true & top_pred] = "red"
    point_colors[top_true & ~top_pred] = "lightcoral"
    point_colors[bottom_true & bottom_pred] = "blue"
    point_colors[bottom_true & ~bottom_pred] = "lightskyblue"

    plt.figure(figsize=(7, 6), dpi=180)
    plt.scatter(true_rank, pred_rank, c=point_colors.tolist(), s=18, alpha=0.85)
    plt.axhline(y=100, color="red", linestyle="--", linewidth=2)
    plt.axhline(y=878, color="blue", linestyle="--", linewidth=2)
    plt.xlabel("rank in $X_2$")
    plt.ylabel("rank in $X_2'$")
    plt.title("precision@100 rank scatter")
    plt.xlim(1, n_genes)
    plt.ylim(1, n_genes)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()


def _plot_precision_deg_at_100(deg_true: np.ndarray, deg_pred: np.ndarray, save_path: Path) -> None:
    true_rank = _rank_desc(deg_true)
    pred_rank = _rank_desc(deg_pred)
    n_genes = len(true_rank)

    point_colors = np.full(n_genes, "gray", dtype=object)
    top_true = true_rank <= 100
    top_pred = pred_rank <= 100
    bottom_true = true_rank > 878
    bottom_pred = pred_rank > 878
    point_colors[top_true & top_pred] = "red"
    point_colors[top_true & ~top_pred] = "lightcoral"
    point_colors[bottom_true & bottom_pred] = "blue"
    point_colors[bottom_true & ~bottom_pred] = "lightskyblue"

    plt.figure(figsize=(7, 6), dpi=180)
    plt.scatter(true_rank, pred_rank, c=point_colors.tolist(), s=18, alpha=0.85)
    plt.axhline(y=100, color="red", linestyle="--", linewidth=2)
    plt.axhline(y=878, color="blue", linestyle="--", linewidth=2)
    plt.xlabel("rank in DEG ($X_2-X_1$)")
    plt.ylabel("rank in DEG ($X_2‘-X_1$)")
    plt.title("precision_DEG@100 rank scatter")
    plt.xlim(1, n_genes)
    plt.ylim(1, n_genes)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()


def train_transigen(args: argparse.Namespace) -> None:
    print(args)
    torch.set_num_threads(1)

    random_seed = args.seed
    if random_seed is None:
        # Keep seed within [0, 2**32 - 1] for NumPy compatibility.
        random_seed = int(np.random.SeedSequence().generate_state(1, dtype=np.uint32)[0])
    args.seed = random_seed
    print(f"Using seed: {random_seed}")
    setup_seed(random_seed)
    dev = _resolve_device(args.dev)
    timestamp = datetime.now().strftime("%y%m%d-%H%M")
    _ds = str(args.dataset_name).strip().replace("/", "_").replace("\\", "_") or "dataset"
    run_name = f"{_ds}_{timestamp}"
    output_dir = Path(args.output_root) / run_name
    output_dir.mkdir(parents=True, exist_ok=True)

    dataset_dir = Path(args.training_data_root) / f"{args.dataset_name}_set"
    split_dir = dataset_dir / "split"
    processed_data_path = dataset_dir / "processed_data.h5"
    mid2smi_path = dataset_dir / "mid2smi.pickle"
    with open(mid2smi_path, "rb") as f:
        mid2smi = pickle.load(f)

    feat_arg = str(args.molecule_feature).strip()
    emb_candidate = Path(feat_arg).expanduser()
    if emb_candidate.is_file() and emb_candidate.suffix.lower() == ".pickle":
        emb_file = emb_candidate.resolve()
        feat_type = emb_file.stem
    elif feat_arg in {"KPGT", "ECFP4"}:
        emb_file = dataset_dir / f"{feat_arg}_emb.pickle"
        if not emb_file.is_file():
            raise FileNotFoundError(
                f"Auto embedding file not found: {emb_file}. "
                f"Please place `{feat_arg}_emb.pickle` under `{dataset_dir}`."
            )
        feat_type = feat_arg
    else:
        raise ValueError(f"--molecule_feature must be .pickle path or KPGT / ECFP4, got {feat_arg!r}")

    with emb_file.open("rb") as f:
        _smi2emb_probe = pickle.load(f)
    _sample = next(iter(_smi2emb_probe.values()))
    features_dim = int(np.asarray(_sample).reshape(-1).shape[0])
    molecule_emb_root = emb_file.parent

    train_pair = load_from_HDF(str(split_dir / "train_id.h5"))
    valid_pair = load_from_HDF(str(split_dir / "valid_id.h5"))
    test_pair = load_from_HDF(str(split_dir / "test_id.h5"))
    processed_data = load_from_HDF(str(processed_data_path))
    max_available_index = len(processed_data["x1"]) - 1
    max_required_index = int(
        max(
            train_pair["LINCS_index"].max(),
            valid_pair["LINCS_index"].max(),
            test_pair["LINCS_index"].max(),
        )
    )
    if max_required_index > max_available_index:
        raise ValueError(
            "Split files reference LINCS indices outside processed_data range. "
            f"max_required_index={max_required_index}, max_available_index={max_available_index}. "
            "Please provide a full processed_data.h5 or use split files matched to the current processed_data."
        )

    print("=============== dataset split ================")
    print(
        "train",
        len(set(train_pair["cid"])),
        len(train_pair["mid"]),
        len(set(train_pair["mid"])),
    )
    print(
        "valid",
        len(set(valid_pair["cid"])),
        len(valid_pair["mid"]),
        len(set(valid_pair["mid"])),
    )
    print(
        "test",
        len(set(test_pair["cid"])),
        len(test_pair["mid"]),
        len(set(test_pair["mid"])),
    )

    data_root = dataset_dir
    processed_data_file = "processed_data.h5"
    train_set = TranSiGenDataset(
        LINCS_index=train_pair["LINCS_index"],
        mol_feature_type=feat_type,
        mid=train_pair["mid"],
        cid=train_pair["cid"],
        data_root=data_root,
        processed_data_file=processed_data_file,
        molecule_emb_file=emb_file.name,
        molecule_emb_root=molecule_emb_root,
    )
    valid_set = TranSiGenDataset(
        LINCS_index=valid_pair["LINCS_index"],
        mol_feature_type=feat_type,
        mid=valid_pair["mid"],
        cid=valid_pair["cid"],
        data_root=data_root,
        processed_data_file=processed_data_file,
        molecule_emb_file=emb_file.name,
        molecule_emb_root=molecule_emb_root,
    )
    test_set = TranSiGenDataset(
        LINCS_index=test_pair["LINCS_index"],
        mol_feature_type=feat_type,
        mid=test_pair["mid"],
        cid=test_pair["cid"],
        data_root=data_root,
        processed_data_file=processed_data_file,
        molecule_emb_file=emb_file.name,
        molecule_emb_root=molecule_emb_root,
    )

    train_loader = torch.utils.data.DataLoader(
        dataset=train_set,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=False,
        num_workers=2,
        worker_init_fn=seed_worker,
    )
    valid_loader = torch.utils.data.DataLoader(
        dataset=valid_set,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=False,
        num_workers=2,
        worker_init_fn=seed_worker,
    )
    test_loader = torch.utils.data.DataLoader(
        dataset=test_set,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=False,
        num_workers=2,
        worker_init_fn=seed_worker,
    )

    best_epoch = None
    model = TranSiGen(
        n_genes=978,
        n_latent=args.n_latent,
        n_en_hidden=[1200],
        n_de_hidden=[800],
        features_dim=features_dim,
        features_embed_dim=args.molecule_feature_embed_dim,
        init_w=True,
        beta=args.beta,
        device=dev,
        dropout=args.dropout,
        path_model=f"{output_dir.as_posix()}/",
        random_seed=random_seed,
    )
    model_dict = model.state_dict()

    if args.use_pretrained_vae:
        print("===== load vae for x1 and x2 =====")
        vae1_dict = _load_pretrained_state_dict(Path(args.pretrained_vae1_path))
        vae2_dict = _load_pretrained_state_dict(Path(args.pretrained_vae2_path))
        for key in model_dict:
            if key in vae1_dict:
                model_dict[key] = vae1_dict[key]
            if key in vae2_dict:
                model_dict[key] = vae2_dict[key]
        model.load_state_dict(model_dict)
    model.to(dev)

    if args.train_flag:
        epoch_hist, best_epoch = model.train_model(
            train_loader=train_loader,
            test_loader=valid_loader,
            n_epochs=args.n_epochs,
            learning_rate=args.learning_rate,
            weight_decay=args.weight_decay,
            save_model=True,
            patient_cons=args.patient_cons,
            patience_turn=args.patience_turn,
        )
        epoch_result = pd.DataFrame.from_dict(epoch_hist)
        epoch_result["epoch"] = list(range(len(epoch_result)))
        epoch_result.to_csv(output_dir / f"epoch{args.n_epochs}_lr{args.learning_rate}.csv", index=False)

    best_main_src = output_dir / "best_model.pt"
    if not best_main_src.exists():
        # Fallback: keep current model when validation loss becomes NaN.
        torch.save(model, best_main_src)
    loss_items = ["loss", "mse_x1", "mse_x2", "mse_pert", "kld_x1", "kld_x2", "kld_pert"]
    metrics_items = ["pearson", "rmse", "precision100"]
    model_loaded = torch.load(best_main_src, map_location="cpu", weights_only=False)
    model_loaded.dev = torch.device(dev)
    model_loaded.to(dev)
    test_loss_dict, test_metrics_dict, test_metrics_dict_ls = model_loaded.test_model(
        loader=test_loader,
        loss_item=loss_items,
        metrics_func=metrics_items,
    )
    x1_test, x2_test, _, _, x2_pred_test, _, _, _, _ = model_loaded.predict_profile(test_loader)
    # Plot with one point per landmark gene (978 points) using test-set mean expression.
    x2_test_mean = np.mean(x2_test, axis=0)
    x2_pred_test_mean = np.mean(x2_pred_test, axis=0)
    deg_test_mean = np.mean(x2_test - x1_test, axis=0)
    deg_pred_test_mean = np.mean(x2_pred_test - x1_test, axis=0)
    pred_pearson_rmse_path = output_dir / "pred_pearson_rmse.png"
    precision_pred100_path = output_dir / "precision_pred@100.png"
    precision_deg100_path = output_dir / "precision_DEG@100.png"
    _plot_pred_pearson_rmse(x2_test_mean, x2_pred_test_mean, pred_pearson_rmse_path)
    _plot_precision_at_100(x2_test_mean, x2_pred_test_mean, precision_pred100_path)
    _plot_precision_deg_at_100(deg_test_mean, deg_pred_test_mean, precision_deg100_path)

    def _fmt_metrics(loss_dict: dict, metrics_dict: dict) -> str:
        merged = {**loss_dict, **metrics_dict}
        ordered_keys = [
            "loss",
            "mse_x1",
            "mse_x2",
            "mse_pert",
            "kld_x1",
            "kld_x2",
            "kld_pert",
            "x1_rec_pearson",
            "x2_rec_pearson",
            "x2_pred_pearson",
            "DEG_rec_pearson",
            "DEG_pred_pearson",
            "x2_pred_rmse",
            "DEG_pred_rmse",
            "x2_pred_negative_precision100",
            "x2_pred_positive_precision100",
            "DEG_pred_negative_precision100",
            "DEG_pred_positive_precision100",
        ]
        lines = []
        for key in ordered_keys:
            if key in merged:
                lines.append(f"- {key}: `{merged[key]:.6f}`")
        return "\n".join(lines)

    summary_path = output_dir / f"TrainLog_{run_name}.md"
    best_epoch_text = "N/A" if best_epoch is None else str(best_epoch + 1)
    summary = (
        f"# Training Summary ({run_name})\n\n"
        f"- Device: `{dev}`\n"
        f"- Seed: `{random_seed}`\n"
        f"- Molecule feature: `{feat_type}`\n"
        f"- Use pretrained VAE: `{args.use_pretrained_vae}`\n"
        f"- Dataset: `{args.dataset_name}_set`\n"
        f"- Epochs: `{args.n_epochs}`\n"
        f"- Batch size: `{args.batch_size}`\n"
        f"- Learning rate: `{args.learning_rate}`\n"
        f"- Early stop patient_cons: `{args.patient_cons}`\n"
        f"- Early stop patience_turn: `{args.patience_turn}`\n"
        f"- Best epoch: `{best_epoch_text}`\n\n"
        f"## Dataset Stats\n"
        f"- Train samples: `{len(train_pair['mid'])}`\n"
        f"- Valid samples: `{len(valid_pair['mid'])}`\n"
        f"- Test samples: `{len(test_pair['mid'])}`\n\n"
        f"## Test set performance\n"
        f"{_fmt_metrics(test_loss_dict, test_metrics_dict)}\n\n"
        f"## Output Files\n"
        f"- `{best_main_src}`\n"
        f"- `{pred_pearson_rmse_path}`\n"
        f"- `{precision_pred100_path}`\n"
        f"- `{precision_deg100_path}`\n"
    )
    summary_path.write_text(summary, encoding="utf-8")

    print("Training artifacts saved:")
    print(f"- Best model: {best_main_src}")
    print(f"- Summary: {summary_path}")

    if args.eval_metric:
        df_rec = pd.DataFrame.from_dict(test_metrics_dict_ls)
        df_rec["canonical_smiles"] = [mid2smi[int(m)] for m in df_rec["mid"]]
        metric_path = output_dir / "test_reconstruction_result_all_samples.csv"
        df_rec.to_csv(metric_path, index=False)
        print(f"- Evaluation metrics: {metric_path}")


if __name__ == "__main__":
    train_transigen(parse_args())
