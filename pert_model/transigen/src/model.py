from __future__ import annotations

from collections import defaultdict
import copy

import numpy as np
import torch
import torch.nn.functional as F
from torch import nn, optim

from .utils import get_metric_func


class TranSiGen(torch.nn.Module):
    def __init__(self, n_genes, n_latent, n_en_hidden, n_de_hidden, features_dim, features_embed_dim, **kwargs):
        super().__init__()
        self.n_genes = n_genes
        self.n_latent = n_latent
        self.n_en_hidden = n_en_hidden
        self.n_de_hidden = n_de_hidden
        self.features_dim = features_dim
        self.feat_embed_dim = features_embed_dim
        self.init_w = kwargs.get("init_w", False)
        self.beta = kwargs.get("beta", 0.05)
        self.path_model = kwargs.get("path_model", "trained_model/")
        self.dev = kwargs.get("device", torch.device("cpu"))
        self.dropout = kwargs.get("dropout", 0.3)

        encoder = [nn.Linear(self.n_genes, self.n_en_hidden[0]), nn.BatchNorm1d(self.n_en_hidden[0]), nn.ReLU(), nn.Dropout(self.dropout)]
        if len(n_en_hidden) > 1:
            for i in range(len(n_en_hidden) - 1):
                encoder.extend([nn.Linear(self.n_en_hidden[i], self.n_en_hidden[i + 1]), nn.BatchNorm1d(self.n_en_hidden[i + 1]), nn.ReLU(), nn.Dropout(self.dropout)])
        self.encoder_x2 = nn.Sequential(*encoder)
        self.encoder_x1 = nn.Sequential(*copy.deepcopy(encoder))
        self.mu_z2 = nn.Linear(self.n_en_hidden[-1], self.n_latent)
        self.logvar_z2 = nn.Linear(self.n_en_hidden[-1], self.n_latent)
        self.mu_z1 = nn.Linear(self.n_en_hidden[-1], self.n_latent)
        self.logvar_z1 = nn.Linear(self.n_en_hidden[-1], self.n_latent)

        if len(n_de_hidden) == 0:
            decoder = [nn.Linear(self.n_latent, self.n_genes)]
        else:
            decoder = [nn.Linear(self.n_latent, self.n_de_hidden[0]), nn.BatchNorm1d(self.n_de_hidden[0]), nn.ReLU(), nn.Dropout(self.dropout)]
            if len(n_de_hidden) > 1:
                for i in range(len(self.n_de_hidden) - 1):
                    decoder.extend([nn.Linear(self.n_de_hidden[i], self.n_de_hidden[i + 1]), nn.BatchNorm1d(self.n_de_hidden[i + 1]), nn.ReLU(), nn.Dropout(self.dropout)])
            decoder.extend([nn.Linear(self.n_de_hidden[-1], self.n_genes), nn.ReLU()])
        self.decoder_x2 = nn.Sequential(*decoder)
        self.decoder_x1 = nn.Sequential(*copy.deepcopy(decoder))

        if self.feat_embed_dim is None:
            feat_in = self.features_dim
        else:
            feat_embeddings = [nn.Linear(self.features_dim, self.feat_embed_dim[0]), nn.BatchNorm1d(self.feat_embed_dim[0]), nn.ReLU(), nn.Dropout(self.dropout)]
            if len(self.feat_embed_dim) > 1:
                for i in range(len(self.feat_embed_dim) - 1):
                    feat_embeddings.extend([nn.Linear(self.feat_embed_dim[i], self.feat_embed_dim[i + 1]), nn.BatchNorm1d(self.feat_embed_dim[i + 1]), nn.ReLU(), nn.Dropout(self.dropout)])
            self.feat_embeddings = nn.Sequential(*feat_embeddings)
            feat_in = self.feat_embed_dim[-1]
        self.mu_z2Fz1 = nn.Linear(self.n_latent + feat_in, self.n_latent)
        self.logvar_z2Fz1 = nn.Linear(self.n_latent + feat_in, self.n_latent)

        if self.init_w:
            for module in [self.encoder_x1, self.decoder_x1, self.encoder_x2, self.decoder_x2, self.mu_z1, self.logvar_z1, self.mu_z2, self.logvar_z2]:
                module.apply(self._init_weights)

    @staticmethod
    def _init_weights(layer):
        if isinstance(layer, nn.Linear):
            torch.nn.init.xavier_uniform_(layer.weight)

    def sample_latent(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std, device=std.device)
        return eps * std + mu

    def encode_x1(self, x):
        y = self.encoder_x1(x)
        mu, logvar = self.mu_z1(y), self.logvar_z1(y)
        return self.sample_latent(mu, logvar), mu, logvar

    def encode_x2(self, x):
        y = self.encoder_x2(x)
        mu, logvar = self.mu_z2(y), self.logvar_z2(y)
        return self.sample_latent(mu, logvar), mu, logvar

    def decode_x1(self, z):
        return self.decoder_x1(z)

    def decode_x2(self, z):
        return self.decoder_x2(z)

    def forward(self, x1, features):
        z1, mu1, logvar1 = self.encode_x1(x1)
        x1_rec = self.decode_x1(z1)
        feat_embed = self.feat_embeddings(features) if self.feat_embed_dim is not None else features
        z1_feat = torch.cat([z1, feat_embed], dim=1)
        mu_pred, logvar_pred = self.mu_z2Fz1(z1_feat), self.logvar_z2Fz1(z1_feat)
        z2_pred = self.sample_latent(mu_pred, logvar_pred)
        x2_pred = self.decode_x2(z2_pred)
        return x1_rec, mu1, logvar1, x2_pred, mu_pred, logvar_pred, z2_pred

    def loss(self, x1_train, x1_rec, mu1, logvar1, x2_train, x2_rec, mu2, logvar2, x2_pred, mu_pred, logvar_pred):
        mse_x1 = F.mse_loss(x1_rec, x1_train, reduction="sum")
        mse_x2 = F.mse_loss(x2_rec, x2_train, reduction="sum")
        mse_pert = F.mse_loss(x2_pred - x1_train, x2_train - x1_train, reduction="sum")
        kld_x1 = -0.5 * torch.sum(1.0 + logvar1 - mu1.pow(2) - logvar1.exp())
        kld_x2 = -0.5 * torch.sum(1.0 + logvar2 - mu2.pow(2) - logvar2.exp())
        kld_pert = -0.5 * torch.sum(1.0 + (logvar_pred - logvar2) - ((mu_pred - mu2).pow(2) + logvar_pred.exp()) / logvar2.exp())
        total = mse_x1 + mse_x2 + mse_pert + self.beta * (kld_x1 + kld_x2 + kld_pert)
        return total, mse_x1, mse_x2, mse_pert, kld_x1, kld_x2, kld_pert

    def train_model(
        self,
        learning_rate,
        weight_decay,
        n_epochs,
        train_loader,
        test_loader,
        save_model=True,
        metrics_func=None,
        patient_cons=10,
        patience_turn=30,
    ):
        epoch_hist = defaultdict(list)
        optimizer = optim.Adam(self.parameters(), lr=learning_rate, weight_decay=weight_decay)
        loss_item = ["loss", "mse_x1", "mse_x2", "mse_pert", "kld_x1", "kld_x2", "kld_pert"]
        best_value = np.inf
        best_epoch = 0
        valid_loss_increase_cons = 0
        valid_loss_increase_turn = 0
        prev_valid_loss = None
        for epoch in range(n_epochs):
            self.train()
            for x1_train, x2_train, features, _, _, _ in train_loader:
                x1_train = x1_train.to(self.dev)
                x2_train = x2_train.to(self.dev)
                features = features.to(self.dev)
                if x1_train.shape[0] == 1:
                    continue
                optimizer.zero_grad()
                x1_rec, mu1, logvar1, x2_pert, mu_pred, logvar_pred, _ = self.forward(x1_train, features)
                z2, mu2, logvar2 = self.encode_x2(x2_train)
                x2_rec = self.decode_x2(z2)
                loss, *_ = self.loss(x1_train, x1_rec, mu1, logvar1, x2_train, x2_rec, mu2, logvar2, x2_pert, mu_pred, logvar_pred)
                loss.backward()
                optimizer.step()

            train_dict, train_metrics_dict, _ = self.test_model(loader=train_loader, loss_item=loss_item, metrics_func=metrics_func)
            test_dict, test_metrics_dict, _ = self.test_model(loader=test_loader, loss_item=loss_item, metrics_func=metrics_func)
            for k, v in train_dict.items():
                epoch_hist["train_" + k].append(v)
            for k, v in train_metrics_dict.items():
                epoch_hist["train_" + k].append(v)
            for k, v in test_dict.items():
                epoch_hist["valid_" + k].append(v)
            for k, v in test_metrics_dict.items():
                epoch_hist["valid_" + k].append(v)
            current_valid_loss = test_dict["loss"]
            if prev_valid_loss is not None and current_valid_loss > prev_valid_loss:
                valid_loss_increase_cons += 1
                valid_loss_increase_turn += 1
            else:
                valid_loss_increase_cons = 0
            prev_valid_loss = current_valid_loss
            if test_dict["loss"] < best_value:
                best_value = test_dict["loss"]
                best_epoch = epoch
                if save_model:
                    torch.save(self, self.path_model + "best_model.pt")
            print(
                f"[Epoch {epoch + 1}/{n_epochs}] "
                f"train_loss={train_dict.get('loss', float('nan')):.6f} "
                f"valid_loss={test_dict.get('loss', float('nan')):.6f} "
                f"best_epoch={best_epoch + 1} "
                f"best_valid_loss={best_value:.6f}"
            )
            if valid_loss_increase_cons >= patient_cons:
                print(
                    f"Early stopping: valid_loss increased for {valid_loss_increase_cons} "
                    f"consecutive epochs (patient_cons={patient_cons})."
                )
                break
            if valid_loss_increase_turn >= patience_turn:
                print(
                    f"Early stopping: valid_loss increased {valid_loss_increase_turn} "
                    f"times in total (patience_turn={patience_turn})."
                )
                break
        return epoch_hist, best_epoch

    def test_model(self, loader, loss_item=None, metrics_func=None):
        test_dict = defaultdict(float)
        metrics_dict_all = defaultdict(float)
        metrics_dict_all_ls = defaultdict(list)
        test_size = 0
        self.eval()
        with torch.no_grad():
            for x1_data, x2_data, mol_features, mid, cid, sig in loader:
                x1_data = x1_data.to(self.dev)
                x2_data = x2_data.to(self.dev)
                mol_features = mol_features.to(self.dev)
                cid = np.array(list(cid))
                sig = np.array(list(sig))
                test_size += x1_data.shape[0]
                x1_rec, mu1, logvar1, x2_pred, mu_pred, logvar_pred, _ = self.forward(x1_data, mol_features)
                z2, mu2, logvar2 = self.encode_x2(x2_data)
                x2_rec = self.decode_x2(z2)
                loss_ls = self.loss(x1_data, x1_rec, mu1, logvar1, x2_data, x2_rec, mu2, logvar2, x2_pred, mu_pred, logvar_pred)
                if loss_item is not None:
                    for loss_i, k in enumerate(loss_item):
                        test_dict[k] += float(loss_ls[loss_i].item())
                if metrics_func is not None:
                    metrics_dict, metrics_dict_ls = self.eval_x_reconstruction(x1_data, x1_rec, x2_data, x2_rec, x2_pred, metrics_func=metrics_func)
                    for k in metrics_dict:
                        metrics_dict_all[k] += metrics_dict[k]
                        metrics_dict_all_ls[k] += metrics_dict_ls[k]
                    metrics_dict_all_ls["mid"] += list(mid.cpu().numpy())
                    metrics_dict_all_ls["cid"] += list(cid)
                    metrics_dict_all_ls["sig"] += list(sig)
        if test_size > 0:
            for k in test_dict:
                test_dict[k] /= test_size
            for k in metrics_dict_all:
                metrics_dict_all[k] /= test_size
        return test_dict, metrics_dict_all, metrics_dict_all_ls

    def predict_profile(self, loader):
        self.eval()
        with torch.no_grad():
            for x1_data, x2_data, mol_features, mid, cid, sig in loader:
                cid = np.array(list(cid))
                sig = np.array(list(sig))
                x1_data = x1_data.to(self.dev)
                x2_data = x2_data.to(self.dev)
                mol_features = mol_features.to(self.dev)
                x1_rec, _, _, x2_pred, _, _, z2_pred = self.forward(x1_data, mol_features)
                z2, _, _ = self.encode_x2(x2_data)
                x2_rec = self.decode_x2(z2)
                try:
                    x1_array = torch.cat([x1_array, x1_data], dim=0)
                    x1_rec_array = torch.cat([x1_rec_array, x1_rec], dim=0)
                    x2_array = torch.cat([x2_array, x2_data], dim=0)
                    x2_rec_array = torch.cat([x2_rec_array, x2_rec], dim=0)
                    x2_pred_array = torch.cat([x2_pred_array, x2_pred], dim=0)
                    z2_pred_array = torch.cat([z2_pred_array, z2_pred], dim=0)
                    mid_array = torch.cat([mid_array, mid], dim=0)
                    cid_array = np.concatenate((cid_array, cid), axis=0)
                    sig_array = np.concatenate((sig_array, sig), axis=0)
                except UnboundLocalError:
                    x1_array = x1_data.clone()
                    x1_rec_array = x1_rec.clone()
                    x2_array = x2_data.clone()
                    x2_rec_array = x2_rec.clone()
                    x2_pred_array = x2_pred.clone()
                    z2_pred_array = z2_pred.clone()
                    mid_array = mid.clone()
                    cid_array = cid.copy()
                    sig_array = sig.copy()
        return (
            x1_array.cpu().numpy().astype(float),
            x2_array.cpu().numpy().astype(float),
            x1_rec_array.cpu().numpy().astype(float),
            x2_rec_array.cpu().numpy().astype(float),
            x2_pred_array.cpu().numpy().astype(float),
            z2_pred_array.cpu().numpy().astype(float),
            mid_array.cpu().numpy().astype(float),
            cid_array,
            sig_array,
        )

    def predict_profile_for_x1(self, loader):
        self.eval()
        with torch.no_grad():
            for x1_data, mol_features, mid, cid in loader:
                cid = np.array(list(cid))
                x1_data = x1_data.to(self.dev)
                mol_features = mol_features.to(self.dev)
                _, _, _, x2_pred, _, _, _ = self.forward(x1_data, mol_features)
                try:
                    x2_pred_array = torch.cat([x2_pred_array, x2_pred], dim=0)
                    mid_array = torch.cat([mid_array, mid], dim=0)
                    cid_array = np.concatenate((cid_array, cid), axis=0)
                except UnboundLocalError:
                    x2_pred_array = x2_pred.clone()
                    mid_array = mid.clone()
                    cid_array = cid.copy()
        return x2_pred_array.cpu().numpy().astype(float), mid_array.cpu().numpy().astype(float), cid_array

    def eval_x_reconstruction(self, x1, x1_rec, x2, x2_rec, x2_pred, metrics_func=("pearson",)):
        x1_np = x1.detach().cpu().numpy().astype(float)
        x2_np = x2.detach().cpu().numpy().astype(float)
        x1_rec_np = x1_rec.detach().cpu().numpy().astype(float)
        x2_rec_np = x2_rec.detach().cpu().numpy().astype(float)
        x2_pred_np = x2_pred.detach().cpu().numpy().astype(float)
        deg_np = x2_np - x1_np
        deg_rec_np = x2_rec_np - x1_np
        deg_pred_np = x2_pred_np - x1_np

        metrics_dict = defaultdict(float)
        metrics_dict_ls = defaultdict(list)

        def _to_scalar(metric_val):
            if isinstance(metric_val, (int, float, np.floating)):
                return float(metric_val)
            if isinstance(metric_val, tuple):
                return float(np.mean(metric_val))
            return float(metric_val)

        for m in metrics_func:
            metric_fn = get_metric_func(m)
            if m == "precision100":
                for i in range(x1_np.shape[0]):
                    for prefix, t_arr, p_arr in (
                        ("x1_rec", x1_np[i, :], x1_rec_np[i, :]),
                        ("x2_rec", x2_np[i, :], x2_rec_np[i, :]),
                        ("x2_pred", x2_np[i, :], x2_pred_np[i, :]),
                        ("DEG_rec", deg_np[i, :], deg_rec_np[i, :]),
                        ("DEG_pred", deg_np[i, :], deg_pred_np[i, :]),
                    ):
                        neg_p, pos_p = metric_fn(t_arr, p_arr)
                        metrics_dict[f"{prefix}_negative_precision100"] += neg_p
                        metrics_dict[f"{prefix}_positive_precision100"] += pos_p
                        metrics_dict_ls[f"{prefix}_negative_precision100"].append(neg_p)
                        metrics_dict_ls[f"{prefix}_positive_precision100"].append(pos_p)
            else:
                for i in range(x1_np.shape[0]):
                    x1_val = _to_scalar(metric_fn(x1_np[i, :], x1_rec_np[i, :]))
                    metrics_dict["x1_rec_" + m] += x1_val
                    metrics_dict["x2_rec_" + m] += _to_scalar(metric_fn(x2_np[i, :], x2_rec_np[i, :]))
                    metrics_dict["x2_pred_" + m] += _to_scalar(metric_fn(x2_np[i, :], x2_pred_np[i, :]))
                    metrics_dict["DEG_rec_" + m] += _to_scalar(metric_fn(deg_np[i, :], deg_rec_np[i, :]))
                    metrics_dict["DEG_pred_" + m] += _to_scalar(metric_fn(deg_np[i, :], deg_pred_np[i, :]))
                    metrics_dict_ls["x1_rec_" + m].append(x1_val)
                    metrics_dict_ls["x2_rec_" + m].append(_to_scalar(metric_fn(x2_np[i, :], x2_rec_np[i, :])))
                    metrics_dict_ls["x2_pred_" + m].append(_to_scalar(metric_fn(x2_np[i, :], x2_pred_np[i, :])))
                    metrics_dict_ls["DEG_rec_" + m].append(_to_scalar(metric_fn(deg_np[i, :], deg_rec_np[i, :])))
                    metrics_dict_ls["DEG_pred_" + m].append(_to_scalar(metric_fn(deg_np[i, :], deg_pred_np[i, :])))
        return metrics_dict, metrics_dict_ls
