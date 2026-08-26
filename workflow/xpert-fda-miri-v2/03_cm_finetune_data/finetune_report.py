#!/usr/bin/env python3
"""
Step 03 — finetune-data composition report (HTML).

Generates a self-contained HTML report of the finetune data built by this
step, showing the cell types and drug types it contains:

  finetune_data report
    - collapsible drug list grouped by the author L1000 `moa` (mechanism of
      action, from pert_model/xpert/metadata/l1000_sdst_drugs_info_8432.csv);
      compounds without a recorded moa go under "Unclassified".
    - 2D PCA and UMAP of cell embeddings: the L1000 pretrain aggregated per
      cell line (mean embedding) + every finetune sample, colored by cell
      type with the CM finetune points split by dataset_id.

The cell embedding is the pretrained XPert `cell_Embeddings` output — binned
expression (n_bins=128 over [0, 15.002]) embedded via the pretrained bin
table, summed with the per-gene PPI vector
(PPI_gene_vector_128d.npy -> projected to 256 d), LayerNorm, then mean-pooled
over the 978 landmarks. That is exactly the "binned + PPI embedding" the XPert
model builds for a cell.

Output:
  03_cm_finetune_data/output/finetune_report.html   (self-contained)

Usage:
    python 03_cm_finetune_data/finetune_report.py --data <finetune.h5ad>
        --name <workflow label> [--seed 0]
"""

import argparse
import base64
import io
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import (WORKFLOW_ROOT, XPERT_ROOT, CM_DATA_ROOT, L1000_PRETRAIN_H5AD,  # noqa: E402
                   L1000_FULL_H5AD)

STEP_DIR = WORKFLOW_ROOT / "03_cm_finetune_data"
OUT_DIR = STEP_DIR / "output"
HTML_OUT = OUT_DIR / "finetune_report.html"
DEFAULT_DATA = OUT_DIR / "cm_union.h5ad"
PRETRAIN_PTH = WORKFLOW_ROOT / "04_xpert_pretrain" / "pretrain_mdmt_full_50_epoch.pth"
PPI_PATH = XPERT_ROOT / "metadata" / "PPI_gene_vector_128d.npy"
MOA_TABLE = XPERT_ROOT / "metadata" / "l1000_sdst_drugs_info_8432.csv"
SMI_MAPPING = WORKFLOW_ROOT / "input" / "all_drugs_idx2smi.npy"
SMI_MAPPING = WORKFLOW_ROOT / "input" / "all_drugs_idx2smi.npy"

N_BINS = 128
BIN_MIN, BIN_MAX = 0.0, 15.002
HIDDEN = 256
MAX_L1000_LINE_SAMPLES = 40000   # cap per-line sample budget for the embedding


def load_bin_edges():
    return np.linspace(BIN_MIN, BIN_MAX, N_BINS - 1)


def load_cell_emb_weights():
    """Pretrained cell_Embeddings weights: bin table, PPI (per-gene), LayerNorm."""
    import torch

    ck = torch.load(PRETRAIN_PTH, map_location="cpu")["model_state_dict"]

    def p(k):
        return ck["cell_emb." + k].numpy()

    return {
        "word_emb": p("word_embeddings.weight").astype(np.float64),   # (128, 256)
        "ppi": p("pretrained_embed.weight").astype(np.float64),       # (978, 256)
        "gamma": p("LayerNorm.gamma").astype(np.float64),             # (256,)
        "beta": p("LayerNorm.beta").astype(np.float64),               # (256,)
    }


def cell_embed_chunk(binned, w):
    """Faithful mean-pooled cell embedding for a chunk of binned samples."""
    we = w["word_emb"][binned]                # (n, 978, 256)
    emb = we + w["ppi"][None, :, :]
    u = emb.mean(-1, keepdims=True)
    s = ((emb - u) ** 2).mean(-1, keepdims=True)
    ln = (emb - u) / np.sqrt(s + 1e-12)
    ln = w["gamma"] * ln + w["beta"]
    return ln.mean(axis=1)                    # (n, 256)


def compute_embeddings(adata, w, chunk=512):
    """Per-sample cell embeddings for an AnnData (bins X on the fly)."""
    edges = load_bin_edges()
    out = np.zeros((adata.n_obs, HIDDEN), dtype=np.float32)
    for i in range(0, adata.n_obs, chunk):
        block = adata.X[i:i + chunk]
        if hasattr(block, "toarray"):
            block = block.toarray()
        xb = np.digitize(np.asarray(block).astype(np.float64), edges).astype(np.int64)
        out[i:i + chunk] = cell_embed_chunk(xb, w).astype(np.float32)
    return out


def l1000_per_line_embeddings(w):
    """Mean cell embedding per L1000 cell line (chunked accumulation, all rows)."""
    import scanpy as sc

    a = sc.read_h5ad(L1000_FULL_H5AD, backed="r")
    edges = load_bin_edges()
    sums = defaultdict(lambda: np.zeros(HIDDEN, dtype=np.float64))
    counts = defaultdict(int)
    for i in range(0, a.n_obs, 256):
        idx = slice(i, i + 256)
        block = a.X[idx]
        if hasattr(block, "toarray"):
            block = block.toarray()
        cell = a.obs["cell_iname"][idx]
        xb = np.digitize(np.asarray(block, dtype=np.float32).astype(np.float64),
                         edges).astype(np.int64)
        emb = cell_embed_chunk(xb, w)
        for j, name in enumerate(cell):
            sums[name] += emb[j]
            counts[name] += 1
    return (pd.Series(counts).sort_index(),
            pd.DataFrame({k: v / counts[k] for k, v in sums.items()}).T)


def load_finetune_drugs(adata, manifest=None):
    """Per-drug (pert_idx) summary: name, canonical smiles, n samples, datasets."""
    idx2smi = np.load(SMI_MAPPING, allow_pickle=True).item()
    from rdkit import Chem

    def canon(smi):
        mol = Chem.MolFromSmiles(smi)
        return None if mol is None else Chem.MolToSmiles(mol, isomericSmiles=False)

    smi2idx = {}
    for i, s in idx2smi.items():
        c = canon(s)
        if c is not None:
            smi2idx.setdefault(c, int(i))

    rows = []
    for pi, grp in adata.obs.groupby("pert_idx"):
        pi = int(pi)
        smi = None
        if "smiles" in grp.columns:
            smis = grp["smiles"].dropna().unique()
            smis = [s for s in smis if s]
            smi = smis[0] if smis else None
        rows.append({
            "pert_idx": pi,
            "compound_name": _cmap_name(pi) or grp["compound_name"].iloc[0],
            "smiles": smi or "",
            "n_samples": int(len(grp)),
            "datasets": ",".join(sorted(grp["dataset_id"].astype(str).unique())),
        })
    df = pd.DataFrame(rows).sort_values("pert_idx").reset_index(drop=True)
    return df, smi2idx


def _cmap_name(idx):
    try:
        import scanpy as sc
        a = sc.read_h5ad(L1000_FULL_H5AD, backed="r")
        names = a.obs[a.obs["pert_idx"].astype(int) == idx]["cmap_name"].dropna().unique()
        return str(names[0]) if len(names) else None
    except Exception:
        return None


def moa_by_drug_idx(smi2idx):
    """drug_idx -> author moa(s) from the author L1000 drug-info table."""
    from rdkit import Chem

    def canon(smi):
        mol = Chem.MolFromSmiles(smi)
        return None if mol is None else Chem.MolToSmiles(mol, isomericSmiles=False)

    t = pd.read_csv(MOA_TABLE)
    out = defaultdict(list)
    n_matched = 0
    for _, r in t.iterrows():
        if pd.isna(r["moa"]) or pd.isna(r["canonical_smiles"]):
            continue
        c = canon(str(r["canonical_smiles"]))
        if c is not None and c in smi2idx:
            out[smi2idx[c]].append(str(r["moa"]))
            n_matched += 1
    for k in out:
        seen, uniq = set(), []
        for m in out[k]:
            for part in m.replace("|", ",").split(","):
                part = part.strip()
                if part and part not in seen:
                    seen.add(part)
                    uniq.append(part)
        out[k] = uniq
    print(f"  moa table: {n_matched} rows joined to {len(out)} drug indices")
    return dict(out)


def build_drug_html(drugs, moa_map):
    """Collapsible drug list grouped by moa (details/summary)."""
    grouped = defaultdict(list)
    for _, r in drugs.iterrows():
        moas = moa_map.get(int(r["pert_idx"]))
        key = "Unclassified" if not moas else "; ".join(moas)
        grouped[key].append(r)
    n_total = int(drugs["n_samples"].sum())
    blocks = []
    for moa in sorted(grouped, key=lambda k: (k == "Unclassified", k.lower())):
        rows = grouped[moa]
        n_drugs = len(rows)
        n_samples = int(sum(r["n_samples"] for r in rows))
        items = "\n".join(
            f'      <li><code>idx {int(r["pert_idx"])}</code> — '
            f'<b>{r["compound_name"]}</b> · {int(r["n_samples"])} samples'
            f' · {r["datasets"]}</li>'
            for r in rows)
        blocks.append(
            f'    <details>\n'
            f'      <summary><b>{moa}</b> '
            f'<span class="muted">({n_drugs} drugs · {n_samples} samples)</span></summary>\n'
            f'      <ul>\n{items}\n      </ul>\n'
            f'    </details>\n')
    return f'<h2>Drugs in the finetune data ({len(drugs)} drugs · {n_total} samples)</h2>\n' \
           f'  <p class="muted">Grouped by the author L1000 <code>moa</code>. ' \
           f'Drugs without a recorded moa are under <i>Unclassified</i>.</p>\n' + "".join(blocks)


def png_base64(fig):
    import matplotlib.pyplot as plt

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def make_embedding_figs(emb_df):
    """Two 2D plots (PCA, UMAP) of the cell embeddings, colored by group."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from sklearn.decomposition import PCA
    import umap

    labels = emb_df["group"].values
    X = emb_df[[f"e{i}" for i in range(HIDDEN)]].values.astype(np.float64)

    fig, ax = plt.subplots(figsize=(8, 6))
    pca = PCA(n_components=2, random_state=0)
    xy = pca.fit_transform(X)
    _scatter(ax, xy, labels)
    ax.set_title(f"Cell-embedding PCA — {emb_df.attrs.get('workflow', '')}")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
    pca_png = png_base64(fig)

    fig, ax = plt.subplots(figsize=(8, 6))
    red = umap.UMAP(n_components=2, random_state=0, n_neighbors=30, min_dist=0.1)
    xy = red.fit_transform(X)
    _scatter(ax, xy, labels)
    ax.set_title(f"Cell-embedding UMAP — {emb_df.attrs.get('workflow', '')}")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    umap_png = png_base64(fig)

    return pca_png, umap_png, pca, red


def _scatter(ax, xy, labels):
    import matplotlib.pyplot as plt

    uniq = pd.unique(labels)
    n_l1000 = sum(1 for l in labels if not l.startswith("CM-"))
    cmap = plt.get_cmap("tab20")
    colors = {l: cmap(i % 20) for i, l in enumerate(uniq)}
    legend_l1000 = []
    for l in uniq:
        mask = labels == l
        if l.startswith("CM-"):
            marker, size = "s", 12
        else:
            marker, size = "o", 26
        ax.scatter(xy[mask, 0], xy[mask, 1], s=size, c=[colors[l]],
                   marker=marker, alpha=0.85, edgecolors="none", label=l)
    ax.set_aspect("auto")
    handles, labels_used = ax.get_legend_handles_labels()
    # collapse L1000 lines into one legend entry + keep CM entries
    keep = [i for i, l in enumerate(labels_used) if l.startswith("CM-")]
    keep += [i for i, l in enumerate(labels_used) if not l.startswith("CM-")][:1]
    handles = [handles[i] for i in keep]
    labels_used = [labels_used[i] for i in keep]
    if n_l1000:
        ax.legend(handles, labels_used + [f"{n_l1000} L1000 cell lines"],
                  fontsize="small", markerscale=0.8, loc="best")


def main():
    parser = argparse.ArgumentParser(description="Finetune-data composition report")
    parser.add_argument("--data", default=str(DEFAULT_DATA),
                        help="Finetune h5ad to report on")
    parser.add_argument("--name", default="", help="Workflow label for titles")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    import anndata as ad
    import matplotlib
    matplotlib.use("Agg")

    data_path = Path(args.data)
    if not data_path.exists():
        print(f"Missing finetune h5ad: {data_path} (run the 03 build first)")
        sys.exit(1)

    print(f"Report: {data_path}")
    print("[1/4] load finetune h5ad + pretrained cell embedding weights")
    adata = ad.read_h5ad(data_path)
    w = load_cell_emb_weights()

    print("[2/4] cell embeddings (L1000 per line + finetune samples)")
    line_counts, line_emb = l1000_per_line_embeddings(w)
    fin_emb = compute_embeddings(adata, w)
    print(f"  L1000: {len(line_emb)} cell lines | finetune: {fin_emb.shape[0]} samples")

    fin_labels = np.array([
        f"CM-{ds}" for ds in adata.obs["dataset_id"].astype(str).values])
    line_names = [str(n) for n in line_emb.index]
    emb_df = pd.DataFrame(line_emb.values.astype(np.float64),
                          columns=[f"e{i}" for i in range(HIDDEN)])
    emb_df["group"] = line_names
    fin_df = pd.DataFrame(fin_emb.astype(np.float64),
                          columns=[f"e{i}" for i in range(HIDDEN)])
    fin_df["group"] = fin_labels
    emb_df = pd.concat([emb_df, fin_df], ignore_index=True)
    emb_df.attrs["workflow"] = args.name

    print("[3/4] 2D PCA + UMAP")
    pca_png, umap_png, _, _ = make_embedding_figs(emb_df)

    print("[4/4] drug list by moa + HTML")
    drugs, smi2idx = load_finetune_drugs(adata)
    moa_map = moa_by_drug_idx(smi2idx)
    drug_html = build_drug_html(drugs, moa_map)

    n_cm_types = adata.obs["dataset_id"].nunique() if "dataset_id" in adata.obs else 1
    cell_summary = (
        f"<li><b>{adata.n_obs}</b> finetune samples, <b>{adata.n_vars}</b> landmark "
        f"genes, <b>{adata.obs['pert_idx'].nunique()}</b> drugs</li>"
        f"<li>Cell type(s): {', '.join(sorted(adata.obs['cell_iname'].astype(str).unique()))}"
        f" ({n_cm_types} dataset groups)</li>"
        f"<li>L1000 reference: {len(line_emb)} cell lines, {int(line_counts.sum())} "
        f"profiles (mean embedding per line)</li>")

    html = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>Finetune-data composition — {args.name or data_path.parent.name}</title>
<style>
  body {{ font-family: -apple-system, 'Segoe UI', Roboto, sans-serif; margin: 2rem auto;
         max-width: 1100px; color: #222; }}
  h1 {{ border-bottom: 2px solid #333; padding-bottom: .3rem; }}
  h2 {{ margin-top: 2rem; }}
  .muted {{ color: #666; font-weight: normal; }}
  details {{ margin: .3rem 0; }}
  summary {{ cursor: pointer; }}
  ul {{ margin: .3rem 0 1rem 1.2rem; }}
  li {{ margin: .15rem 0; }}
  code {{ background: #f2f2f2; padding: 0 .25rem; border-radius: 3px; }}
  img {{ max-width: 100%; height: auto; border: 1px solid #ddd; margin: .5rem 0; }}
  .row {{ display: flex; gap: 1rem; flex-wrap: wrap; }}
  .col {{ flex: 1 1 45%; }}
  .legend {{ font-size: .85rem; color: #555; }}
</style></head><body>
<h1>Finetune-data composition — {args.name or data_path.parent.name}</h1>
<p class="muted">Generated by <code>03_cm_finetune_data/finetune_report.py</code>. '
 'Cell embeddings are the pretrained XPert <code>cell_Embeddings</code> output '
 '(binned expression + PPI gene vectors, mean-pooled over the 978 landmarks).</p>
<h2>Cells</h2>
<ul>{cell_summary}</ul>
<div class="row">
  <div class="col"><h3>PCA</h3><img src="data:image/png;base64,{pca_png}" alt="PCA"></div>
  <div class="col"><h3>UMAP</h3><img src="data:image/png;base64,{umap_png}" alt="UMAP"></div>
</div>
<p class="legend">Squares = finetune samples (grouped by dataset), circles = L1000 cell lines '
 '(one point per line, mean embedding). Finetune samples are colored by dataset; the L1000 '
 'lines are collapsed into a single legend entry.</p>
{drug_html}
</body></html>"""

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    HTML_OUT.write_text(html)
    print(f"Report -> {HTML_OUT} ({HTML_OUT.stat().st_size / 1024:.0f} KiB)")


if __name__ == "__main__":
    main()
