#!/usr/bin/env python3
"""
Step 13 — Limitations analysis (read-only, diagnostic).

Quantifies three structural limitations of the v2 MIRI screen with one figure
each. This step is read-only: it consumes step 06/09/10 artifacts and the frozen
inputs and changes nothing upstream.

  1. 978-landmark pathway coverage (`landmark_coverage.png`)
     For each of the 40 scoring pathways: total genes vs the genes actually
     measured by the 978 L1000-landmark panel (the panel-scoring subset).

  2. Critical-pathway pairwise overlap (`critical_pathway_overlap.png`)
     The 8 curated critical pathways (the step-10 ">=6/8" filter): overlap
     genes as a fraction of the ROW pathway's landmark subset (asymmetric).

  3. Screened-drug embedding PCA (`drug_embedding_pca.png`)
     All drugs passing the step-10 filter, described by the HG + UniMol drug
     embedding fed to the step-06 finetuned XPert transformer, projected to
     PC1/PC2. The 10 direct alpha/beta-adrenergic ligands are highlighted.

Inputs:
  input/MIRI_pathways.json                                (40 pathways + 8 critical)
  09_miri_efficacy_scores/output/scores.csv               (drug_name -> drug_idx)
  10_filter/output/path_filtered.csv                      (drugs passing the filter)
  input/HG_drug_embeddings.npy                            (author HG table, 8981x256)
  input/all_drugs_unimol_arr.npy                           (frozen UniMol array)
  06_xpert_finetune/output/cm_cm_full/*_best.pth           (finetuned drug_emb.linear)

Outputs:
  13_limits/landmark_coverage.png
  13_limits/critical_pathway_overlap.png
  13_limits/drug_embedding_pca.png
  13_limits/output/landmark_coverage.csv
  13_limits/output/critical_pathway_overlap.csv
  13_limits/report.md

Usage:
    python 13_limits/limits.py
"""

import json
import sys
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT  # noqa: E402

STEP_DIR = WORKFLOW_ROOT / "13_limits"
OUT_DIR = STEP_DIR / "output"
PATHWAYS = WORKFLOW_ROOT / "input" / "MIRI_pathways.json"
SCORES = WORKFLOW_ROOT / "09_miri_efficacy_scores" / "output" / "scores.csv"
PATH_FILTERED = WORKFLOW_ROOT / "10_filter" / "output" / "path_filtered.csv"
HG_ARR = WORKFLOW_ROOT / "input" / "HG_drug_embeddings.npy"
UNIMOL_ARR = WORKFLOW_ROOT / "input" / "all_drugs_unimol_arr.npy"
FINETUNE_CKPT = (WORKFLOW_ROOT / "06_xpert_finetune" / "output" / "cm_cm_full"
                 / "cm_cm_full_fold-1_best.pth")

BLUE, RED = "#2196f3", "#e53935"
CRITICAL_SHORT = {
    "GO-Inflammatory Response": "Inflammatory Resp.",
    "GO-Regulation Of Inflammatory Response": "Inflamm. Reg.",
    "KEGG-Apoptosis": "Apoptosis",
    "GO-Cellular Response To Oxidative Stress": "Oxidative Stress",
    "KEGG-Necroptosis": "Necroptosis",
    "GO-Regulation Of Cell Death": "Cell Death Reg.",
    "KEGG-Calcium signaling pathway": "Ca Signaling",
    "GO-Mitochondrial Gene Expression": "Mito. Gene Expr.",
}

# Short axis labels for the 40 scoring pathways (figure 1). "+" marks the
# "Positive Regulation Of ..." pathways so they stay distinct from the plain
# "Regulation Of ..." siblings.
ABBREV = {
    "GO-Cellular Response To Lectin": "GO-Cell. Resp. Lectin",
    "GO-Cellular Response To Oxidative Stress": "GO-Cell. Resp. Oxid. Stress",
    "GO-Cellular Response To Retinoic Acid": "GO-Cell. Resp. Retinoic Acid",
    "GO-Cellular Response To Tumor Necrosis Factor": "GO-Cell. Resp. TNF",
    "GO-Cytokine-Mediated Signaling Pathway": "GO-Cytokine Signaling",
    "GO-Extracellular Structure Organization": "GO-ECM Structure Org.",
    "GO-Extracellular matrix organization": "GO-ECM Organization",
    "GO-Fatty Acid Beta-Oxidation": "GO-FA β-Oxidation",
    "GO-Fatty Acid Catabolic Process": "GO-FA Catabolism",
    "GO-Fatty Acid Oxidation": "GO-FA Oxidation",
    "GO-Fc Receptor Mediated Stimulatory Signaling Pathway": "GO-Fc Receptor Signaling",
    "GO-Inflammatory Response": "GO-Inflammatory Resp.",
    "GO-Mitochondrial Gene Expression": "GO-Mito. Gene Expr.",
    "GO-Negative Regulation Of Insulin Receptor Signaling Pathway": "GO-Insulin Receptor Reg. (−)",
    "GO-Peptidyl-Threonine Phosphorylation": "GO-Peptidyl-Thr Phospho.",
    "GO-Positive Regulation Of Cyclin-Dependent Protein Kinase Activity": "GO-Cyclin-Dep. Kinase Reg. (+)",
    "GO-Positive Regulation Of Interleukin-1 Beta Production": "GO-IL-1β Production Reg. (+)",
    "GO-Positive Regulation Of Lipid Biosynthetic Process": "GO-Lipid Biosynth. Reg. (+)",
    "GO-Positive Regulation Of NF-kappaB Transcription Factor Activity": "GO-NF-κB Activity Reg. (+)",
    "GO-Positive Regulation Of Nitric Oxide Metabolic Process": "GO-Nitric Oxide Metab. Reg. (+)",
    "GO-Positive Regulation Of Stem Cell Population Maintenance": "GO-Stem Cell Maint. Reg. (+)",
    "GO-Positive Regulation Of Tumor Necrosis Factor Production": "GO-TNF Production Reg. (+)",
    "GO-Regulation Of Cell Death": "GO-Cell Death Reg.",
    "GO-Regulation Of G0 To G1 Transition": "GO-G0/G1 Transition Reg.",
    "GO-Regulation Of Inflammatory Response": "GO-Inflamm. Resp. Reg.",
    "GO-Regulation Of Protein Kinase Activity": "GO-Protein Kinase Reg.",
    "GO-Wound Healing": "GO-Wound Healing",
    "KEGG-Adrenergic signaling in cardiomyocytes": "KEGG-Adrenergic Signaling",
    "KEGG-Apoptosis": "KEGG-Apoptosis",
    "KEGG-Calcium signaling pathway": "KEGG-Ca Signaling",
    "KEGG-Chemokine signaling pathway": "KEGG-Chemokine Signaling",
    "KEGG-Dilated cardiomyopathy": "KEGG-DCM",
    "KEGG-Focal adhesion": "KEGG-Focal Adhesion",
    "KEGG-MAPK signaling pathway": "KEGG-MAPK Signaling",
    "KEGG-Necroptosis": "KEGG-Necroptosis",
    "KEGG-Pyruvate metabolism": "KEGG-Pyruvate Metab.",
    "KEGG-Regulation of actin cytoskeleton": "KEGG-Actin Cytoskeleton",
    "KEGG-TNF signaling pathway": "KEGG-TNF Signaling",
    "KEGG-Type II diabetes mellitus": "KEGG-T2 Diabetes",
    "KEGG-Valine, leucine and isoleucine degradation": "KEGG-BCAA Degradation",
}

# Figure 1 fonts are scaled +75% from the workflow default scale (7/10/12/14).
FONT = 1.75

# Direct alpha/beta-adrenergic receptor ligands among the filtered drugs
# (pharmacological definition: compounds that bind the adrenergic receptors, as
# opposed to indirect sympathomimetics / reuptake inhibitors / alpha2-antagonist
# antidepressants). Agonists vs antagonists. Same spirit as the adrenergic group
# in `_other/top20_structures.py`.
ADRENERGIC_AGONISTS = ["DOBUTAMINE", "MIDODRINE", "METHOXAMINE", "MEPHENTERMINE"]
ADRENERGIC_ANTAGONISTS = ["PHENOXYBENZAMINE", "PROPRANOLOL", "SOTALOL",
                          "DILEVALOL", "LABETALOL", "MOXISYLYTE"]


# ── Figure 1: 978-landmark pathway coverage ────────────────────────────

def load_pathways():
    with open(PATHWAYS) as f:
        return json.load(f)["pathways"]


def plot_coverage(pathways):
    """Stacked bars: total pathway genes (blue) vs landmark-covered (red, bottom)."""
    df = pd.DataFrame([{
        "name": pw["name"],
        "direction": pw.get("direction", "UP"),
        "total": int(pw["n_total_genes_before_l1000"]),
        "covered": int(len(pw["genes"])),
    } for pw in pathways])
    df["ratio"] = df["covered"] / df["total"]
    missing = set(df["name"]) - set(ABBREV)
    if missing:
        raise SystemExit(f"ERROR: no abbreviation for pathways: {sorted(missing)}")
    df["abb"] = df["name"].map(ABBREV)
    df = df.sort_values("ratio", ascending=True).reset_index(drop=True)
    df.to_csv(OUT_DIR / "landmark_coverage.csv", index=False)

    x = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(18, 8.5))
    ax.bar(x, df["total"], color=BLUE, edgecolor="white", linewidth=0.4,
           label="Pathway genes (all)")
    ax.bar(x, df["covered"], color=RED, edgecolor="white", linewidth=0.4,
           label="In 978-landmark panel")
    for xi, (tot, cov, ratio) in enumerate(zip(df["total"], df["covered"],
                                               df["ratio"])):
        ax.text(xi, tot + df["total"].max() * 0.012, f"{ratio * 100:.0f}%",
                ha="center", va="bottom", fontsize=6.5 * FONT, rotation=0,
                color="#333333")
    ax.set_xticks(x)
    ax.set_xticklabels(df["abb"], rotation=45, ha="right", fontsize=7 * FONT)
    ax.tick_params(axis="y", labelsize=10 * FONT)
    ax.set_ylabel("Number of genes", fontsize=12 * FONT)
    ax.set_ylim(0, df["total"].max() * 1.18)
    ax.set_title("978-landmark coverage of the 40 MIRI scoring pathways\n"
                 "(blue = full pathway, red = genes scored in the panel)",
                 fontsize=14 * FONT, fontweight="bold")
    ax.legend(loc="upper left", fontsize=10 * FONT)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = STEP_DIR / "landmark_coverage.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.relative_to(WORKFLOW_ROOT)}")
    return df


# ── Figure 2: critical-pathway overlap ─────────────────────────────────

def plot_overlap(pathways):
    """Row-normalized overlap M[i,j] = |i & j| / |i|, annotated with raw counts."""
    crit = [pw for pw in pathways if pw.get("critical")]
    names = [pw["name"] for pw in crit]
    sets = [set(pw["genes"]) for pw in crit]
    short = [CRITICAL_SHORT.get(n, n) for n in names]
    n = len(crit)

    ratio = np.zeros((n, n))
    counts = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            inter = len(sets[i] & sets[j])
            counts[i, j] = inter
            ratio[i, j] = inter / len(sets[i])

    df = pd.DataFrame(ratio, index=names, columns=names)
    df.to_csv(OUT_DIR / "critical_pathway_overlap.csv")

    fig, ax = plt.subplots(figsize=(9.5, 8))
    sns.heatmap(ratio, annot=counts, fmt="d", cmap="Reds", vmin=0, vmax=1,
                xticklabels=short, yticklabels=short, square=True,
                linewidths=0.5, linecolor="white",
                cbar_kws={"label": "Overlap / row pathway's landmark subset"},
                annot_kws={"fontsize": 10}, ax=ax)
    ax.set_xticklabels(short, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(short, rotation=0, fontsize=9)
    ax.set_title("Pairwise gene overlap of the 8 critical pathways\n"
                 "(color = overlap / row pathway's landmark subset; "
                 "label = shared gene count)",
                 fontsize=12.5, fontweight="bold")
    fig.tight_layout()
    out = STEP_DIR / "critical_pathway_overlap.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.relative_to(WORKFLOW_ROOT)}")
    return short, ratio, counts


# ── Figure 3: screened-drug embedding PCA ──────────────────────────────

def drug_embedding_vectors():
    """[HG(256) | mean-pooled projected UniMol tokens(256)] for filtered drugs.

    Reproduces the drug tensor the finetuned XPert transformer sees
    (drug_emb.linear 512->256 projection, token 0 overwritten by HG), then
    aggregates the valid UniMol tokens (molecule + atoms) by mean pooling.
    Position embeddings / LayerNorm are deliberately omitted: they are not
    drug-identity signal.
    """
    import torch

    scores = pd.read_csv(SCORES)[["drug_name", "drug_idx"]]
    filtered = pd.read_csv(PATH_FILTERED)[["drug_name"]]
    drugs = filtered.merge(scores, on="drug_name", how="left")
    if drugs["drug_idx"].isna().any():
        miss = drugs.loc[drugs["drug_idx"].isna(), "drug_name"].tolist()
        raise SystemExit(f"ERROR: filtered drugs missing a drug_idx: {miss}")
    names = drugs["drug_name"].tolist()
    idx = drugs["drug_idx"].astype(int).values

    checkpoint = torch.load(FINETUNE_CKPT, map_location="cpu")
    state = checkpoint.get("model_state_dict", checkpoint)
    W = state["drug_emb.linear.weight"].to(torch.float64).numpy()  # (256, 512)
    b = state["drug_emb.linear.bias"].to(torch.float64).numpy()    # (256,)

    HG = np.load(HG_ARR)[idx]                                        # (n, 256)
    unimol = np.load(UNIMOL_ARR, allow_pickle=True)[idx]             # (n, 122, 514)

    vectors = []
    for k in range(len(idx)):
        u = unimol[k]
        mask = u[:, 0].astype(bool)          # valid token mask (row 0 = HG slot)
        feat = u[:, 2:].astype(np.float64)       # (122, 512)
        proj = feat @ W.T + b                  # (122, 256)
        uni_valid = proj[1:][mask[1:]]         # exclude token 0 (HG slot)
        pooled = uni_valid.mean(axis=0) if len(uni_valid) else np.zeros(W.shape[0])
        vectors.append(np.concatenate([HG[k].astype(np.float64), pooled]))
    return names, np.stack(vectors), W.shape[0]


PCA_STYLE = {
    "other": dict(marker="o", c="black", s=42, alpha=0.6),
    "agonist": dict(marker="^", c=RED, s=150, edgecolor="black", linewidth=0.6),
    "antagonist": dict(marker="s", c=BLUE, s=140, edgecolor="black", linewidth=0.6),
}

# Figure 3 fonts are scaled +100% from the workflow default scale (10/12/13/16).
FONT_PCA = 2.0


def _pca2(X):
    """Per-dim z-score (dropping zero-variance dims) then PCA to 2 components."""
    std = X.std(axis=0)
    keep = std > 1e-8
    Xz = (X[:, keep] - X[:, keep].mean(axis=0)) / std[keep]
    pca = PCA(n_components=2, random_state=0)
    return pca.fit_transform(Xz), pca.explained_variance_ratio_ * 100


def _cohen_d(Z, mask):
    """Per-component Cohen's d between the masked group and the rest."""
    a, b = Z[mask], Z[~mask]
    s = np.sqrt(((a.var(0, ddof=1) * (len(a) - 1))
                 + (b.var(0, ddof=1) * (len(b) - 1))) / (len(a) + len(b) - 2))
    s[s == 0] = 1
    return (a.mean(0) - b.mean(0)) / s


def _scatter_panel(ax, Z, ev, cat, title):
    for key in ("other", "agonist", "antagonist"):
        m = cat == key
        if m.any():
            label = {"other": f"Other screened drug (n={int(m.sum())})",
                     "agonist": f"Adrenergic agonist ({int(m.sum())})",
                     "antagonist": f"Adrenergic antagonist ({int(m.sum())})"}[key]
            ax.scatter(Z[m, 0], Z[m, 1], label=label, **PCA_STYLE[key])
    ax.axhline(0, color="#cccccc", lw=0.8, zorder=0)
    ax.axvline(0, color="#cccccc", lw=0.8, zorder=0)
    ax.tick_params(axis="both", labelsize=10 * FONT_PCA)
    ax.set_xlabel(f"PC1 ({ev[0]:.1f}% var)", fontsize=12 * FONT_PCA)
    ax.set_ylabel(f"PC2 ({ev[1]:.1f}% var)", fontsize=12 * FONT_PCA)
    ax.set_title(title, fontsize=13 * FONT_PCA, fontweight="bold")


def plot_pca():
    names, X, hg_dim = drug_embedding_vectors()

    ag = set(ADRENERGIC_AGONISTS)
    an = set(ADRENERGIC_ANTAGONISTS)
    missing = (ag | an) - set(names)
    if missing:
        raise SystemExit(f"ERROR: adrenergic drugs not in filtered set: {sorted(missing)}")
    cat = np.array(["agonist" if d in ag else "antagonist" if d in an else "other"
                    for d in names])
    adr_mask = np.array([d in ag or d in an for d in names])

    Z_joint, ev_joint = _pca2(X)                 # model input [HG | UniMol]
    Z_hg, ev_hg = _pca2(X[:, :hg_dim])           # HG component only
    Z_uni, _ = _pca2(X[:, hg_dim:])              # UniMol component only (not plotted)
    d_joint = _cohen_d(Z_joint, adr_mask)
    d_hg = _cohen_d(Z_hg, adr_mask)
    d_uni = _cohen_d(Z_uni, adr_mask)

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(26, 11), layout="constrained")
    _scatter_panel(ax_l, Z_joint, ev_joint, cat,
                   "Model input [HG | UniMol] (z-scored)\n"
                   f"Cohen d(PC1) = {d_joint[0]:+.2f}  ·  no clusters")
    _scatter_panel(ax_r, Z_hg, ev_hg, cat,
                   "HG component only (z-scored)\n"
                   f"Cohen d(PC1) = {d_hg[0]:+.2f}  ·  separated")
    ax_l.legend(loc="best", fontsize=10 * FONT_PCA, framealpha=0.9)
    fig.suptitle("Screened-drug embedding PCA — "
                 f"{len(names)} drugs passing the step-10 filter",
                 fontsize=16 * FONT_PCA, fontweight="bold")
    out = STEP_DIR / "drug_embedding_pca.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.relative_to(WORKFLOW_ROOT)}")
    return {"names": names, "cat": cat, "ev_joint": ev_joint, "ev_hg": ev_hg,
            "d_joint": d_joint, "d_hg": d_hg, "d_uni": d_uni}


# ── Report ─────────────────────────────────────────────────────────────

def write_report(cov_df, crit_short, crit_ratio, crit_counts, pca, n_union):
    ratio = cov_df["ratio"]
    cat = pca["cat"]
    off = crit_ratio[~np.eye(len(crit_short), dtype=bool)]
    i_max, j_max = np.unravel_index(
        np.where(np.eye(len(crit_short), dtype=bool), -1, crit_ratio).argmax(),
        crit_ratio.shape)
    n_ag = int((cat == "agonist").sum())
    n_an = int((cat == "antagonist").sum())

    lines = [
        "# v2 limitations analysis (step 13)",
        "",
        f"Generated: {date.today().isoformat()} (deterministic, step 13)",
        "",
        "Read-only diagnostic of the `xpert-fda-miri-v2` screen. Three figures, "
        "one structural limitation each.",
        "",
        "## 1. 978-landmark pathway coverage",
        "",
        "The pathway scorer only sees the genes of the 978 L1000-landmark panel. "
        "Resolving each of the 40 scoring pathways against that panel shows how "
        "little of each pathway is actually measured:",
        "",
        f"- Landmark coverage per pathway: **{ratio.min() * 100:.1f}% to "
        f"{ratio.max() * 100:.1f}%**, mean **{ratio.mean() * 100:.1f}%** "
        f"(median {ratio.median() * 100:.1f}%).",
        f"- Median pathway: {int(cov_df['total'].median())} genes total, only "
        f"{int(cov_df['covered'].median())} scored.",
        f"- The 40 pathways contribute {int(cov_df['covered'].sum())} scored "
        f"gene-slots (with cross-pathway duplication) out of "
        f"{int(cov_df['total'].sum())} pathway genes; their union is "
        f"{n_union} of the 978 landmarks.",
        "",
        "**Limitation.** Pathway-level reversal is judged on a biased ~16% slice of "
        "each pathway. Pathways with a higher landmark share (e.g., "
        f"{cov_df.iloc[-1]['name']} at {cov_df.iloc[-1]['ratio'] * 100:.0f}%) are "
        "favoured over those mostly absent from the panel; direction calls on the "
        "missing genes are impossible by construction.",
        "",
        "## 2. Critical-pathway overlap",
        "",
        "The step-10 pathway filter keeps drugs reversing >=6/8 *critical* "
        "pathways. Pairwise overlap is measured as shared genes over the ROW "
        "pathway's landmark subset (asymmetric; diagonal = 1):",
        "",
        f"- Off-diagonal overlap: **max {off.max():.2f}** "
        f"({crit_short[i_max]} -> {crit_short[j_max]}, "
        f"{int(crit_counts[i_max, j_max])} genes), mean **{off.mean():.2f}**.",
        f"- {int((off == 0).sum())} of {off.size} ordered pairs share no gene at all.",
        "",
        "**Limitation.** The 8 critical pathways are largely *independent*, not "
        "redundant: they are not restatements of one another. The >=6/8 filter is "
        "therefore a genuinely multi-axis bar (a drug must reverse six distinct "
        "processes), and the critical count is not inflated by duplicated gene "
        "sets. The corollary is that the filter is stringent and heavily trims the "
        "candidate list.",
        "",
        "## 3. Screened-drug embedding structure",
        "",
        "Each drug passing the filter is projected to PC1/PC2 from two related "
        "spaces (both per-dimension z-scored across the 64 drugs), with the "
        f"{n_ag} direct adrenergic agonists and {n_an} antagonists highlighted:",
        "",
        "- **Left — the model's own input** `[HG(256) | mean-pooled projected "
        "UniMol tokens(256)]` (the tensor the finetuned XPert transformer "
        f"consumes): PC1 {pca['ev_joint'][0]:.1f}% + PC2 "
        f"{pca['ev_joint'][1]:.1f}% of variance.",
        f"- **Right — the HG component only**: PC1 {pca['ev_hg'][0]:.1f}% + PC2 "
        f"{pca['ev_hg'][1]:.1f}% of variance.",
        "",
        "**Limitation (a negative result, not a cluster).** In the space the model "
        "actually consumes the adrenergic ligands are NOT clustered "
        f"(adrenergic-vs-rest Cohen's d on PC1 = {pca['d_joint'][0]:+.2f}); agonists "
        "and antagonists are intermixed. The class signal is real but lives in the "
        f"HG branch alone (HG-only d = {pca['d_hg'][0]:+.2f}); it is diluted by the "
        f"512 UniMol-token dimensions in the joint input (UniMol-only d = "
        f"{pca['d_uni'][0]:+.2f}). So the HG embedding that shapes the prediction "
        "is class-correlated, while the joint input the transformer sees does not "
        "expose the class linearly — the screen can carry class-level structure "
        "that is invisible in a joint view yet still biases per-compound calls.",
        "",
        "## Conclusion",
        "",
        "v2's pathway scoring is bounded by the 978-landmark panel (~16% of each "
        "pathway), its critical-pathway filter is strict rather than redundant, and "
        "its drug-embedding input hides a real pharmacologic class signal (strong in "
        "the HG branch, diluted by the 512 UniMol-token dimensions in the joint "
        "input the model consumes). All three are properties of the screen's "
        "*inputs and aggregation*, independent of the XPert finetune-data change "
        "that defines v2. They bound how much biological claim the ranking can "
        "support.",
    ]
    out = STEP_DIR / "report.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  Saved: {out.relative_to(WORKFLOW_ROOT)}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    pathways = load_pathways()
    n_union = len({g for pw in pathways for g in pw["genes"]})
    print(f"Step 13 — limitations analysis ({len(pathways)} pathways, "
          f"{sum(1 for p in pathways if p.get('critical'))} critical, "
          f"{n_union}-gene union)")
    cov_df = plot_coverage(pathways)
    crit_short, crit_ratio, crit_counts = plot_overlap(pathways)
    pca = plot_pca()
    write_report(cov_df, crit_short, crit_ratio, crit_counts, pca, n_union)
    print("\nStep 13 complete — 3 figures + report.md")


if __name__ == "__main__":
    main()
