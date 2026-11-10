"""Render the two exploratory figures from the published score tables."""

import os

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.patheffects as effects
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


COLORS = {
    "no_change": "#757575", "delta_mean": "#a3a3a3", "XPert_base": "#ffffff",
    "XPert_G2000": "#a6cee3", "XPert_G5000": "#66a6d3",
    "XPert_G10000": "#327db8", "XPert_G18164": "#11518b",
    "XPert_MLP512": "#a9d8a4", "XPert_MLP1024": "#72b674",
    "XPert_MLP2048": "#38944c", "XPert_MLP4096": "#15623b",
}
PANELS = (("auprc_full", "AUPRC@full"), ("auprc_g", "AUPRC@G"), ("pcc", "PCC"),
          ("spearman", "Spearman"), ("rmse", "RMSE"), ("r2", "R²"),
          ("training_minutes", "Training minutes"))
TITLE = "XPert Round 1 • exploratory tenv test • fold 0, seed 17"
OUTLINE = [effects.Stroke(linewidth=3.6, foreground="#333333"), effects.Normal()]


def display_name(model):
    return model.replace("XPert_", "")


def draw_metrics(metrics, path):
    fig, axes = plt.subplots(2, 4, figsize=(18.5, 10.4), constrained_layout=True)
    locations = np.arange(len(metrics))
    labels = metrics.model.map(display_name)
    for axis, (metric, title) in zip(axes.flat, PANELS):
        axis.axvline(0, linewidth=.6, color="#777777")
        for row, item in enumerate(metrics.to_dict("records")):
            value = item[metric]
            if pd.notna(value):
                model = item["model"]
                axis.barh(row, value, height=.64, color=COLORS[model], edgecolor="#333333", linewidth=.65)
                text = f"{value:.2f}" if metric == "training_minutes" else f"{value:.3f}"
                if metric != "training_minutes":
                    text += f"   {int(item[f'{metric}_defined_groups'])}/{int(item['groups'])}"
                if item["status"] == "incomplete":
                    text += "   incomplete"
            else:
                text = "incomplete" if item["status"] == "incomplete" else "undefined"
            axis.text(.98, row, text, transform=axis.get_yaxis_transform(), va="center", ha="right",
                      fontsize=7.7, bbox={"facecolor": "white", "edgecolor": "none", "alpha": .78, "pad": .8})
        axis_labels = labels
        if metric == "auprc_full":
            axis_labels = [f"{label} ({coverage:.1%})" for label, coverage in zip(labels, metrics.cm_coverage, strict=True)]
        axis.set_yticks(locations, axis_labels)
        axis.set_ylim(len(metrics) - .48, -.55)
        axis.grid(axis="x", alpha=.24)
        axis.set_axisbelow(True)
        axis.set_title(title, fontsize=13, weight="semibold")
        axis.tick_params(labelsize=8.5)
        axis.spines[["top", "right"]].set_visible(False)
    notes = axes.flat[-1]
    notes.axis("off")
    size = int(metrics.shared_gene_count.iloc[0])
    provisional = bool(metrics.continuous_provisional.iloc[0])
    notes.text(.05, .87, "Metric definitions", fontsize=13, weight="semibold", transform=notes.transAxes)
    notes.text(.05, .72,
               f"Continuous scores use {size} shared genes.\n"
               + ("These scores are provisional.\n" if provisional else "")
               + "Each metric shows defined / total groups.\n"
               + "Undefined means zero defined groups.\n"
               + "AUPRC@full labels show CM coverage.\n"
               + "Training time covers both fit stages.",
               va="top", fontsize=10, linespacing=1.55, transform=notes.transAxes)
    fig.suptitle(TITLE, fontsize=17, weight="semibold")
    fig.savefig(path, dpi=160, metadata={"Software": "Stage 2 Step 09"})
    plt.close(fig)


def model_style(model):
    return {"color": COLORS[model], "linestyle": {"no_change": "--", "delta_mean": ":"}.get(model, "-"),
            "linewidth": 1.8, "path_effects": OUTLINE if model == "XPert_base" else None}


def draw_pr(metrics, curves, path):
    fig, axes = plt.subplots(2, 2, figsize=(14, 8.9), sharex="col",
                             gridspec_kw={"height_ratios": [3.5, 1]})
    peak_precision = float(curves.precision.max()) if not curves.empty else 1.
    precision_limit = min(1.04, max(.05, .05 * np.ceil(1.05 * peak_precision / .05)))
    legend = []
    for item in metrics.to_dict("records"):
        model = item["model"]
        label = display_name(model)
        if item["status"] == "incomplete":
            label += " (incomplete)"
        elif pd.isna(item["auprc_full"]) and pd.isna(item["auprc_g"]):
            label += " (undefined)"
        legend.append(Line2D([], [], label=label, **model_style(model)))

    for index, (metric, title) in enumerate((("full", "AUPRC@full"), ("g", "AUPRC@G"))):
        panel, support = axes[:, index]
        subset = curves[curves.metric.eq(metric)]
        for model in metrics.model:
            rows = subset[subset.model.eq(model)]
            if rows.empty:
                continue
            style = model_style(model)
            if len(rows) == 1:
                style = {**style, "marker": ".", "markersize": 5}
            panel.step(rows.recall, rows.precision, where="pre", **style)
            support.step(rows.recall, rows.contributing_groups, where="pre", **style)
        panel.set_title(title, fontsize=14, weight="semibold")
        panel.set_ylim(0, precision_limit)
        panel.set_ylabel("Macro precision")
        undefined = metrics.loc[metrics.status.eq("complete") & metrics[f"auprc_{metric}"].isna(), "model"]
        if not undefined.empty:
            names = ", ".join(undefined.map(display_name))
            panel.text(.98, .98, f"Undefined for {names}", ha="right", va="top",
                       transform=panel.transAxes, fontsize=8, wrap=True)
        support.set_ylabel("Groups at recall")
        support.set_xlabel("Recall")
        support.set_xlim(0, 1)
        support.set_ylim(bottom=0)
        for axis in (panel, support):
            axis.grid(alpha=.22)
            axis.set_axisbelow(True)
            axis.spines[["top", "right"]].set_visible(False)

    fig.suptitle(TITLE, fontsize=17, weight="semibold", y=.98)
    fig.legend(handles=legend, loc="lower center", ncol=6, bbox_to_anchor=(.5, .038), frameon=False)
    fig.text(.5, .015,
             "Precision at each recall averages only the eligible groups that reach that recall.",
             ha="center", fontsize=9)
    fig.subplots_adjust(left=.08, right=.98, top=.91, bottom=.18, wspace=.24, hspace=.16)
    fig.savefig(path, dpi=160, metadata={"Software": "Stage 2 Step 09"})
    plt.close(fig)


def render(scores, destination, names):
    metrics = pd.read_csv(scores / "metrics.csv")
    pr = pd.read_csv(scores / "pr.csv")
    draw_metrics(metrics, destination / names[0])
    draw_pr(metrics, pr, destination / names[1])
