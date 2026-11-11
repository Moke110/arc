"""The common Stage 2 PR curve and AUPRC presentation."""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.hashsalt'] = 'stage2-significance-v1'
import numpy as np
import pandas as pd
from stage2.evaluation.metrics import METRICS


def render(root):
    data = pd.read_csv(root / 'per_run.csv')
    curves = pd.read_csv(root / 'pr.csv')
    if data.role.eq('test').any():
        data = data[data.role.eq('test')]
        curves = curves[curves.role.eq('test')]
    strata = data[['dataset', 'regime', 'role']].drop_duplicates().values.tolist()
    models = data.model.drop_duplicates().tolist()
    colors = {m: plt.get_cmap('tab20')(i % 20) for i, m in enumerate(models)}
    for m in ('no_change', 'delta_mean'):
        colors[m] = '#777777' if m == 'no_change' else '#aaaaaa'
    fig, axes = plt.subplots(len(strata), 3, figsize=(16, 4.5 * len(strata)), squeeze=False)
    for row, (dataset, regime, role) in enumerate(strata):
        source = curves[curves.dataset.eq(dataset) & curves.regime.eq(regime) & curves.role.eq(role)]
        for axis, metric in zip(axes[row], METRICS):
            for model in models:
                points = source[source.model.eq(model) & source.metric.eq(metric)]
                if points.empty:
                    continue
                axis.step(points.recall, points.precision, where='pre', color=colors[model],
                          label=model, linewidth=1.6, marker='.' if len(points) == 1 else None)
                end = points.iloc[-1]
                axis.plot(end.recall, end.precision, '.', color=colors[model], markersize=5)
            axis.set_title(f"PR curve{metric.removeprefix('AUPRC')}\n{dataset} · {regime} · {role}")
            axis.set(xlabel='Recall', ylabel='Precision', xlim=(0, 1), ylim=(0, 1.02))
            axis.grid(alpha=.22)
            zoom = axis.inset_axes([.23, .35, .73, .55])
            for line in axis.lines[::2]:
                zoom.step(line.get_xdata(), line.get_ydata(), where='pre', color=line.get_color(), linewidth=1.3)
            zoom_top = .25 if metric == 'AUPRC_positive' else .08
            zoom.set(xlim=(0, 1), ylim=(0, zoom_top))
            zoom.tick_params(labelsize=7)
            zoom.set_title(f'Precision 0–{zoom_top:g}', fontsize=8)
            zoom.grid(alpha=.18)
            axis.spines[['top', 'right']].set_visible(False)
    legend = {}
    for axis in axes[:, 0]:
        handles, labels = axis.get_legend_handles_labels()
        legend.update(zip(labels, handles))
    fig.legend(legend.values(), legend.keys(), loc='lower center', bbox_to_anchor=(.5, .045), ncol=min(7, len(models)), frameon=False, fontsize=9)
    fig.suptitle('PR curve', fontsize=18)
    fig.text(.5, .015, 'Fixed-group pooled retrieval; direction-matched significance on 18,164 CM genes; endpoints mark attainable recall.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .18, 1, .96))
    for suffix in ('png', 'svg'):
        fig.savefig(root / f'PR_curve.{suffix}', dpi=150, metadata={'Creator': 'Stage 2 evaluation', 'Date': '2026-11-11'} if suffix == 'svg' else {'Software': 'Stage 2 evaluation'})
    plt.close(fig)
    extra = [m for m in ('PCC', 'Spearman', 'RMSE', 'R2') if m in data]
    metrics = [*METRICS, *extra]
    metric_rows = (len(metrics) + 2) // 3
    height = max(4.2, .3 * len(models) + 1.6)
    fig, axes = plt.subplots(len(strata) * metric_rows, 3, figsize=(16, height * len(strata) * metric_rows), squeeze=False)
    for row, (dataset, regime, role) in enumerate(strata):
        subset = data[data.dataset.eq(dataset) & data.regime.eq(regime) & data.role.eq(role)]
        local_models = subset.model.drop_duplicates().tolist()
        for i in range(metric_rows * 3):
            axis = axes[row * metric_rows + i // 3, i % 3]
            if i >= len(metrics):
                axis.axis('off')
                continue
            metric = metrics[i]
            means = subset.groupby('model', sort=False)[metric].mean().reindex(local_models)
            axis.barh(local_models, means, color=[colors[m] for m in local_models])
            for index, value in enumerate(means):
                if np.isfinite(value):
                    axis.text(.98, index, f'{value:.4f}', va='center', ha='right', transform=axis.get_yaxis_transform(), fontsize=8,
                              bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .75})
            axis.invert_yaxis()
            axis.set_title(f'{metric}\n{dataset} · {regime} · {role}')
            axis.grid(axis='x', alpha=.22)
            axis.set_axisbelow(True)
            axis.spines[['top', 'right']].set_visible(False)
    fig.suptitle('AUPRC', fontsize=18)
    fig.text(.5, .015, 'Exact group AP, equally averaged within run; run means shown. Undefined groups are excluded and counted in the tables.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .04, 1, .95))
    for suffix in ('png', 'svg'):
        fig.savefig(root / f'AUPRC.{suffix}', dpi=150, metadata={'Creator': 'Stage 2 evaluation', 'Date': '2026-11-11'} if suffix == 'svg' else {'Software': 'Stage 2 evaluation'})
    plt.close(fig)
