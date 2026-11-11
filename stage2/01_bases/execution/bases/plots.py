"""Standalone tenv retrieval panels with the exact plotted numeric source."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from .metrics import METRICS, PR_FIELDS


def render_retrieval(runs, output, *, source='tenv', diagnostic=False):
    """Write a three-direction PR panel and replayable CSV/JSON.

    Each run supplies label, equal-group macro summary and fixed-group pooled
    curves from score_groups. Pooled curve geometry is never labelled as its AP.
    No Tahoe significance or PR publication is supported.
    """
    if source != 'tenv':
        raise ValueError('significance retrieval panels are tenv-only')
    panel_title = ('Diagnostic · measured delta vs tenv truth (not model performance)' if diagnostic
             else 'Tenv · fixed-group pooled precision–recall')
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    output = Path(output)
    if output.suffix != '.png' or not runs:
        raise ValueError('a PNG output and at least one named run are required')
    labels = [str(run['label']) for run in runs]
    if len(labels) != len(set(labels)) or any(not x for x in labels):
        raise ValueError('run labels must be unique and nonempty')
    output.parent.mkdir(parents=True,exist_ok=True)
    rows, summaries = [], []
    fig, axes = plt.subplots(1,3,figsize=(13.5,4.8),sharex=True,sharey=True,layout='constrained')
    titles = ('Combined', 'Upregulated', 'Downregulated')
    colors = plt.get_cmap('tab10').colors
    try:
        for index, run in enumerate(runs):
            label = str(run['label'])
            summaries.append({'label':label, **{name:run['summary'].get(name) for name in METRICS}})
            for name, ax in zip(METRICS,axes):
                curve = run['curves'][name].copy()
                if not set(PR_FIELDS) <= set(curve):
                    raise ValueError('pooled curve schema differs')
                value = run['summary'].get(name)
                if value is not None and (not np.isfinite(value) or not 0 <= value <= 1):
                    raise ValueError('invalid macro AP')
                if len(curve):
                    if (not np.isfinite(curve[['recall','precision']]).all().all()
                            or not curve[['recall','precision']].ge(0).all().all()
                            or not curve[['recall','precision']].le(1).all().all()
                            or not curve.recall.is_monotonic_increasing
                            or curve.eligible_groups.nunique() != 1
                            or not curve.contributing_groups.eq(curve.eligible_groups).all()):
                        raise ValueError('expected fixed-population monotonic PR source')
                    ax.step(curve.recall,curve.precision,where='post',color=colors[index%len(colors)],
                            linewidth=1.8,label=f'{label} · macro AP {value:.3f}' if value is not None else label)
                    rows.append(curve.assign(label=label,metric=name))
                elif value is not None:
                    raise ValueError('eligible AP has no PR curve')
        for ax,title in zip(axes,titles):
            ax.set(title=title,xlabel='Recall',xlim=(0,1.01),ylim=(0,1.04))
            ax.grid(alpha=.18)
            if ax.get_legend_handles_labels()[0]:
                ax.legend(loc='upper center',bbox_to_anchor=(.5,-.23),fontsize=8,frameon=False)
            else:
                ax.text(.5,.5,'No eligible measured groups',ha='center',va='center',transform=ax.transAxes)
        axes[0].set_ylabel('Precision')
        fig.suptitle(panel_title,fontsize=13)
        fig.savefig(output,dpi=160)
    finally:
        plt.close(fig)
    frame = pd.concat(rows,ignore_index=True) if rows else pd.DataFrame(columns=[*PR_FIELDS,'label','metric'])
    frame.to_csv(output.with_suffix('.csv'),index=False)
    report = {'schema':'stage2-indexed-pr-panel-v1','source':source,'diagnostic':bool(diagnostic),'title':panel_title,'summaries':summaries,
              'curve_semantics':'fixed eligible groups pooled at shared magnitude thresholds; displayed AP is exact equal-group macro',
              'files':{output.name:hashlib.sha256(output.read_bytes()).hexdigest(),
                       output.with_suffix('.csv').name:hashlib.sha256(output.with_suffix('.csv').read_bytes()).hexdigest()}}
    output.with_suffix('.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    return report


def replay_retrieval(report, output):
    """Replay the panel solely from its sealed CSV and macro summaries."""
    report = Path(report)
    value = json.loads(report.read_text())
    if value.get('schema') != 'stage2-indexed-pr-panel-v1' or value.get('source') != 'tenv':
        raise ValueError('unknown tenv PR panel source')
    source = report.with_suffix('.csv')
    if hashlib.sha256(source.read_bytes()).hexdigest() != value['files'].get(source.name):
        raise ValueError('PR numeric source seal differs')
    frame = pd.read_csv(source)
    runs = []
    for summary in value['summaries']:
        label = summary['label']
        curves = {name:frame.loc[frame.label.eq(label)&frame.metric.eq(name),list(PR_FIELDS)] for name in METRICS}
        runs.append({'label':label,'summary':summary,'curves':curves})
    return render_retrieval(runs,output,source=value['source'],diagnostic=value.get('diagnostic',False))
