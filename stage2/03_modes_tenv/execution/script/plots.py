"""All three figures regenerate exclusively from checked evaluator artifacts."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from artifacts import rows
from metrics import METRICS


def render(output):
    output=Path(output)
    summary=rows(output/'summary.csv'); runs=rows(output/'per_run.csv')
    models=[r['model'] for r in summary]
    labels=[f"{r['model']} ({r['completed_runs']}/9)" + (' [deterministic copies]' if r['deterministic_reference']=='True' else '') for r in summary]
    figure,axes=plt.subplots(2,4,figsize=(19,10),constrained_layout=True)
    for metric,axis in zip(METRICS,axes.flat):
        for x,row in enumerate(summary):
            mean=float(row[f'{metric}_mean']); sd=float(row[f'{metric}_sample_sd'])
            if np.isfinite(mean):
                axis.bar(x,mean,color='lightsteelblue',yerr=sd if np.isfinite(sd) else None,capsize=3)
            else: axis.text(x,.05,'undefined',rotation=90,ha='center',transform=axis.get_xaxis_transform())
            dots=[float(r[metric]) for r in runs if r['model']==row['model']]
            axis.scatter(x+np.linspace(-.13,.13,len(dots)),dots,s=13,color='black',zorder=3)
        axis.set_title(metric); axis.set_xticks(range(len(models)),labels,rotation=35,ha='right',fontsize=8)
    figure.suptitle('Nine fold × seed run samples: mean ± sample SD; incomplete entries cannot rank\nRun variation does not represent independent biological replication')
    figure.savefig(output/'metrics.png',dpi=150); plt.close(figure)
    bands=rows(output/'pr_bands.csv')
    figure,axes=plt.subplots(1,len(models),figsize=(max(6,5*len(models)),5),squeeze=False,constrained_layout=True)
    for model,label,axis in zip(models,labels,axes.flat):
        for direction,color in [('positive','red'),('negative','blue')]:
            selected=[r for r in bands if r['model']==model and r['direction']==direction]
            recall=np.array([float(r['recall']) for r in selected]); mean=np.array([float(r['mean']) for r in selected])
            low=np.array([float(r['min']) for r in selected]); high=np.array([float(r['max']) for r in selected])
            axis.plot(recall,mean,color=color,label=direction); axis.fill_between(recall,low,high,color=color,alpha=.15)
        axis.set(xlabel='Recall',ylabel='Precision',title=label,xlim=(0,1),ylim=(0,1.01)); axis.legend()
    figure.suptitle('Group macro PR within run → run mean and min–max band; tied thresholds, no synthetic perfect endpoint')
    figure.savefig(output/'pr_panels.png',dpi=150); plt.close(figure)
    source=np.load(output/'pr_sources.npz',allow_pickle=False)
    ids=source['identifiers'].tolist()
    figure,axis=plt.subplots(figsize=(8,6),constrained_layout=True)
    for selected in rows(output/'best_runs.csv'):
        if not selected['fold']:
            axis.plot([],[],label=f"{selected['model']}: undefined selection"); continue
        key=f"{selected['model']}:{selected['fold']}:{selected['seed']}"
        curves=source['curves'][ids.index(key)]
        # Both directions must be defined for best-run eligibility.
        curve=(curves[0]+curves[1])/2
        axis.plot(source['recall'],curve,label=f"{selected['model']} (fold {selected['fold']}, seed {selected['seed']})")
    axis.set(xlabel='Recall',ylabel='Mean directional precision',xlim=(0,1),ylim=(0,1.01),
             title='Post-hoc best test run per entry\nDifferent entries may select different folds; nine-run summaries are primary')
    axis.legend(fontsize=8); figure.savefig(output/'best_run.png',dpi=150); plt.close(figure)
