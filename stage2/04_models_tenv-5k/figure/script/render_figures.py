"""Render complete recorded plot values; no model, input or metric access."""
import argparse
import importlib.util
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

# Pure public palette module reads no payload (including Bases __init__).
PALETTE=Path(__file__).resolve().parents[3]/'01_bases/output/api/bases/palette.py'
spec=importlib.util.spec_from_file_location('stage_palette',PALETTE)
palette=importlib.util.module_from_spec(spec);spec.loader.exec_module(palette)
HUES=dict(PRnet='blue',TranSiGen='orange',MLP_Morgan='green',Mean='gray',Mean_cell='purple',Mean_drug='magenta')
ORDER=list(HUES)


def metric_annotation(row):
    """All status/count text comes from the tracked plot-value row."""
    if not bool(row.applicable):return 'N/A: '+str(row.reason)
    return ('n='+str(int(row.n_parents))+'; undefined='+str(int(row.n_undefined))+
            '; truth constant='+str(int(row.n_constant_truth))+
            '; pred constant='+str(int(row.n_constant_prediction)))


def render(data,out):
    data=Path(data);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':140})
    values=pd.read_csv(data/'metrics_values.csv')
    fig,axes=plt.subplots(2,4,figsize=(17,7),constrained_layout=True)
    for i,scenario in enumerate(('cold-drug','warm-start')):
        for j,metric in enumerate(('pcc','spearman','rmse','r2')):
            ax=axes[i,j];rows=values[(values.scenario==scenario)&(values.metric==metric)].set_index('method')
            methods=[m for m in ORDER if m in rows.index]
            ax.bar(range(len(methods)),[rows.loc[m,'value'] for m in methods],color=[palette.color(HUES[m],6) for m in methods],edgecolor=[palette.color(HUES[m],2) for m in methods])
            labels=[]
            for k,m in enumerate(methods):
                row=rows.loc[m]
                labels.append(m+'\n'+('N/A' if not bool(row.applicable) else 'n='+str(int(row.n_parents))+' u='+str(int(row.n_undefined))+' t='+str(int(row.n_constant_truth))+' p='+str(int(row.n_constant_prediction))))
                if not bool(row.applicable):ax.text(k,0,'N/A\n'+str(row.reason),rotation=90,ha='center',va='bottom',fontsize=6)
            ax.set_xticks(range(len(methods)),labels,rotation=50,ha='right',fontsize=7);ax.set_title(scenario+' · Δ '+metric.upper())
            ax.axhline(float(rows.zero.iloc[0]),color=palette.color('gray',4),lw=.6)
            ax.grid(axis='y',alpha=.15)
    fig.suptitle('Equal-weight parent metrics; n=defined parents, u=undefined; t=constant truth parents, p=constant prediction seed-parent vectors. IDs/reasons in CSV')
    fig.savefig(out/'metrics.png');plt.close(fig)
    pairs=pd.read_csv(data/'paired_differences_values.csv')
    fig,axes=plt.subplots(2,4,figsize=(19,10),constrained_layout=True)
    for i,scenario in enumerate(('cold-drug','warm-start')):
        for j,metric in enumerate(('pcc','spearman','rmse','r2')):
            ax=axes[i,j];rows=pairs[(pairs.scenario==scenario)&(pairs.metric==metric)]
            for k,row in enumerate(rows.itertuples()):
                color=palette.color(HUES[row.candidate],2)
                ax.plot([row.low,row.high],[k,k],color=color,lw=2);ax.plot(row.estimate,k,'o',color=color,ms=4)
            ax.set_yticks(range(len(rows)),[r.candidate+' − '+r.reference+' [pred const '+str(int(r.candidate_constant_prediction))+'/'+str(int(r.reference_constant_prediction))+']' for r in rows.itertuples()],fontsize=7)
            if len(rows):ax.axvline(float(rows.zero.iloc[0]),color=palette.color('gray',4),lw=.7)
            excluded=int(rows.n_excluded_constant.iloc[0]) if len(rows) else 0
            unavailable=str(rows.unavailable_reference.iloc[0]) if len(rows) and pd.notna(rows.unavailable_reference.iloc[0]) else ''
            ax.set_title(scenario+' · Δ '+metric.upper()+'; excluded='+str(excluded),fontsize=9)
            if unavailable:ax.set_xlabel('N/A: '+unavailable,fontsize=7)
            ax.grid(axis='x',alpha=.15)
    fig.suptitle('Paired drug bootstrap: candidate minus reference, 95% CI (errors: negative is better)')
    fig.savefig(out/'paired_differences.png');plt.close(fig)
    curves=pd.read_csv(data/'learning_curve_values.csv')
    fig,axes=plt.subplots(2,3,figsize=(17,9),constrained_layout=True)
    for i,scenario in enumerate(('cold-drug','warm-start')):
        for j,method in enumerate(('PRnet','TranSiGen','MLP_Morgan')):
            ax=axes[i,j];validation_ax=ax.twinx()
            validation_ax.spines['right'].set_visible(True)
            rows=curves[(curves.scenario==scenario)&(curves.method==method)]
            for index,(_,run) in enumerate(rows.groupby(['fold','seed'])):
                # Every recorded epoch is plotted; shades identify the nine fold/seed runs.
                validation_ax.plot(run.epoch,run.validation_mse,color=palette.color('green',index%7),alpha=.65,lw=.8)
                ax.plot(run.epoch,run.train_loss,color=palette.color('blue',index%7),alpha=.35,lw=.6,ls='--')
            ax.set_title(scenario+' · '+method);ax.set_xlabel('epoch')
            ax.set_ylabel('native train loss (dashed)',color=palette.color('blue',2))
            validation_ax.set_ylabel('validation MSE (solid)',color=palette.color('green',2))
            ax.tick_params(axis='y',labelcolor=palette.color('blue',2))
            validation_ax.tick_params(axis='y',labelcolor=palette.color('green',2))
            ax.grid(alpha=.15)
    fig.suptitle('Independent y-scales: native training loss (blue, left) and validation MSE (green, right). Native losses are not comparable across models')
    fig.savefig(out/'learning_curve.png');plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--figure-data',type=Path,required=True);p.add_argument('--out-dir',type=Path,required=True)
    a=p.parse_args();render(a.figure_data,a.out_dir)
