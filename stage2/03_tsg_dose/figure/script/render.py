"""Render exactly three round1 figures from plotting CSVs alone."""
import argparse
from collections import defaultdict
import csv
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

METHODS=['no-dose','scaling','input-numeric','latent-numeric','input-embedding','latent-embedding','film','selected']
SCENARIOS=['cold-drug','held-dose']
METRICS=['pcc','spearman','rmse','r2']
LABELS={'no-dose':'No dose','scaling':'Scaling','input-numeric':'Input numeric','latent-numeric':'Latent numeric','input-embedding':'Input embedding','latent-embedding':'Latent embedding','film':'FiLM','selected':'Selected procedure'}
TITLES={'pcc':'Δ PCC ↑','spearman':'Δ Spearman ↑','rmse':'Δ RMSE ↓','r2':'Δ R² ↑'}

def rows(path):
    with path.open(newline='') as stream:return list(csv.DictReader(stream))

def number(value):return float(value) if value not in ['',None] else float('nan')

def save(fig,out,name):
    fig.savefig(out/name,dpi=180,bbox_inches='tight',metadata={'Software':'ARC round1 CSV renderer'});plt.close(fig)

def metrics(data,out,style):
    points=rows(data/'round1_runs.csv');summaries=rows(data/'round1_summary.csv')
    fig,axes=plt.subplots(2,4,figsize=(21,10))
    for i,scenario in enumerate(SCENARIOS):
        methods=METHODS if i==0 else METHODS[:-1]
        for j,metric in enumerate(METRICS):
            ax=axes[i,j];labels=[]
            for x,method in enumerate(methods):
                summary=next(r for r in summaries if (r['scenario'],r['method'],r['metric'])==(scenario,method,metric))
                color=style[method]['line'];mean=number(summary['mean']);sd=number(summary['sd'])
                if math.isfinite(mean):
                    ax.bar(x,mean,color=style[method]['fill'],edgecolor=color,lw=1.3,width=.72)
                    if math.isfinite(sd):ax.errorbar(x,mean,yerr=sd,color=color,capsize=4,lw=1.5)
                else:ax.text(x,0,'undefined',ha='center',va='bottom',rotation=90,fontsize=7)
                part=[r for r in points if (r['scenario'],r['method'],r['metric'])==(scenario,method,metric) and math.isfinite(number(r['value']))]
                ax.scatter([float(r['x']) for r in part],[float(r['value']) for r in part],color=color,s=19,zorder=3,edgecolors='white',linewidths=.3)
                labels.append(LABELS[method]+f"\nn={summary['n_defined']}/{summary['n_total']}, u={summary['n_undefined']}")
            ax.set_xticks(range(len(methods)),labels,rotation=48,ha='right',fontsize=8)
            ax.set_title(scenario+' · '+TITLES[metric],fontsize=11);ax.axhline(0,color='.65',lw=.6);ax.grid(axis='y',alpha=.15)
    fig.suptitle('Round1 · equal-weight fold × seed run mean ± sample SD (ddof=1)',fontsize=16)
    fig.text(.5,.018,'One dot = one test run; n = finite / total runs, u = undefined. Finite dose metrics → parent mean → run mean; no zero substitution.\nSD describes non-independent run spread, not a CI or biological replication. Held-dose folds omit different doses. Historical parent-equal OOF inference is separate.',ha='center',fontsize=10)
    fig.subplots_adjust(left=.045,right=.99,top=.92,bottom=.23,hspace=.75,wspace=.27)
    save(fig,out,'round1_metrics.png')

def learning(data,out,style,runstyles):
    groups=defaultdict(list)
    for row in rows(data/'round1_curves.csv'):groups[(row['scenario'],row['method'],int(row['fold']),int(row['seed']))].append(row)
    for part in groups.values():part.sort(key=lambda r:int(r['epoch']))
    fig,axes=plt.subplots(2,8,figsize=(40,12))
    for i,scenario in enumerate(SCENARIOS):
        for j,method in enumerate(METHODS[:-1]):
            ax=axes[i,j];right=ax.twinx();right.spines['right'].set_visible(True)
            for index,(fold,seed) in enumerate((f,s) for f in range(3) for s in [17,42,2026]):
                part=groups[(scenario,method,fold,seed)];epochs=[int(r['epoch']) for r in part];colors=runstyles[index]
                ax.plot(epochs,[float(r['train_loss']) for r in part],ls='--',color=colors['train'],alpha=.8,lw=.8)
                right.plot(epochs,[float(r['validation_mse']) for r in part],color=colors['validation'],alpha=.85,lw=.85)
            ax.set_title(scenario+' · '+LABELS[method],fontsize=12)
            ax.set_xlabel('Epoch');ax.set_ylabel('Native train loss · dashed',color=runstyles[0]['train'],fontsize=9)
            right.set_ylabel('Validation MSE · solid',color=runstyles[0]['validation'],fontsize=9)
            ax.tick_params(labelsize=8);right.tick_params(labelsize=8);ax.grid(alpha=.12)
        ax=axes[i,7]
        for method in METHODS[:-1]:
            seen=False
            for (s,m,fold,seed),part in groups.items():
                if (s,m)!=(scenario,method):continue
                ax.plot([int(r['epoch']) for r in part],[float(r['validation_mse']) for r in part],color=style[method]['line'],alpha=.5,lw=.9,label=LABELS[method] if not seen else None);seen=True
        ax.set_yscale('log');ax.set_title(scenario+' · all methods',fontsize=12);ax.set_xlabel('Epoch');ax.set_ylabel('Validation MSE · log scale');ax.legend(fontsize=9,loc='upper right');ax.grid(alpha=.15)
    handles=[]
    for index,(fold,seed) in enumerate((f,s) for f in range(3) for s in [17,42,2026]):
        handles.append(Line2D([0],[0],color=runstyles[index]['train'],ls='--',label=f'Fold {fold} / seed {seed} · train'))
        handles.append(Line2D([0],[0],color=runstyles[index]['validation'],label=f'Fold {fold} / seed {seed} · validation'))
    fig.legend(handles=handles,loc='lower center',ncol=9,fontsize=10,bbox_to_anchor=(.5,.018))
    fig.suptitle('Round1 · every epoch and fold/seed run · native loss (blue dashed, left) and validation MSE (green solid, right)',fontsize=18)
    fig.text(.5,.105,'Native panels use independent y-scales; shades identify the same fold/seed in both scenarios. No curve averaging or final/deployment panels.',ha='center',fontsize=12)
    fig.subplots_adjust(left=.025,right=.98,top=.92,bottom=.20,hspace=.36,wspace=.67)
    save(fig,out,'round1_learning_curves.png')

def doses(data,out,style):
    values=rows(data/'round1_dose.csv');fig,axes=plt.subplots(1,4,figsize=(21,5.8))
    for ax,metric in zip(axes,METRICS):
        for method in METHODS:
            part=sorted([r for r in values if (r['method'],r['metric'])==(method,metric)],key=lambda r:float(r['dose_um']))
            ax.plot([float(r['dose_um']) for r in part],[number(r['value']) for r in part],'-o',color=style[method]['line'],lw=2 if method=='selected' else 1.3,label=LABELS[method],ms=4)
        ax.set_xscale('log');ax.set_xticks([.2,1,10],['0.2','1','10']);ax.set_xlabel('Dose (µM)');ax.set_title(TITLES[metric]);ax.grid(alpha=.15)
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=8,fontsize=9,bbox_to_anchor=(.5,.02))
    fig.suptitle('Round1 · cold-drug test performance by measured dose',fontsize=16)
    fig.text(.5,.115,'Saved parent-equal, equal-seed OOF dose summaries; selected procedure uses locked cold-drug fold choices. No held-dose selection.',ha='center',fontsize=10)
    fig.subplots_adjust(top=.86,bottom=.27,wspace=.25)
    save(fig,out,'round1_metrics_cold_drug_dose.png')

def render(data,out):
    out.mkdir(parents=True,exist_ok=True);plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    style_rows=rows(data/'plot_style.csv');style={r['name']:r for r in style_rows if r['role']=='method'};runstyles={int(r['name']):r for r in style_rows if r['role']=='run'}
    metrics(data,out,style);learning(data,out,style,runstyles);doses(data,out,style)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--data',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();render(args.data,args.output)
