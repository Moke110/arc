"""Reproduce the six report panels from sibling CSV files only."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

METHODS=['no-dose','scaling','input-numeric','latent-numeric','input-embedding','latent-embedding','film','selected']
COLORS=dict(zip(METHODS,['#78838d','#e6a23c','#3985ad','#64a6bd','#865b9a','#b179ad','#3c9d77','#1b3248']))


def rows(path):
    with path.open(newline='') as stream:return list(csv.DictReader(stream))


def save(fig,out,name):
    fig.savefig(out/(name+'.png'),dpi=180,bbox_inches='tight');plt.close(fig)


def render(data,out):
    out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    primary=rows(data/'primary_values.csv')
    paired=rows(data/'primary_paired.csv')
    fig,axes=plt.subplots(1,2,figsize=(12,4),gridspec_kw={'width_ratios':[2,1]})
    axes[0].bar([r['method'] for r in primary],[float(r['pcc']) for r in primary],color=[COLORS[r['method']] for r in primary])
    axes[0].tick_params(axis='x',rotation=40);axes[0].set_ylabel('Parent-equal Δ PCC');axes[0].set_title('Cold-drug OOF: fixed configurations')
    point=next(r for r in paired if r['metric']=='pcc');estimate=float(point['estimate']);low=float(point['low']);high=float(point['high'])
    axes[1].errorbar(estimate,0,xerr=[[estimate-low],[high-estimate]],fmt='o',color=COLORS['selected'],capsize=5)
    axes[1].axvline(0,color='#9aa1a8',lw=1);axes[1].set_yticks([0],['Selected − no-dose']);axes[1].set_xlabel('Δ PCC difference');axes[1].set_title('95% conditional paired interval')
    fig.tight_layout();save(fig,out,'primary')
    for table,name,title in [('primary_dose_values.csv','primary_dose','Cold-drug per-dose OOF'),('held_dose_values.csv','held_dose','Held-dose: fixed comparison, no reselection')]:
        values=rows(data/table);fig,ax=plt.subplots(figsize=(10,4))
        for method in METHODS:
            part=sorted([r for r in values if r['method']==method],key=lambda r:float(r['dose_um']))
            if part:
                label=method;highlight=False
                if name=='held_dose':
                    selected=[str(f) for f in range(3) if part[0].get(f'primary_fold_{f}_selected_method')==method]
                    tags=(['primary fold '+','.join(selected)] if selected else [])+(['final'] if part[0].get('final_selected_method')==method else [])
                    if tags:label+=' ('+'; '.join(tags)+')';highlight=True
                ax.plot([float(r['dose_um']) for r in part],[float(r['pcc']) for r in part],'-o',label=label,color=COLORS[method],lw=2.2 if highlight else 1.2)
                if name=='held_dose':
                    defined=[r for r in part if r.get('pcc_low') not in [None,''] and r.get('pcc_high') not in [None,'']]
                    ax.vlines([float(r['dose_um']) for r in defined],[float(r['pcc_low']) for r in defined],[float(r['pcc_high']) for r in defined],color=COLORS[method],alpha=.55,lw=1.2)
        ax.set_xscale('log');ax.set_xticks([.2,1,10],['0.2\nlower extrapolation' if name=='held_dose' else '0.2','1\ninterpolation' if name=='held_dose' else '1','10\nupper extrapolation' if name=='held_dose' else '10'])
        ax.set_xlabel('Dose (µM)');ax.set_ylabel('Parent-equal Δ PCC');ax.set_title(title+('\nDescriptive 95% parent-block intervals; three measured doses' if name=='held_dose' else ''));ax.legend(fontsize=8,ncol=2);fig.tight_layout();save(fig,out,name)
    validation=rows(data/'final_validation_values.csv');fig,ax=plt.subplots(figsize=(11,4))
    for method in METHODS:
        part=[r for r in validation if r['scenario']=='final' and r['method']==method]
        if part:ax.scatter([method]*len(part),[float(r['pcc']) for r in part],color=COLORS[method],s=32)
    ax.tick_params(axis='x',rotation=35);ax.set_ylabel('Validation parent-equal Δ PCC');ax.set_title('Deployment selection: three seed scores, no independent test');fig.tight_layout();save(fig,out,'final_validation')
    curves=rows(data/'learning_curves_all.csv');fig,axes=plt.subplots(1,3,figsize=(15,4),sharey=True)
    for ax,scenario in zip(axes,['cold-drug','held-dose','final']):
        groups=defaultdict(list)
        for row in curves:
            if row['scenario']==scenario:groups[(row['method'],row['fold'],row['seed'])].append(row)
        seen=set()
        for (method,fold,seed),part in groups.items():
            part.sort(key=lambda r:int(r['epoch']))
            ax.plot([int(r['epoch']) for r in part],[float(r['validation_mse']) for r in part],color=COLORS[method],alpha=.35,lw=.7,label=method if method not in seen else None);seen.add(method)
        ax.set_title(scenario);ax.set_xlabel('Epoch');ax.set_yscale('log')
    axes[0].set_ylabel('Complete validation expression MSE');axes[-1].legend(fontsize=7);fig.tight_layout();save(fig,out,'learning_curves')
    resources=rows(data/'resources_all.csv');fig,ax=plt.subplots(figsize=(10,4))
    for method in METHODS:
        part=[r for r in resources if r['method']==method]
        if part:ax.scatter([int(r['parameter_count'])/1e6 for r in part],[float(r['elapsed_seconds']) for r in part],label=method,color=COLORS[method],alpha=.65,s=24)
    ax.set_xlabel('Parameters (millions)');ax.set_ylabel('Fit elapsed seconds');ax.set_title('All 147 fixed fits');ax.legend(fontsize=8,ncol=2);fig.tight_layout();save(fig,out,'resources')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--data',type=Path,default=Path(__file__).resolve().parents[1]/'data');parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]);args=parser.parse_args();render(args.data,args.output)
