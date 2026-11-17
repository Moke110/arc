"""Reproduce Stage2-style metric and learning-curve figures from CSV files only."""
import argparse
import csv
import math
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


METRICS=['pcc','spearman','rmse','r2']
SCENARIOS=['cold-drug','held-dose','final']
LABELS={'no-dose':'No dose','scaling':'Scaling','input-numeric':'Input numeric','latent-numeric':'Latent numeric','input-embedding':'Input embedding','latent-embedding':'Latent embedding','film':'FiLM','selected':'Selected procedure'}


def number(value):
    return float(value) if value not in [None,''] else float('nan')


def metric_figures(data,out,fills):
    values=rows(data/'metrics_values.csv')
    fig,axes=plt.subplots(2,4,figsize=(19,8.5),constrained_layout=True)
    for i,scenario in enumerate(SCENARIOS[:2]):
        for j,metric in enumerate(METRICS):
            ax=axes[i,j];lookup={r['method']:r for r in values if r['scenario']==scenario and r['metric']==metric}
            methods=[m for m in METHODS if m in lookup]
            ax.bar(range(len(methods)),[number(lookup[m]['value']) for m in methods],color=[fills[m] for m in methods],edgecolor=[COLORS[m] for m in methods])
            labels=[]
            for k,method in enumerate(methods):
                r=lookup[method];labels.append(LABELS[method]+'\n'+f"n={r['n_defined']} u={r['n_undefined']} t={r['n_constant_truth']} p={r['n_constant_prediction']}")
                if not math.isfinite(number(r['value'])):ax.text(k,0,'Undefined',rotation=90,ha='center',va='bottom',fontsize=6)
            ax.set_xticks(range(len(methods)),labels,rotation=50,ha='right',fontsize=7)
            ax.set_title(scenario+' · Δ '+metric.upper());ax.axhline(0,color='#949494',lw=.6);ax.grid(axis='y',alpha=.15)
    fig.suptitle('Equal-weight parent metrics (three seed metrics averaged); n=defined parents, u=undefined; t=constant truth parent–dose vectors, p=constant prediction seed–parent–dose vectors')
    save(fig,out,'metrics')
    primary=rows(data/'primary_paired.csv');descriptive=rows(data/'primary_descriptive.csv')
    fig,axes=plt.subplots(2,4,figsize=(19,8),constrained_layout=True)
    for i,(part,title) in enumerate([(primary,'Primary selected procedure − no dose'),([r for r in descriptive if r.get('reference')=='no-dose'],'Descriptive fixed candidate − no dose')]):
        for j,metric in enumerate(METRICS):
            ax=axes[i,j];points=[r for r in part if r['metric']==metric]
            for k,r in enumerate(points):
                candidate=r.get('candidate','selected');color=COLORS[candidate]
                ax.plot([number(r['low']),number(r['high'])],[k,k],color=color,lw=2);ax.plot(number(r['estimate']),k,'o',color=color,ms=4)
            ax.set_yticks(range(len(points)),[LABELS.get(r.get('candidate','selected'),r.get('candidate','selected'))+' − No dose' for r in points],fontsize=7)
            ax.axvline(0,color='#949494',lw=.7);ax.set_title(title+'\nΔ '+metric.upper(),fontsize=9);ax.grid(axis='x',alpha=.15)
            if metric=='rmse':ax.set_xlabel('Negative is better',fontsize=8)
    fig.suptitle('Paired parent bootstrap: candidate minus reference, 95% intervals; selections fixed (method-selection uncertainty omitted)')
    save(fig,out,'paired_differences')


def learning_figures(data,out,styles):
    curves=rows(data/'learning_curves_all.csv')
    groups=defaultdict(list)
    for row in curves:groups[(row['scenario'],row['method'],row['fold'],row['seed'])].append(row)
    for part in groups.values():part.sort(key=lambda r:int(r['epoch']))
    def panel(ax,scenario,method):
        validation_ax=ax.twinx();validation_ax.spines['right'].set_visible(True)
        runs=sorted([(key,part) for key,part in groups.items() if key[:2]==(scenario,method)],key=lambda item:(int(item[0][2]),int(item[0][3])))
        for index,(key,part) in enumerate(runs):
            style=styles[index];epochs=[int(r['epoch']) for r in part]
            ax.plot(epochs,[number(r['train_loss']) for r in part],color=style['train'],alpha=.45,lw=.65,ls='--')
            validation_ax.plot(epochs,[number(r['validation_mse']) for r in part],color=style['validation'],alpha=.7,lw=.85)
        ax.set_title(scenario+' · '+LABELS[method],fontsize=9);ax.set_xlabel('epoch',fontsize=8)
        ax.set_ylabel('native train loss (dashed)',color=styles[0]['train'],fontsize=7)
        validation_ax.set_ylabel('validation MSE (solid)',color=styles[0]['validation'],fontsize=7)
        ax.tick_params(axis='y',labelcolor=styles[0]['train'],labelsize=7);validation_ax.tick_params(axis='y',labelcolor=styles[0]['validation'],labelsize=7)
        ax.tick_params(axis='x',labelsize=7);ax.grid(alpha=.15)
    caption='Independent y-scales: native summed VAE train loss (blue, left) and validation expression MSE (green, right); every epoch and fold/seed run shown'
    fig,axes=plt.subplots(3,7,figsize=(33,12),constrained_layout=True)
    for i,scenario in enumerate(SCENARIOS):
        for j,method in enumerate(METHODS[:-1]):panel(axes[i,j],scenario,method)
    fig.suptitle(caption);save(fig,out,'learning_curve')
    for scenario in SCENARIOS:
        fig,axes=plt.subplots(2,4,figsize=(19,9),constrained_layout=True)
        for ax,method in zip(axes.flat,METHODS[:-1]):panel(ax,scenario,method)
        axes.flat[-1].axis('off');axes.flat[-1].text(.05,.8,'Blue dashed: native train loss\nGreen solid: validation MSE\nSeparate y-scales; no curve averaging\nShades identify fold/seed runs\nFinal = internal validation only',va='top',fontsize=11)
        fig.suptitle(caption,fontsize=11);save(fig,out,'learning_curve_'+scenario)


def render(data,out):
    out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    style=rows(data/'plot_style.csv');fills={r['name']:r['fill'] for r in style if r['role']=='method'}
    COLORS.update({r['name']:r['line'] for r in style if r['role']=='method'})
    runstyles={int(r['name']):r for r in style if r['role']=='run'}
    metric_figures(data,out,fills);learning_figures(data,out,runstyles)
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
