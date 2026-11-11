"""ADR 0060 native-axis scoring and replayable 24-run publication.

score_run consumes selected-checkpoint deltas, native observed deltas and the
canonical Labels.truth callback. well_loss must be the full test well-weighted
XPert objective, evaluated with the selected checkpoint (never group MSE).
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(next(p for p in Path(__file__).resolve().parents if (p / '.git').exists())))
from stage2.evaluation.metrics import METRICS as RETRIEVAL, PooledPR, retrieval
from round_metrics import continuous

MODELS = tuple(f'XPert_{kind}{size}' for kind in ('G', 'MLP') for size in (512, 1024, 2048))
METRICS = ('test_loss', 'pcc', 'spearman', 'rmse', 'r2', *RETRIEVAL)
ENTRIES = {'checkpoints', 'logs', 'figure_data', 'r1_learning_curves.png', 'r1_metrics.png', 'r1_PR_curves.png'}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def score_run(cases, prediction, observed_delta, covered_columns, truth, well_loss):
    """Bounded group scoring; truth is canonical Labels.truth(group, case_ids).

    Missing contrasts raise through the label callback and cannot become zero.
    Untestable contrasts and zero-positive directions remain explicitly undefined.
    Arrays may be memmaps; only one group's native rows are read at a time.
    """
    if (prediction.shape != observed_delta.shape or prediction.shape[0] != len(cases)
            or prediction.shape[1] != len(covered_columns) or cases.case_id.duplicated().any()
            or not np.isfinite(well_loss)):
        raise ValueError('invalid native-axis test inputs or test objective')
    pools = {m: PooledPR() for m in RETRIEVAL}
    rows = []
    for group, indices in cases.groupby('group_id', sort=False).indices.items():
        sums = [np.zeros(prediction.shape[1], dtype=np.float64) for _ in range(2)]
        for start in range(0, len(indices), 32):
            select = indices[start:start + 32]
            for total, matrix in zip(sums, (prediction, observed_delta), strict=True):
                values = matrix[select]
                if not np.isfinite(values).all():
                    raise ValueError('nonfinite native test expression')
                total += values.sum(axis=0, dtype=np.float64)
        predicted, observed = (total / len(indices) for total in sums)
        labels = truth(group, cases.iloc[indices].case_id.tolist())
        curves = retrieval(predicted, labels, covered_columns) if labels is not None else dict.fromkeys(RETRIEVAL)
        row = {'group_id': str(group), 'wells': len(indices), **continuous(predicted, observed)}
        for metric, curve in curves.items():
            row[metric] = None if curve is None else curve.ap
            row[f'{metric}_status'] = ('defined' if curve is not None else
                                      'untestable' if labels is None else 'zero_significant')
            pools[metric].add(curve)
        rows.append(row)
    metrics = {'test_loss': float(well_loss)}
    counts = {}
    for metric in METRICS[1:]:
        values = [r[metric] for r in rows if r[metric] is not None]
        metrics[metric] = float(np.mean(values)) if values else None
        counts[metric] = {'defined': len(values), 'undefined': len(rows) - len(values)}
    return {'metrics': metrics, 'groups': rows, 'counts': counts,
            'coverage': {'covered_genes': len(covered_columns), 'groups': len(rows), 'wells': len(cases)},
            'curves': {m: pool.frame().to_dict('records') for m, pool in pools.items()}}


def summarize_pr(curves):
    """Equal four-run interpolation, bounded by every attained endpoint."""
    if len(curves) != 4 or any(not curve for curve in curves):
        return []
    arrays = [pd.DataFrame(curve)[['recall', 'precision']].to_numpy(float) for curve in curves]
    for values in arrays:
        if (not np.isfinite(values).all() or np.any(np.diff(values[:, 0]) < 0)
                or np.any(values < 0) or np.any(values > 1)):
            raise ValueError('invalid pooled curve')
    endpoint = min(values[-1, 0] for values in arrays)
    grid = np.linspace(0, endpoint, 201) if endpoint else np.array([0.])
    samples = []
    for values in arrays:
        # Keep the last precision for each recall coordinate. This convention,
        # including the first precision at recall zero, is frozen in source data.
        unique = pd.DataFrame(values, columns=['recall','precision']).groupby('recall', sort=True).last()
        samples.append(np.interp(grid, unique.index, unique.precision))
    samples = np.asarray(samples)
    return [{'recall':float(x), 'mean':float(y), 'minimum':float(lo), 'maximum':float(hi),
             'contributors':4} for x,y,lo,hi in zip(grid,samples.mean(0),samples.min(0),samples.max(0),strict=True)]


def render(source, destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    destination = Path(destination)
    styles = ('-', '--', ':', '-.')
    fig, axes = plt.subplots(2,3,figsize=(16,9), constrained_layout=True)
    for model, ax in zip(MODELS, axes.flat, strict=True):
        for slot, run in enumerate(r for r in source['runs'] if r['model'] == model):
            offset = 0
            for phase in ('tahoe','tenv'):
                record = run['phases'][phase]
                points = record['history']
                x = np.array([p['update'] for p in points])
                train = [p['train_loss'] if 'train_loss' in p else float(np.dot(p['train_terms'], [.2,.003,.2,1.])) for p in points]
                ax.plot(x+offset,train,color='royalblue',ls=styles[slot],alpha=.8,lw=1)
                ax.plot(x+offset,[p['validation_loss'] for p in points],color='forestgreen',ls=styles[slot],lw=1)
                selected = next(p for p in points if p['update'] == record['best_update'])
                ax.plot(selected['update']+offset,selected['validation_loss'],'o',mfc='none',mec='forestgreen',ms=5)
                if phase == 'tahoe':
                    offset = int(x[-1]) + max(1,int(x[-1])*.04)
                    ax.axvline(offset,color='grey',ls='--',lw=.6,alpha=.35)
        ax.set(title=model.removeprefix('XPert_'), xlabel='Parameter updates (tenv shifted after Tahoe)', ylabel='Weighted objective')
        ax.grid(alpha=.15)
    legend = [Line2D([0],[0],color='royalblue',label='Window-mean training'),Line2D([0],[0],color='forestgreen',label='Full validation')]
    legend += [Line2D([0],[0],color='black',ls=s,label=f'fold {f}, seed {seed}') for s,(f,seed) in zip(styles,((0,17),(0,42),(1,17),(1,42)),strict=True)]
    fig.legend(handles=legend,loc='outside lower center',ncol=6)
    fig.savefig(destination/'r1_learning_curves.png',dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(2,4,figsize=(17,8),constrained_layout=True)
    for metric, ax in zip(METRICS,axes.flat,strict=True):
        for index,model in enumerate(MODELS):
            values = [r['metrics'][metric] for r in source['runs'] if r['model']==model]
            valid = [v for v in values if v is not None]
            if len(valid)==4:
                mean = np.mean(valid)
                ax.errorbar(index,mean,yerr=[[mean-min(valid)],[max(valid)-mean]],fmt='s',color='black',capsize=4)
            ax.scatter(index+np.linspace(-.12,.12,len(valid)),valid,color='steelblue',s=17)
            if len(valid)<4:
                ax.text(index,.98,f'{len(valid)}/4',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=8,color='firebrick')
        ax.set(title=metric,xticks=range(6),xticklabels=[m.removeprefix('XPert_') for m in MODELS])
        ax.tick_params(axis='x',labelrotation=35)
        ax.grid(axis='y',alpha=.2)
    fig.suptitle('Native-axis continuous metrics; axes differ across models. Group macro AUPRC; well-level test loss.\nSquares: equal four-run mean; bars: min–max; dots: individual runs. Undefined summaries are withheld.',fontsize=11)
    fig.savefig(destination/'r1_metrics.png',dpi=160)
    plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(15,9),constrained_layout=True)
    for model,ax in zip(MODELS,axes.flat,strict=True):
        for metric,color in (('AUPRC_negative','royalblue'),('AUPRC_positive','firebrick')):
            points=source['pr_summary'][model][metric]
            if points:
                table=pd.DataFrame(points)
                ax.plot(table.recall,table['mean'],color=color,label=metric.removeprefix('AUPRC_'))
                ax.fill_between(table.recall,table.minimum,table.maximum,color=color,alpha=.17)
            else:
                ax.text(.05,.9 if color=='royalblue' else .8,f'{metric}: undefined run',transform=ax.transAxes,color=color,fontsize=8)
        ax.set(title=model.removeprefix('XPert_'),xlabel='Recall on full CM truth',ylabel='Pooled precision',xlim=(0,1),ylim=(0,1.02))
        ax.grid(alpha=.2)
        if ax.lines: ax.legend()
    fig.suptitle('Four-run pooled PR mean and range; summaries stop at the minimum attained recall.\nPooled displayed area differs from group-macro AUPRC.',fontsize=11)
    fig.savefig(destination/'r1_PR_curves.png',dpi=160)
    plt.close(fig)


def publish(bundles, destination):
    """Publish a complete 24-run matrix to a fresh staging destination.

    Each bundle merges score_run output with model/fold/seed, checkpoint,
    status='complete' and phases from the selected run.json. The caller retains
    scratch prediction arrays elsewhere. Existing publications are never erased.
    """
    destination=Path(destination)
    expected={(m,f,s) for m in MODELS for f in (0,1) for s in (17,42)}
    identities=[(r['model'],r['fold'],r['seed']) for r in bundles]
    if len(identities)!=24 or set(identities)!=expected:
        raise ValueError('incomplete or duplicate 24-run matrix')
    runs=sorted(bundles,key=lambda r:(MODELS.index(r['model']),r['fold'],r['seed']))
    for run in runs:
        if run['status']!='complete' or not Path(run['checkpoint']).is_file():
            raise ValueError('incomplete selected checkpoint')
        if set(run['metrics'])!=set(METRICS):
            raise ValueError('missing test metrics')
        for value in run['metrics'].values():
            if value is not None and not np.isfinite(value):
                raise ValueError('nonfinite metric')
        for phase in ('tahoe','tenv'):
            history=run['phases'][phase]['history']
            if not history or run['phases'][phase]['best_update'] not in [p['update'] for p in history]:
                raise ValueError('missing phase selection or history')
            if any(not np.isfinite(p['validation_loss']) for p in history):
                raise ValueError('nonfinite training history')
        digest=sha(run['checkpoint'])
        if run.get('checkpoint_sha256',digest)!=digest:
            raise ValueError('selected checkpoint hash differs')
    summaries={m:{metric:summarize_pr([r['curves'].get(metric,[]) for r in runs if r['model']==m])
                  for metric in ('AUPRC_positive','AUPRC_negative')} for m in MODELS}
    destination.mkdir(parents=True,exist_ok=False)
    for entry in ('checkpoints','logs','figure_data'): (destination/entry).mkdir()
    stored=[]
    for run in runs:
        name=f"{run['model']}-fold-{run['fold']}-seed-{run['seed']}"
        item=dict(run,checkpoint=f'checkpoints/{name}.pt',checkpoint_sha256=sha(run['checkpoint']))
        shutil.copyfile(run['checkpoint'],destination/item['checkpoint'])
        write(destination/f'logs/{name}.json',item)
        stored.append(item)
    source={'schema':'step09-rerun-publication-v1','runs':stored,'pr_summary':summaries,
            'configuration':{'run_weighting':'equal; four required for mean/range',
                             'continuous_axis':'native model axis; cross-axis values are not paired-gene comparisons',
                             'test_loss_unit':'weighted native well-level author objective',
                             'continuous_unit':'native drug-dose group-mean delta; equal defined groups',
                             'AUPRC_unit':'exact canonical AP per group, equal defined groups on full-CM truth',
                             'PR_unit':'fixed-group pooled per run; four-run linear interpolation to minimum endpoint',
                             'recall_zero':'first measured precision; duplicate recall uses last precision',
                             'curve_coordinates':'raw phase update and exposure in histories; display shifts tenv after Tahoe',
                             'loss_weights':[.2,.003,.2,1.]}}
    write(destination/'figure_data/source.json',source)
    pd.DataFrame([{'model':r['model'],'fold':r['fold'],'seed':r['seed'],**r['metrics'],
                   'checkpoint_sha256':r['checkpoint_sha256']} for r in stored]).to_csv(destination/'figure_data/per_run.csv',index=False)
    aggregate=[]
    for model in MODELS:
        for metric in METRICS:
            values=[r['metrics'][metric] for r in stored if r['model']==model]
            valid=[v for v in values if v is not None]
            aggregate.append({'model':model,'metric':metric,'defined_runs':len(valid),'undefined_runs':4-len(valid),
                              'mean':float(np.mean(valid)) if len(valid)==4 else None,
                              'minimum':min(valid) if len(valid)==4 else None,'maximum':max(valid) if len(valid)==4 else None})
    pd.DataFrame(aggregate).to_csv(destination/'figure_data/summary.csv',index=False)
    render(source,destination)
    write(destination/'figure_data/figures.json',{str(p.relative_to(destination)):sha(p) for p in [*destination.glob('*.png'), *[destination/'figure_data'/name for name in ('source.json','per_run.csv','summary.csv')]]})
    return source


def replay(destination):
    """Regenerate all three PNGs from figure data and verify selected weights."""
    import tempfile
    destination=Path(destination)
    if {p.name for p in destination.iterdir()}!=ENTRIES:
        raise ValueError('publication top-level entries differ')
    source=json.loads((destination/'figure_data/source.json').read_text())
    expected=json.loads((destination/'figure_data/figures.json').read_text())
    identities={(r['model'],r['fold'],r['seed']) for r in source['runs']}
    if (len(source['runs'])!=24 or identities!={(m,f,s) for m in MODELS for f in (0,1) for s in (17,42)}
            or len(list((destination/'checkpoints').iterdir()))!=24):
        raise ValueError('incomplete retained checkpoint matrix')
    for name,digest in expected.items():
        if sha(destination/name)!=digest: raise ValueError('published source or figure changed')
    for run in source['runs']:
        if sha(destination/run['checkpoint'])!=run['checkpoint_sha256']:
            raise ValueError('selected checkpoint changed')
        name=Path(run['checkpoint']).stem
        if json.loads((destination/f'logs/{name}.json').read_text())!=run:
            raise ValueError('training log differs from figure provenance')
    with tempfile.TemporaryDirectory() as tmp:
        render(source,Path(tmp))
        if any(sha(Path(tmp)/name)!=digest for name,digest in expected.items() if name.endswith('.png')):
            raise ValueError('figure source replay differs')
    return True


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('publish','replay','score-publish'))
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--bundles',type=Path)
    parser.add_argument('--plan',type=Path)
    parser.add_argument('--work-root',type=Path)
    parser.add_argument('--predictions-root',type=Path)
    parser.add_argument('--loss-root',type=Path)
    parser.add_argument('--labels-root',type=Path)
    args=parser.parse_args()
    if args.command=='replay': replay(args.output)
    elif args.command=='score-publish':
        if any(getattr(args,k) is None for k in ('plan','work_root','predictions_root','loss_root')):
            parser.error('score-publish requires --plan --work-root --predictions-root --loss-root')
        publish(score_matrix(args.plan,args.work_root,args.predictions_root,args.loss_root,args.labels_root),args.output)
    else:
        if args.bundles is None: parser.error('publish requires --bundles JSON')
        publish(json.loads(args.bundles.read_text()),args.output)




def score_matrix(plan_path, work_root, predictions_root, loss_root, labels_root=None):
    """Verify selected run/prediction/target bindings and score a complete plan.

    Paths: work-root/fits/ID/{run.json,best.pt}, predictions-root/ID/
    {manifest.json,cases.csv,genes.csv,delta.npy}, loss-root/ID.json (the
    xpert.py evaluate-loss output). Significance labels use canonical Labels.
    """
    from stage2.evaluation.evaluate import Labels
    from rerun_queue import verify_plan
    plan=json.loads(Path(plan_path).read_text())
    verify_plan(plan)
    labels=Labels('tenv',labels_root)
    bundles=[]
    for item in plan['runs']:
        identity=item['id']
        fit=Path(work_root)/'fits'/identity
        checkpoint=fit/'best.pt'
        run=json.loads((fit/'run.json').read_text())
        prediction=Path(predictions_root)/identity
        manifest=json.loads((prediction/'manifest.json').read_text())
        loss=json.loads((Path(loss_root)/f'{identity}.json').read_text())
        digest=sha(checkpoint)
        predictor=Path(item['predictor'])
        evaluator=Path(item['evaluator'])
        pred_manifest=json.loads((predictor/'manifest.json').read_text())
        eval_manifest=json.loads((evaluator/'manifest.json').read_text())
        if (run['status']!='complete' or run['policy']['kind']!='formal'
                or any(run[key]!=item[key] for key in ('model','fold','seed','predictor_sha256'))
                or run['phases']['tenv']['checkpoint_sha256']!=digest
                or manifest['checkpoint_sha256']!=digest or loss['checkpoint_sha256']!=digest
                or manifest.get('role')!='test' or 'qualification_limit' in manifest
                or manifest['run_key']!=loss['run_key']
                or any(run.get(key)!=value for key,value in manifest['run_key'].items())
                or set(manifest['files'])!={'cases.csv','genes.csv','delta.npy'}
                or any(manifest['run_key'][key]!=item[key] for key in ('model','fold','seed','predictor_sha256'))
                or loss['predictor_sha256']!=item['predictor_sha256']
                or loss['evaluator_sha256']!=item['evaluator_sha256']):
            raise ValueError(f'{identity}: selected run/prediction/test loss differs')
        for name,digest_file in manifest['files'].items():
            if name not in ('cases.csv','genes.csv','delta.npy') or sha(prediction/name)!=digest_file:
                raise ValueError('prediction payload differs')
        prefix=f"tenv/fold_{item['fold']}"
        def bound(root,record,relative):
            path=root/relative
            if path.resolve()!=path.absolute() or sha(path)!=record['files'][relative]:
                raise ValueError(f'changed or unsafe evaluator input: {relative}')
            return path
        cases=pd.read_csv(prediction/'cases.csv',keep_default_na=False)
        native_cases=pd.read_csv(bound(predictor,pred_manifest,f'{prefix}/test.csv'),keep_default_na=False)
        if not cases.equals(native_cases[cases.columns]):
            raise ValueError('prediction cases differ from frozen test membership')
        genes=pd.read_csv(prediction/'genes.csv',keep_default_na=False)
        native_genes=pd.read_csv(bound(predictor,pred_manifest,'genes.csv'),keep_default_na=False)
        if not genes.equals(native_genes[genes.columns]):
            raise ValueError('prediction gene axis differs')
        columns=pd.Index(labels.genes.hgnc_id).get_indexer(genes.hgnc_id)
        if np.any(columns<0): raise ValueError('rerun axis must remain within CM')
        treated=np.load(bound(evaluator,eval_manifest,f'{prefix}/test_expression.npy'),mmap_mode='r')
        controls=np.load(bound(predictor,pred_manifest,'tenv/controls.npy'),mmap_mode='r')
        if (loss['test_expression_sha256']!=eval_manifest['files'][f'{prefix}/test_expression.npy']
                or loss['test_wells']!=len(cases) or loss['gene_count']!=len(genes)):
            raise ValueError('test objective cohort differs')
        class NativeDelta:
            shape=treated.shape
            def __getitem__(self,rows):
                return treated[rows]-controls[cases.control_index.to_numpy()[rows]]
        scored=score_run(cases,np.load(prediction/'delta.npy',mmap_mode='r'),NativeDelta(),columns,
                         labels.truth,loss['test_loss'])
        bundles.append({**scored,'model':item['model'],'fold':item['fold'],'seed':item['seed'],
                        'status':'complete','checkpoint':str(checkpoint),'checkpoint_sha256':digest,
                        'phases':run['phases'],'training_record':run,'provenance':{'plan_sha256':sha(plan_path),
                        'run_sha256':sha(fit/'run.json'),'prediction_manifest_sha256':sha(prediction/'manifest.json'),
                        'test_loss_sha256':sha(Path(loss_root)/f'{identity}.json'),
                        'significance_sha256':labels.manifest_digest,'predictor_sha256':item['predictor_sha256'],
                        'evaluator_sha256':item['evaluator_sha256']}})
    labels.verify()
    return bundles


if __name__=='__main__': main()
