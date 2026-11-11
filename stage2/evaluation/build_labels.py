"""Build bounded, checkpointed Stage 2 DESeq2 labels from native counts."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import h5py
from anndata.io import read_elem
import numpy as np
import pandas as pd

REPO = next(path for path in Path(__file__).resolve().parents if (path / '.git').exists())
STEP09 = REPO / 'stage2/09_xpert_optimize'
sys.path.insert(0, str(STEP09 / 'execution/script'))
from build_significance import fit_shards, write_statistic
sys.path.insert(0, str(STEP09 / 'execution/01_corpus'))
from build import raw_batches


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def membership(ids):
    return hashlib.sha256(json.dumps(sorted(ids), separators=(',', ':')).encode()).hexdigest()


def significance_gene_family(step09=STEP09):
    """Read the preserved CM-plus-HDGFL3 family across Round 1 replacement."""
    active = step09 / 'output/round1/labels/test/genes.csv'
    archived = step09 / 'execution/history/round1-exploratory/round1/labels/test/genes.csv'
    return pd.read_csv(active if active.is_file() else archived)


def prepare(dataset):
    native = REPO / f'stage2/02_dataset_corpus/output/corpus/v2/{dataset}.h5ad'
    genes = significance_gene_family()
    with h5py.File(native) as f:
        obs, var = read_elem(f['obs']), read_elem(f['var'])
    if dataset == 'tenv':
        cases = pd.read_csv(REPO / 'stage2/03_modes_tenv/output/input/wells.csv')
        columns = var.reset_index().set_index('hgnc_id').index.get_indexer(genes.hgnc_id)
    else:
        release = REPO / 'stage2/04_pdiff_hyperparameter/execution/output/tahoe-smoke-ticket68/release'
        cases = pd.concat([pd.read_csv(release / f'fold_{i}/test_cases.csv') for i in range(3)])
        axis = pd.read_csv(REPO / 'stage2/03_modes_tenv/output/input/genes.csv')
        extra = var.loc[var.index == 'HDGFL3', 'ensembl_id']
        if len(extra) != 1:
            raise ValueError('expected HDGFL3 native Tahoe gene')
        ensembl = [*axis.ensembl_id, extra.iloc[0]]
        columns = var.set_index('ensembl_id').index.get_indexer(ensembl)
    cases = cases.drop_duplicates('case_id').copy()
    cases['source_row'] = obs.index.get_indexer(cases.case_id)
    if (cases.source_row < 0).any() or cases.case_id.duplicated().any():
        raise ValueError('invalid native joins')
    rows = obs.iloc[cases.source_row]
    if rows.is_control.any() or not np.array_equal(rows.dose_um, cases.dose_um):
        raise ValueError('treated metadata mismatch')
    cases['plate'] = rows.plate.to_numpy()
    if dataset == 'tahoe':
        cases['cell_line_id'] = rows.cell_line_id.to_numpy()
    return native, genes, obs, columns, cases


def shard_inputs(path, matrix, obs, columns, cases, dataset):
    records, samples = [], []
    with (path / 'counts.bin').open('wb') as stream:
        for group, part in cases.groupby('group_id', sort=True):
            plates = part.plate.unique()
            mask = obs.is_control & obs.plate.isin(plates)
            if dataset == 'tenv':
                mask &= obs.compound_name.str.upper().eq('DMSO') & ~obs.is_blank
            else:
                if part.cell_line_id.nunique() != 1:
                    raise ValueError('mixed cell lines in Tahoe contrast')
                mask &= obs.drug.eq('DMSO_TF') & obs.cell_line_id.eq(part.cell_line_id.iloc[0])
            controls = np.flatnonzero(mask)
            if not len(controls) or set(obs.iloc[controls].plate) != set(plates):
                raise ValueError(f'missing matched controls: {group}')
            batch = pd.concat([part[['case_id','source_row','plate']].assign(condition='treated'),
                pd.DataFrame({'case_id':obs.index[controls], 'source_row':controls,
                              'plate':obs.iloc[controls].plate.to_numpy(),'condition':'control'})]).sort_values('source_row')
            if batch.source_row.duplicated().any():
                raise ValueError('duplicated contrast well')
            records.append({'group_id':group,'n_treated':len(part),'n_controls':len(controls),
                            'n_plates':len(plates),'case_ids_sha256':membership(part.case_id.tolist())})
            samples.append(batch.assign(group_id=group))
            present = columns >= 0
            for _, raw in raw_batches(matrix,batch.source_row.to_numpy(),columns[present]):
                values=np.zeros((len(raw),len(columns)),dtype='<i4')
                if np.any(raw > np.iinfo(np.int32).max): raise ValueError('count overflow')
                values[:,present]=raw
                values.tofile(stream)
    groups=pd.DataFrame(records);groups.to_csv(path/'groups.csv',index=False)
    pd.concat(samples).to_csv(path/'samples.csv',index=False)
    return groups


def build(dataset, output, rscript, limit=None, batch_size=50, reuse=None, upgrade_checkpoint=False):
    final = output
    if (final / 'manifest.json').exists():
        raise FileExistsError(final)
    output = final.with_name(final.name + '.building')
    output.mkdir(parents=True,exist_ok=True)
    native,genes,obs,columns,cases=prepare(dataset)
    identity = {'dataset': dataset, 'native_sha256': digest(native),
                'builder_sha256': digest(Path(__file__)),
                'worker_sha256': digest(STEP09/'execution/script/deseq2_contrasts.R'),
                'case_memberships': membership(cases.case_id.tolist()), 'limit': limit,
                'batch_size': batch_size}
    pin = output / 'identity.json'
    if pin.exists() and json.loads(pin.read_text()) != identity:
        previous = json.loads(pin.read_text())
        unchanged = all(previous.get(k) == v for k,v in identity.items() if k != 'builder_sha256')
        if not upgrade_checkpoint or not unchanged:
            raise ValueError('checkpoint input/code identity differs')
        (output / 'previous-identity.json').write_text(json.dumps(previous, sort_keys=True))
    pin.write_text(json.dumps(identity, sort_keys=True))
    groups=sorted(cases.group_id.unique())
    if limit: groups=groups[:limit]
    reused = None
    if reuse:
        reuse_manifest = json.loads((reuse/'manifest.json').read_text())
        if any(digest(reuse/name)!=value for name,value in reuse_manifest['files'].items()):
            raise ValueError('reused label payload changed')
        pd.testing.assert_frame_equal(genes, pd.read_csv(reuse/'genes.csv'))
        reused = pd.read_csv(reuse/'groups.csv',keep_default_na=False)
        expected = cases.groupby('group_id').case_id.apply(lambda x: membership(x.tolist()))
        if not all(expected[g] == h for g,h in zip(reused.group_id,reused.case_ids_sha256)):
            raise ValueError('reused label memberships changed')
        # Preserve already-started first batch; skip other verified cached groups.
        cached = set(reused.group_id)
        groups = groups[:batch_size] + [g for g in groups[batch_size:] if g not in cached]
    def process_part(start):
        dest=output/f'parts/{start//batch_size:05d}'
        manifest=dest/'manifest.json'
        if not manifest.exists():
            dest.mkdir(parents=True,exist_ok=True)
            with h5py.File(native) as f:
                metadata=shard_inputs(dest,f['X'],obs,columns,cases[cases.group_id.isin(groups[start:start+batch_size])],dataset)
            fit_shards(rscript,[(dest,metadata)],len(genes))
            pd.concat([metadata,pd.read_csv(dest/'results.csv',keep_default_na=False)],axis=1).to_csv(dest/'metadata.csv',index=False)
            (dest/'counts.bin').unlink()
            manifest.write_text(json.dumps({'groups':metadata.group_id.tolist(),'files':{p.name:digest(p) for p in dest.iterdir() if p.is_file() and p.name!='manifest.json'}}))
        saved=json.loads(manifest.read_text())
        if any(digest(dest/name)!=value for name,value in saved['files'].items()):
            raise ValueError('checkpoint payload checksum differs')
        metadata=pd.read_csv(dest/'metadata.csv',keep_default_na=False)
        if metadata.group_id.tolist()!=groups[start:start+batch_size]: raise ValueError('resume mismatch')
        return dest, metadata
    with ThreadPoolExecutor(max_workers=2) as workers:
        shards = list(workers.map(process_part, range(0,len(groups),batch_size)))
    if reused is not None:
        keep = ~reused.group_id.isin(groups)
        rows = np.flatnonzero(keep.to_numpy())
        part = output/'reused';part.mkdir(exist_ok=True)
        metadata = reused.loc[keep].reset_index(drop=True)
        for field in ['padj','log2_fold_change','pvalue']:
            original = np.load(reuse/f'{field}.npy',mmap_mode='r')
            with (part/f'{field}.bin').open('wb') as stream:
                for index in rows: original[index].astype('<f8',copy=False).tofile(stream)
        samples = pd.read_csv(reuse/'contrast_samples.csv')
        samples[samples.group_id.isin(metadata.group_id)].to_csv(part/'samples.csv',index=False)
        shards.append((part,metadata))
    pd.concat([g for _,g in shards]).to_csv(output/'groups.csv',index=False)
    genes.to_csv(output/'genes.csv',index=False)
    for field in ['padj','log2_fold_change','pvalue']:
        write_statistic(shards,field,output/f'{field}.npy',len(genes))
    pd.concat([pd.read_csv(p/'samples.csv') for p,_ in shards]).to_csv(output/'contrast_samples.csv',index=False)
    payload={'schema':'stage2-significance-v1','dataset':dataset,'fdr':.05,'gene_count':len(genes),
             'native_sha256':digest(native),'method':'two-sided DESeq2 Wald; within-contrast BH over CM + HDGFL3',
             'worker_sha256':digest(STEP09/'execution/script/deseq2_contrasts.R'),
             'files':{p.name:digest(p) for p in output.iterdir() if p.is_file() and p.name!='manifest.json'}}
    payload['complete'] = limit is None
    payload['expected_groups'] = int(cases.group_id.nunique())
    payload['reuse_manifest_sha256'] = digest(reuse/'manifest.json') if reuse else None
    payload['identity'] = identity
    (output/'manifest.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    output.rename(final)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset',choices=['tenv','tahoe'],required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--rscript',default='/home/chang/miniconda3/envs/biomni_e1/bin/Rscript')
    parser.add_argument('--limit',type=int)
    parser.add_argument('--reuse',type=Path)
    parser.add_argument('--upgrade-checkpoint',action='store_true')
    args=parser.parse_args();build(args.dataset,args.output,args.rscript,args.limit,reuse=args.reuse,upgrade_checkpoint=args.upgrade_checkpoint)

if __name__=='__main__':main()
