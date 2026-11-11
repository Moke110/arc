"""Bounded real-population HVG isolation and measurement-derived diagnostic QA."""
import json
from pathlib import Path
import tempfile

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem

from .hvg import HVGRelease, build_hvg, population_hash, sha256, verify_hvg
from .metrics import score_groups, summarize_runs
from .plots import render_retrieval, replay_retrieval


def verify_population_isolation(training, hvg, *, progress=None):
    """Intervene on every held-out row through VDS; training values stay exact.

    Each fold gets a temporary virtual h5ad. Only training rows reference
    original X; held-out rows become a legal single-gene CP10K spike. No whole
    expression matrix is copied or materialized. All train statistics and ranks
    must remain bit-exact while descriptive statistics must change.
    """
    training = Path(training).resolve()
    release = HVGRelease(hvg)
    before = sha256(training)
    if before != release.manifest['training_sha256']:
        raise ValueError('HVG isolation training seal differs')
    genes = release.manifest['gene_count']
    if genes < 2:
        raise ValueError('HVG intervention requires at least two genes')
    results = {}
    with h5py.File(training, 'r') as original:
        if not isinstance(original['X'],h5py.Dataset):
            raise ValueError('HVG intervention requires the published dense chunked X dataset')
        obs, var = read_elem(original['obs']), read_elem(original['var'])
        positions = pd.Series(np.arange(len(obs)), index=obs.sample_id.astype(str))
        for fold, assignment in release.cohorts.groupby('fold',sort=True):
            name = f'fold_{int(fold)}'
            train_ids = assignment.loc[assignment.role.eq('train'),'sample_id'].astype(str).tolist()
            heldout_ids = assignment.loc[~assignment.role.eq('train'),'sample_id'].astype(str).tolist()
            train_rows = np.sort(positions.loc[train_ids].to_numpy())
            heldout_rows = np.sort(positions.loc[heldout_ids].to_numpy())
            if not len(heldout_rows):
                raise ValueError('HVG intervention has no held-out rows')
            spike_gene = int(np.argmin(original['X'][heldout_rows[0]]))
            spike = np.zeros(genes,dtype=np.float32)
            spike[spike_gene] = np.log1p(np.float32(10000))
            if np.array_equal(original['X'][heldout_rows[0]],spike):
                raise ValueError('held-out intervention did not change expression')
            with tempfile.TemporaryDirectory(prefix='bases-hvg-isolation-') as directory:
                root = Path(directory)
                with h5py.File(root/'spike.h5','w') as replacement:
                    replacement.create_dataset('value',data=[spike[spike_gene]],dtype='float32')
                layout = h5py.VirtualLayout(shape=original['X'].shape,dtype='float32')
                source = h5py.VirtualSource(str(training),'X',shape=original['X'].shape)
                replacement = h5py.VirtualSource(str(root/'spike.h5'),'value',shape=(1,))
                for row in train_rows:
                    layout[int(row),:] = source[int(row),:]
                for row in heldout_rows:
                    layout[int(row),spike_gene] = replacement[0]
                mutated = root/'heldout-intervention.h5ad'
                with h5py.File(mutated,'w') as file:
                    write_elem(file,'obs',obs)
                    write_elem(file,'var',var)
                    file.create_virtual_dataset('X',layout,fillvalue=np.float32(0))
                build_hvg(mutated,root/'hvg',cohorts=assignment,expected_genes=genes)
                changed = HVGRelease(root/'hvg')
                pd.testing.assert_frame_equal(release.ranking(name),changed.ranking(name),check_exact=True)
                if np.array_equal(release.ranking('descriptive').log_mean,
                                  changed.ranking('descriptive').log_mean,equal_nan=True):
                    raise ValueError('descriptive rank statistics ignored held-out intervention')
                np.testing.assert_array_equal(
                    release.select(name,min(2000,genes),training_sample_ids=train_ids),
                    changed.select(name,min(2000,genes),training_sample_ids=train_ids))
                for population, ids in [('descriptive',train_ids),(name,train_ids[:-1]+heldout_ids[:1])]:
                    try:
                        release.select(population,1,training_sample_ids=ids)
                    except ValueError:
                        pass
                    else:
                        raise ValueError('HVG selector accepted descriptive or held-out population')
            results[name] = {'train_wells':len(train_ids),'heldout_wells':len(heldout_ids),
                             'train_population_sha256':population_hash(train_ids),
                             'train_statistics_bit_exact':True,'descriptive_statistics_changed':True,
                             'selector_rejects_descriptive_and_wrong_population':True,
                             'intervention':'Every held-out X row replaced with single-gene CP10K; exact training X rows referenced through read-only HDF5 VDS'}
            if progress:
                progress({'phase':'hvg_population_isolation','fold':int(fold),'train_wells':len(train_ids),'heldout_wells':len(heldout_ids),'result':'PASS'})
    if sha256(training) != before:
        raise ValueError('training changed during HVG population verification')
    return {'result':'PASS','training_sha256':before,'folds':results}


class _Rows:
    """Bounded prediction-array protocol over a selected h5ad layer population."""
    ndim = 2
    def __init__(self, matrix, rows):
        self.matrix, self.rows = matrix, np.asarray(rows,dtype=int)
        self.shape = (len(rows),matrix.shape[1])

    def __getitem__(self, rows):
        source = self.rows[rows]
        if len(source) and np.any(np.diff(source) <= 0):
            raise ValueError('diagnostic source rows must have unique increasing identity')
        return self.matrix[source]


def verify_hvg_metrics(training, hvg, successor, output, *, expected_groups=1377, progress=None):
    """Actual-data QA only: measured delta ranked against its statistical truth.

    This deliberately evaluates no model. Exact test-role populations across
    the new folds must cover every group once; neither subsets nor duplicated
    contrasts are accepted. SuccessorTruth supplies the accepted v2 binding.
    """
    from .truth_successor import SuccessorTruth, check_successor
    training, hvg, successor, output = map(Path,(training,hvg,successor,output))
    check_successor(successor,training)
    payload = training/'tenv.h5ad'
    release = HVGRelease(hvg)
    if (release.manifest['gene_count'] != 18164
            or sorted(release.cohorts.fold.unique().tolist()) != [0,1,2]):
        raise ValueError('diagnostic requires the full 18164 axis and three accepted folds')
    hvg_replay = verify_hvg(payload,hvg)
    if progress:
        progress({'phase':'hvg_source_replay','result':'PASS'})
    isolation = verify_population_isolation(payload,hvg,progress=progress)
    output.mkdir(parents=True,exist_ok=False)
    runs, summaries, covered_samples, covered_groups, inputs = [], [], set(), set(), None
    with h5py.File(payload,'r') as file:
        obs = read_elem(file['obs']).reset_index(drop=True)
        if 'delta' not in file['layers'] or file['layers/delta'].shape != file['X'].shape:
            raise ValueError('diagnostic requires the complete measured delta layer')
        if len(obs) != 4206:
            raise ValueError('diagnostic requires all 4206 measured treated wells')
        positions = pd.Series(np.arange(len(obs)),index=obs.sample_id.astype(str))
        columns = np.arange(release.manifest['gene_count'])
        for fold, assignment in release.cohorts.groupby('fold',sort=True):
            labels = SuccessorTruth(successor,training,assignment,'test')
            binding = {'truth_manifest_sha256':sha256(labels.root/'manifest.json'),
                       'successor_manifest_sha256':labels.successor_manifest_sha256}
            if inputs is not None and binding != inputs:
                raise ValueError('diagnostic truth binding differs across folds')
            inputs = binding
            samples = assignment.loc[assignment.role.eq('test'),['sample_id','drug_parent','dose_um']].copy()
            samples = samples.assign(_position=samples.sample_id.map(positions)).sort_values('_position')
            if samples._position.isna().any() or samples.sample_id.duplicated().any():
                raise ValueError('diagnostic role population differs from training')
            ids = set(samples.sample_id.astype(str))
            if covered_samples & ids:
                raise ValueError('diagnostic test population repeats across folds')
            expected_ids = set(assignment.loc[assignment.role.eq('test'),'sample_id'].astype(str))
            if ids != expected_ids:
                raise ValueError('diagnostic must cover the complete declared test role')
            predictions = _Rows(file['layers/delta'],samples._position.to_numpy(dtype=int))
            result = score_groups(predictions,samples.drop(columns='_position'),labels,columns)
            keys = {(row.drug_parent,format(float(row.dose_um),'.12g')) for row in result['groups'].itertuples()}
            if covered_groups & keys:
                raise ValueError('diagnostic test contrast repeats across folds')
            covered_groups.update(keys)
            covered_samples.update(ids)
            # Public reader rejects a complete group belonging to a different role.
            forbidden = assignment[assignment.role.eq('train')].groupby(['drug_parent','dose_um'],sort=True)
            (parent,dose), group = next(iter(forbidden))
            try:
                labels.truth(parent,dose,group.sample_id.astype(str).tolist())
            except ValueError:
                pass
            else:
                raise ValueError('diagnostic evaluator accessed training-role truth')
            result['groups'].to_csv(output/f'fold_{int(fold)}_groups.csv',index=False)
            runs.append({'label':f'Diagnostic fold {int(fold)}','summary':result['summary'],'curves':result['curves']})
            summaries.append({'fold':int(fold),**result['summary']})
            if progress:
                progress({'phase':'measurement_delta_diagnostic','fold':int(fold),'samples':len(ids),'groups':len(keys),'result':'PASS'})
        expected = {(str(row.drug_parent),format(float(row.dose_um),'.12g')) for row in obs.itertuples()}
        if (covered_samples != set(obs.sample_id.astype(str)) or covered_groups != expected
                or len(covered_groups) != expected_groups):
            raise ValueError('diagnostic folds do not cover every measured sample and contrast exactly once')
    panel = render_retrieval(runs,output/'diagnostic_pr.png',diagnostic=True)
    replay_retrieval(output/'diagnostic_pr.json',output/'diagnostic_pr_replayed.png')
    if sha256(output/'diagnostic_pr.png') != sha256(output/'diagnostic_pr_replayed.png'):
        raise ValueError('diagnostic PR byte replay differs')
    release.cohorts.to_csv(output/'cohorts.csv',index=False)
    payload_sha = sha256(payload)
    if payload_sha != release.manifest['training_sha256']:
        raise ValueError('diagnostic training source changed')
    report = {'schema':'stage2-hvg-metrics-acceptance-v1','result':'NUMERIC_PASS_VISUAL_PENDING','numeric_result':'PASS','diagnostic':True,
              'interpretation':'Measurement-derived mechanical QA only; no model evaluated and no model performance claim',
              'training_sha256':payload_sha,'hvg_manifest_sha256':sha256(hvg/'manifest.json'),
              **inputs,'hvg_replay':hvg_replay,'population_isolation':isolation,
              'samples':len(covered_samples),'groups':len(covered_groups),'folds':summaries,
              'equal_run_macro':summarize_runs(summaries),'panel':panel,
              'png_byte_replay':True,'visual_qa':'PENDING manual inspection',
              'files':{p.name:sha256(p) for p in sorted(output.iterdir()) if p.is_file()}}
    (output/'acceptance.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    return report
