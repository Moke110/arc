"""Metadata-only binding/structure audit for stage2/01_bases products (#155).

Reads manifests, non-expression products and existing build records only.
Never opens *.h5ad, *.npy, *.npz or *.bin payloads.
"""
import argparse, hashlib, json, csv
from pathlib import Path

_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--step-root', type=Path, default=Path('/home/chang/arc/stage2/01_bases'),
                     help='bases step root to audit (read-only)')
_parser.add_argument('--write', type=Path, help='write the audit JSON to this path')
_args = _parser.parse_args()
STEP = _args.step_root
OUT = STEP / 'output'
REC = STEP / 'execution/build_record'
EXPR = {'.h5ad', '.npy', '.npz', '.bin'}
checks, structure = [], []
manifests = {}

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def add(target, name, status, detail):
    target.append({'name': name, 'status': status, 'detail': detail})

def digest_map(target, base, declared):
    matched, mismatch, absent, skipped, problems = 0, 0, 0, 0, []
    for name, digest in declared.items():
        if Path(name).suffix in EXPR:
            skipped += 1
            continue
        if digest is None:
            continue
        f = base / name
        if not f.is_file():
            absent += 1
            problems.append({'file': name, 'problem': 'declared but absent'})
            continue
        actual = sha(f)
        if actual == digest:
            matched += 1
        else:
            mismatch += 1
            problems.append({'file': name, 'declared': digest, 'actual': actual})
    add(target, f'{base.name}: declared non-expression digests',
        'PASS' if mismatch == 0 and absent == 0 else 'FAIL',
        {'matched': matched, 'mismatched': mismatch, 'absent': absent,
         'expression_skipped': skipped, 'problems': problems[:6]})
    return mismatch + absent

def load(rel):
    p = OUT / rel
    value = json.loads(p.read_text())
    manifests[rel] = sha(p)
    return value

idn = load('identity-v1/manifest.json')
add(checks, 'identity-v1: genes_sha256',
    'PASS' if sha(OUT / 'identity-v1/genes.csv') == idn['genes_sha256'] else 'FAIL',
    {'genes_sha256': idn['genes_sha256']})
digest_map(checks, OUT / 'native-v1/tables',
           {n: v['sha256'] for n, v in idn['inputs'].items()
            if n in ('genes.csv', 'gene_crosswalk.csv')})
digest_map(checks, OUT / 'shared-v1',
           {n: v['sha256'] for n, v in idn['inputs'].items()
            if n in ('tahoe_mapping.csv', 'tenv_mapping.csv')})
_annotation = STEP.parents[1] / 'pert_data/Annotation/source/hgnc_complete_set.txt'
add(checks, 'identity-v1: HGNC annotation input seal',
    'PASS' if sha(_annotation) == idn['inputs']['hgnc_complete_set.txt']['sha256'] else 'FAIL',
    {'path': str(_annotation)})

nat = load('native-v1/tenv_manifest.json')
digest_map(checks, OUT / 'native-v1', nat['files'])

sh = load('shared-v2/manifest.json')
digest_map(checks, OUT / 'shared-v2', dict(sh['files']))
digest_map(checks, OUT / 'native-v1/tables',
           {n: d for n, d in sh['inputs'].items() if n in ('genes.csv', 'gene_crosswalk.csv')})

tk = load('tahoe-counts-v1/manifest.json')
add(checks, 'tahoe-counts-v1: mapping_sha256 -> token_mapping.csv',
    'PASS' if sha(OUT / 'tahoe-counts-v1/token_mapping.csv') == tk['mapping_sha256'] else 'FAIL',
    {'mapping_sha256': tk['mapping_sha256'],
     'source_shards_declared': len(tk['inputs']['source_files']),
     'source_shard_rehash': 'not run (immutable external source)'})

for name in ('training-v1', 'training-v2'):
    digest_map(checks, OUT / name, load(f'{name}/manifest.json')['files'])

tr1 = load('truth-tenv-v1/manifest.json')
digest_map(checks, OUT / 'truth-tenv-v1', dict(tr1['files']))

tr2 = load('truth-tenv-v2/manifest.json')
digest_map(checks, OUT / 'truth-tenv-v1', dict(tr2['source_truth_files']))
add(checks, 'truth-tenv-v2: worker_sha256 -> deseq2_contrasts.R',
    'PASS' if sha(OUT / 'truth-tenv-v2/deseq2_contrasts.R') == tr2['worker_sha256'] else 'FAIL',
    {'worker_sha256': tr2['worker_sha256']})

hv = load('hvg-tenv-v1/manifest.json')
hvg_declared = {'cohorts.csv': hv['cohorts']['sha256']}
for pop in hv['populations'].values():
    hvg_declared[pop['ranking']] = pop['sha256']
digest_map(checks, OUT / 'hvg-tenv-v1', hvg_declared)

def cross(name, status, detail):
    add(checks, f'cross: {name}', status, detail)

def same(label, expected, actual):
    cross(label, 'PASS' if expected == actual else 'FAIL',
          {'expected': expected, 'actual': actual})

add(checks, 'cross: training-v1.inputs.shared_manifest -> shared-v1 manifest seal',
    'EXPECTED_GAP',
    {'declared': load('training-v1/manifest.json')['inputs']['shared_manifest'],
     'detail': 'shared-v1 has no manifest.json: its expression payload and seal were deleted with the '
               'historical release; training-v1 stays a historical seal, not a reconstruction input'})
same('training-v2.inputs.shared_manifest -> shared-v2/manifest.json',
     load('training-v2/manifest.json')['inputs']['shared_manifest'],
     manifests['shared-v2/manifest.json'])
same('truth-v1.provenance.training_manifest_sha256 -> training-v1/manifest.json',
     tr1['provenance']['training_manifest_sha256'], manifests['training-v1/manifest.json'])
same('truth-v1.provenance.native_counts_sha256 -> native manifest tenv.h5ad',
     tr1['provenance']['native_counts_sha256'], nat['files']['tenv.h5ad'])
same('truth-v1.provenance.tenv_training_sha256 -> training-v1 manifest tenv.h5ad',
     tr1['provenance']['tenv_training_sha256'], load('training-v1/manifest.json')['files']['tenv.h5ad'])
same('truth-v1.provenance.worker_sha256 -> truth-v2 deseq2_contrasts.R',
     tr1['provenance']['worker_sha256'], tr2['worker_sha256'])
same('truth-v1 manifest group_count -> groups.csv CSV records',
     tr1['group_count'], len(list(csv.DictReader(open(OUT / 'truth-tenv-v1/groups.csv')))))
same('truth-complete-verify.manifest_sha256 -> truth-v1/manifest.json',
     json.loads((REC / 'truth_complete_verify.json').read_text())['manifest_sha256'],
     manifests['truth-tenv-v1/manifest.json'])
same('truth-successor-acceptance.manifest_sha256 -> truth-v2/manifest.json',
     json.loads((REC / 'truth_successor_acceptance.json').read_text())['manifest_sha256'],
     manifests['truth-tenv-v2/manifest.json'])
same('hvg manifest training_sha256 -> training-v2 tenv.h5ad (manifest-to-manifest)',
     hv['training_sha256'], load('training-v2/manifest.json')['files']['tenv.h5ad'])
same('hvg manifest cohorts sha -> diagnostic cohorts.csv',
     hv['cohorts']['sha256'], sha(OUT / 'hvg-metrics-diagnostic-v1/cohorts.csv'))
acc = json.loads((REC / 'hvg_metrics_actual_acceptance.json').read_text())['acceptance']
same('hvg acceptance.hvg_manifest_sha256 -> hvg-tenv-v1/manifest.json',
     acc['hvg_manifest_sha256'], manifests['hvg-tenv-v1/manifest.json'])
same('hvg acceptance.successor_manifest_sha256 -> truth-v2/manifest.json',
     acc['successor_manifest_sha256'], manifests['truth-tenv-v2/manifest.json'])
same('hvg acceptance.truth_manifest_sha256 -> truth-v1/manifest.json',
     acc['truth_manifest_sha256'], manifests['truth-tenv-v1/manifest.json'])
same('hvg acceptance.training_sha256 -> training-v2 manifest tenv.h5ad',
     acc['training_sha256'], load('training-v2/manifest.json')['files']['tenv.h5ad'])
digest_map(checks, OUT / 'hvg-metrics-diagnostic-v1', acc['files'])

binding_rows = []
for key, item in tr2['bindings'].items():
    path = Path(item['path'])
    row = {'binding': key, 'recorded_path': item['path'],
           'path_is_resolved_absolute': str(path.resolve()) == str(path),
           'exists': path.exists(), 'suffix': path.suffix}
    if path.is_file() and path.suffix not in EXPR:
        row['digest_matches'] = sha(path) == item['sha256']
    elif path.is_file():
        row['digest_matches'] = 'not rehashed (expression payload)'
    binding_rows.append(row)
add(checks, 'truth-tenv-v2: recorded binding paths are live absolute dependencies',
    'PASS' if all(r['path_is_resolved_absolute'] and r['exists'] for r in binding_rows) else 'FAIL',
    binding_rows)

def rows(path):
    return list(csv.DictReader(open(path)))

axis = rows(OUT / 'identity-v1/genes.csv')
add(structure, 'identity axis: 18164 rows, gene_index 0..18163, unique symbol/HGNC',
    'PASS' if (len(axis) == 18164 and [int(r['gene_index']) for r in axis] == list(range(18164))
               and len({r['symbol'] for r in axis}) == 18164
               and len({r['hgnc_id'] for r in axis}) == 18164) else 'FAIL',
    {'rows': len(axis), 'index_5091': axis[5091]['symbol']})
native_axis = rows(OUT / 'native-v1/tables/genes.csv')
add(structure, 'frozen axis order identical in native tables and identity product',
    'PASS' if [r['symbol'] for r in native_axis] == [r['symbol'] for r in axis] else 'FAIL',
    {'rows': len(native_axis),
     'native_tables_equals_shared_v2_bytes': sha(OUT / 'native-v1/tables/genes.csv') == sha(OUT / 'shared-v2/genes.csv')})

groups = rows(OUT / 'truth-tenv-v1/groups.csv')
ids = [g['group_id'] for g in groups]
add(structure, 'truth-v1 groups: 1377 unique group ids',
    'PASS' if len(groups) == 1377 and len(set(ids)) == 1377 else 'FAIL',
    {'rows': len(groups), 'statuses': sorted({g['status'] for g in groups}),
     'sum_treated': sum(int(g['n_treated']) for g in groups),
     'sum_controls': sum(int(g['n_controls']) for g in groups),
     'tested_empty': sum(1 for g in groups if g['status'] == 'eligible' and int(g['n_significant']) == 0),
     'multiplate': sum(1 for g in groups if int(g['n_plates']) > 1)})
group_files = list((OUT / 'truth-tenv-v1/groups').glob('*.json'))
add(structure, 'truth-v1 per-group payloads cover exactly the group table',
    'PASS' if {p.stem for p in group_files} == set(ids) else 'FAIL',
    {'payloads': len(group_files)})
cs = rows(OUT / 'truth-tenv-v1/contrast_samples.csv')
add(structure, 'truth-v1 contrast membership: 4206 treated wells plus controls',
    'PASS' if sum(1 for r in cs if r['condition'] == 'treated') == 4206 else 'FAIL',
    {'rows': len(cs), 'treated': sum(1 for r in cs if r['condition'] == 'treated'),
     'control': sum(1 for r in cs if r['condition'] == 'control')})

coh = rows(OUT / 'hvg-tenv-v1/cohorts.csv')
per_fold = {}
for r in coh:
    per_fold.setdefault(r['fold'], {}).setdefault(r['role'], 0)
    per_fold[r['fold']][r['role']] += 1
roles_per_sample = {}
for r in coh:
    roles_per_sample.setdefault(r['sample_id'], set()).add((r['fold'], r['role']))
complete = all(len(v) == 3 and sum(1 for f, role in v if role == 'test') == 1
               for v in roles_per_sample.values())
add(structure, 'HVG cohorts: every well in exactly 3 folds with exactly one test role',
    'PASS' if complete and len(roles_per_sample) == 4206 else 'FAIL',
    {'samples': len(roles_per_sample), 'rows': len(coh), 'per_fold_role': per_fold})
for pop in ('descriptive', 'fold_0', 'fold_1', 'fold_2'):
    rank_rows = rows(OUT / f'hvg-tenv-v1/{pop}.csv')
    ranks = sorted(int(r['rank']) for r in rank_rows)
    add(structure, f'hvg-tenv-v1/{pop}: rank permutation over 18164 genes',
        'PASS' if ranks == list(range(1, 18165)) else 'FAIL',
        {'rows': len(rank_rows), 'min': ranks[0], 'max': ranks[-1]})

pr = rows(STEP / 'diagnostic_tenv_pr.csv')
add(structure, 'diagnostic PR panel CSV: three fold labels x three metrics',
    'PASS' if (sorted({r['label'] for r in pr}) == ['Diagnostic fold 0', 'Diagnostic fold 1', 'Diagnostic fold 2']
               and {r['metric'] for r in pr} == {'AUPRC', 'AUPRC_positive', 'AUPRC_negative'}) else 'FAIL',
    {'rows': len(pr), 'columns': list(pr[0]), 'labels': sorted({r['label'] for r in pr})})

for name, internal in (('diagnostic_tenv_pr.csv', 'diagnostic_pr.csv'),
                       ('diagnostic_tenv_pr.png', 'diagnostic_pr.png')):
    same(f'step-root {name} byte-identical to hvg-metrics-diagnostic-v1/{internal}',
         sha(OUT / 'hvg-metrics-diagnostic-v1' / internal), sha(STEP / name))

v1m = {r['symbol']: r for r in rows(OUT / 'shared-v1/tahoe_mapping.csv')}
v2m = {r['symbol']: r for r in rows(OUT / 'shared-v2/tahoe_mapping.csv')}
symbol_delta = sum(1 for k in v1m if v1m[k]['native_symbol'] != v2m[k]['native_symbol'])
token_delta = sum(1 for k in v1m if v1m[k]['native_token_id'] != v2m[k]['native_token_id'])
add(structure, 'shared-v1 vs shared-v2 Tahoe mapping: same tokens, annotated symbol vintage differs',
    'PASS' if symbol_delta == 229 and token_delta == 0 else 'FAIL',
    {'rows': len(v1m), 'native_symbol_differences': symbol_delta,
     'native_token_id_differences': token_delta,
     'native_column_redefined': 'v1 = deleted 62710-gene matrix column; v2 = projected 18164 axis column',
     'example': {'symbol': 'ACKR5', 'v1_native_symbol': v1m['ACKR5']['native_symbol'],
                 'v2_native_symbol': v2m['ACKR5']['native_symbol'],
                 'token_id': v2m['ACKR5']['native_token_id']}})
kv1 = rows(OUT / 'shared-v1/tenv_mapping.csv')
kv2 = rows(OUT / 'shared-v2/tenv_mapping.csv')
add(structure, 'shared-v1 vs shared-v2 tenv mapping: byte-identical evidence',
    'PASS' if sha(OUT / 'shared-v1/tenv_mapping.csv') == sha(OUT / 'shared-v2/tenv_mapping.csv') else 'FAIL',
    {'rows': len(kv1)})

payload = {'schema': 'bases-155-metadata-binding-audit-v1',
           'scope': 'manifests, non-expression products and existing build records under '
                    '/home/chang/arc/stage2/01_bases; expression payloads (*.h5ad,*.npy,*.npz,*.bin) never opened',
           'expression_payloads_opened': 0,
           'manifest_sha256': manifests,
           'checks': checks, 'structure': structure}
failed = [c for c in checks + structure if c['status'] not in ('PASS', 'EXPECTED_GAP')]
expected_gaps = [c['name'] for c in checks + structure if c['status'] == 'EXPECTED_GAP']
payload['summary'] = {'checks': len(checks), 'structure': len(structure),
                      'non_pass': [c['name'] for c in failed],
                      'expected_gaps': expected_gaps}
if _args.write:
    _args.write.parent.mkdir(parents=True, exist_ok=True)
    _args.write.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
print(json.dumps(payload['summary'], indent=2))
for c in failed:
    print('NON-PASS:', c['name'], json.dumps(c['detail'])[:800])
