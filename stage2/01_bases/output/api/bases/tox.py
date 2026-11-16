"""Unified dataset-level TOX candidates from independent drug admission (ADR 0072).

Candidates live outside the released disease tree. Every independent dataset
supplies one collapsed DEG table and one rank, regardless of mechanism count.
"""
from pathlib import Path
import gzip
import itertools
import tempfile

import numpy as np
import pandas as pd

from . import cardiotox, deg_rank, disease_deg, tox_admission
from .release import ReleaseError, check_manifest, read_json, sha256, write_json

RECIPE = ('Pool independently admitted small-molecule drugs vs matched vehicle '
          'within original experimental strata; preserve source-specific tests '
          'and donor handling; collapse once per independent dataset with '
          'maximum effect and pass-any-contrast padj<0.05 & abs(log2FC)>=2; '
          'one DEG-only rank per dataset, rank_norm=rank/n_deg; finish each '
          'declared qualifying tier without a DEG-count stopping condition.')


def source_plan(ledger):
    """Resolve the admitted independent experiments before reading measurements."""
    plan, experiments, source_ids = [], set(), set()
    drugs = pd.DataFrame(ledger['drugs'])
    for spec in ledger['sources']:
        if spec['disposition'] != 'admitted':
            continue
        source_id = spec['source']
        experiment = spec.get('same_experiment', spec.get('accession', source_id))
        if experiment in experiments or source_id in source_ids:
            raise ReleaseError('duplicate independent TOX experiment: ' + source_id)
        experiments.add(experiment)
        source_ids.add(source_id)
        selected = drugs[drugs.source.eq(source_id) & drugs.disposition.eq('admitted')].copy()
        if selected.empty:
            raise ReleaseError('admitted TOX source has no positive drug: ' + source_id)
        if not (selected.small_molecule.astype(bool).all()
                and selected.evidence_label.eq('positive').all()
                and selected.parent_inchikey.astype(str).str.len().gt(0).all()):
            raise ReleaseError('TOX admission contradicts independent positive identity')
        if selected.source_id.duplicated().any():
            raise ReleaseError('TOX source repeats a drug code')
        plan.append(dict(spec, drugs=selected.reset_index(drop=True)))
    if not plan:
        raise ReleaseError('no admitted TOX experiments')
    disease_deg.tier_ladder([dict(dataset_id=s['source'], tier=s['tier']) for s in plan])
    return sorted(plan, key=lambda s: (s['tier'], s['source']))


def _gse243674_membership(repo, members, identity):
    rows = []
    columns = pd.read_csv(repo / cardiotox.GSE243674_COUNTS, sep=' ', nrows=0).columns[1:]
    panel = cardiotox._well_meta(columns)
    for time_h in sorted(panel.time_h.unique()):
        part = panel[panel.time_h.eq(time_h)]
        controls = part[part.token.eq('VEH')]
        treated = part[part.token.isin(members)]
        if treated.empty:
            continue
        if controls.empty or not set(treated.donor) <= set(controls.donor):
            raise ReleaseError('GSE243674 lacks matched donor/time vehicle')
        for row in part.itertuples(index=False):
            role = 'treated' if row.token in members else 'vehicle' if row.token == 'VEH' else 'excluded'
            rows.append(dict(sample=row.well, source_id=row.token,
                             parent_inchikey=identity.get(row.token, ''),
                             role=role, donor=str(row.donor), time_h=int(time_h),
                             stratum=f'{time_h}h', matched_controls=sorted(controls.well),
                             reason='admitted exposure' if role == 'treated' else
                             'same-time vehicle' if role == 'vehicle' else 'drug not independently admitted'))
    return rows


def _gse244740_membership(repo, members, identity):
    rows = []
    panel = cardiotox.gse244740_panel(repo / cardiotox.GSE244740_SERIES)
    for plate in sorted(set(panel.loc[panel.compound.isin(members), 'plate'])):
        part = panel[panel.plate.eq(plate)]
        controls = sorted(part.loc[part.compound.eq('DMSO'), 'well'])
        if not controls:
            raise ReleaseError('GSE244740 admitted exposure lacks same-plate DMSO')
        for row in part.itertuples(index=False):
            role = 'treated' if row.compound in members else 'vehicle' if row.compound == 'DMSO' else 'excluded'
            rows.append(dict(sample=row.well, source_id=row.compound,
                             parent_inchikey=identity.get(row.compound, ''),
                             role=role, plate=plate, stratum=f'plate{plate}',
                             time=row.time, dose_index=row.dose_index,
                             matched_controls=controls,
                             reason='admitted exposure' if role == 'treated' else
                             'same-plate DMSO' if role == 'vehicle' else 'drug not independently admitted'))
    return rows


def _dtoxs_ipsc_membership(repo, members, identity):
    rows = []
    strata = cardiotox.dtox_pooled_strata(repo / cardiotox.DTOXS_DEG_TAR, members)
    for (cell, plate, hour), record in strata.items():
        stratum = f'{cell}|plate{plate}|{hour}h'
        for code, wells in [('CTRL', record['controls'])] + list(record['treated'].items()):
            for well in wells:
                rows.append(dict(sample=well, source_id=code,
                                 parent_inchikey=identity.get(code, ''),
                                 role='vehicle' if code == 'CTRL' else 'treated',
                                 stratum=stratum, cell_line=cell, plate=plate, time_h=hour,
                                 matched_controls=record['controls'],
                                 reason='source-pipeline QC-retained well'))
    selected = {row['sample'] for row in rows}
    with gzip.open(repo / cardiotox.DTOXS_COUNTS, 'rt') as handle:
        all_wells = handle.readline().rstrip('\n').split('\t')
    for well in all_wells:
        if well not in selected:
            code = well.split('.')[0]
            rows.append(dict(sample=well, source_id=code, role='excluded',
                             parent_inchikey=identity.get(code, ''),
                             stratum='', matched_controls=[],
                             reason='not independently admitted or absent from source QC-retained matched strata'))
    return rows


def legacy_membership(repo, source, drugs):
    """Audit pooled source wells with the same adapter used for measurements."""
    adapter = _LEGACY_ADAPTERS.get(source)
    if adapter is None:
        raise ReleaseError('not an original TOX source: ' + source)
    members = set(drugs.source_id.astype(str))
    identity = drugs.set_index('source_id').parent_inchikey.to_dict()
    rows = adapter[1](Path(repo), members, identity)
    observed = {r['source_id'] for r in rows if r['role'] == 'treated'}
    missing = sorted(members - observed)
    if missing:
        raise ReleaseError(f'{source}: admitted drugs lack source QC-retained matched exposures: {missing}')
    return rows


def validate_membership(audit):
    """Reject duplicate wells or non-vehicle controls before any candidate write."""
    included = [r for r in audit if r['role'] in ('treated', 'vehicle')]
    by_sample = {}
    for row in included:
        if row['sample'] in by_sample:
            raise ReleaseError('TOX source pools a sample twice')
        by_sample[row['sample']] = row
    if not any(row['role'] == 'treated' for row in included):
        raise ReleaseError('TOX source has no admitted treated sample')
    for row in included:
        if row['role'] != 'treated':
            continue
        controls = row.get('matched_controls', [])
        if not controls:
            raise ReleaseError('TOX treated sample lacks matched vehicle')
        for control in controls:
            other = by_sample.get(control)
            if other is None or other['role'] != 'vehicle' or other['stratum'] != row['stratum']:
                raise ReleaseError('TOX control is not the matched-stratum vehicle')


def assemble_candidate(output, sources, axis, inputs, axis_digest, resolver=None,
                       admission=None, environment=None):
    """Publish supplied source-measurement contrasts through the product seam."""
    output = Path(output)
    root = output / 'tox'
    root.mkdir(parents=True, exist_ok=True)
    axis_size = len(axis)
    resolver = resolver or deg_rank.SymbolResolver(axis)
    records, prepared, measured, outputs, seen = [], [], {}, {}, set()
    for source in sources:
        source_id = source['source']
        if source_id in seen:
            raise ReleaseError('duplicate independent TOX dataset: ' + source_id)
        seen.add(source_id)
        validate_membership(source['audit'])
        table, rank, stats = cardiotox.class_deg(source['contrasts'], resolver,
                                                axis_size=axis_size)
        derived = disease_deg.source_deg_rank(table, axis_size=axis_size)
        full_rank = np.asarray(derived['rank'], dtype=np.int64)
        measured[source_id] = dict(deg=np.flatnonzero(full_rank <= derived['n_deg']),
                                   n_deg=derived['n_deg'], informative=derived['informative'])
        prepared.append(dict(dataset_id=source_id, tier=source['tier'],
                             rank=full_rank, n_deg=derived['n_deg'],
                             informative=derived['informative'],
                             informativeness_reason=derived['reason']))
        deg_name, rank_name, audit_name = (source_id + suffix for suffix in
                                          ('.deg.csv', '.rank.csv', '.audit.json'))
        cardiotox.write_deg_table(root / deg_name, table, axis)
        cardiotox.write_rank_table(root / rank_name, rank)
        write_json(root / audit_name, source['audit'])
        for name in (deg_name, rank_name, audit_name):
            outputs[name] = {'sha256': sha256(root / name)}
        contrast_records = source['contrast_records']
        n_case, n_case_basis = cardiotox._arm_counts(contrast_records, 'n_case')
        n_control, n_control_basis = cardiotox._arm_counts(contrast_records, 'n_control')
        records.append(dict(dataset_id=source_id, category='TOX',
                            source_item=source['item'], tier=source['tier'],
                            evidence_type='bulk', species='Homo sapiens',
                            n_case=n_case, n_control=n_control,
                            n_case_basis=n_case_basis, n_control_basis=n_control_basis,
                            contrasts=contrast_records, n_contrasts=len(contrast_records),
                            n_measured=int(stats['n_measured']),
                            n_significant=int(stats['n_significant']),
                            informative=bool(stats['informative']),
                            informativeness=dict(reason=stats['informativeness_reason'],
                                                  n_distinct_effect=stats['n_distinct_effect']),
                            alignment=stats['alignment'], method=source['method'],
                            member_drugs=sorted(source['member_drugs']),
                            member_source_ids=sorted(source['member_source_ids']),
                            drug_identities=source.get('drug_identities', []),
                            n_small_molecule_drugs=len(source['member_drugs']),
                            deg=deg_name, rank=rank_name, audit=audit_name,
                            inputs=source.get('inputs', {})))
    if not records:
        raise ReleaseError('TOX candidate contains no independently admitted dataset')
    fusion = disease_deg.fuse_disease(prepared, axis_size=axis_size)
    ledger = disease_deg.tier_ledger('TOX', records, measured)
    disease_deg.deg1000_rows(fusion).to_csv(root / 'tox.deg1000.csv', index=False,
                                           float_format='%.17g', lineterminator='\n')
    outputs['tox.deg1000.csv'] = {'sha256': sha256(root / 'tox.deg1000.csv')}
    per_source = {row['dataset_id']: row for row in fusion['sources']}
    catalog = disease_deg.build_catalog(
        'TOX', 'tox', min(record['tier'] for record in records), 'bulk', axis,
        axis_digest, records, axis_size=axis_size, n_deg=fusion['n_deg'],
        deg1000=fusion['deg1000'], tiering=ledger['tiering'],
        tier_sources=ledger['tier_sources'], source_deg=per_source)
    # The shared builder's historical path predates the tables/ directory.
    catalog['axis']['path'] = cardiotox.AXIS
    catalog_sha = disease_deg.write_catalog(root / 'catalog.json', catalog)
    outputs['catalog.json'] = {'sha256': sha256(root / 'catalog.json')}
    if admission is not None:
        write_json(root / 'admission.json', admission)
        outputs['admission.json'] = {'sha256': sha256(root / 'admission.json')}
    manifest = disease_deg.assemble_manifest(
        'TOX', RECIPE, inputs, outputs,
        {'datasets': len(records), 'n_deg': fusion['n_deg'], 'deg1000': fusion['deg1000']},
        catalog_sha)
    manifest['axis'] = dict(path=cardiotox.AXIS, sha256=axis_digest, size=axis_size)
    if environment is not None:
        manifest['environment'] = environment
    write_json(root / 'manifest.json', manifest)
    check_manifest(root, manifest)
    return catalog


def _gse243674_contrasts(repo, members, workdir):
    frame, contrasts, metrics = cardiotox.gse243674_class_contrasts(repo, members, 'TOX', workdir)
    method = dict(name='donor-paired DESeq2', model='~ donor + group',
                  test='Wald; unshrunken log2FC; BH and independent filtering',
                  pooling='all admitted drugs per original time, same-time vehicle')
    inputs = [cardiotox.GSE243674_COUNTS, cardiotox.HGNC, cardiotox.DESEQ2_SCRIPT]
    return frame, contrasts, metrics, method, inputs


def _dtoxs_ipsc_contrasts(repo, members, workdir):
    frame, contrasts, metrics = cardiotox.dtoxs_class_contrasts(
        repo, cardiotox.read_dtoxs_counts(repo), members, 'TOX', workdir)
    method = dict(name='deposition pinned edgeR',
                  test='TMM; trended dispersion; exactTest; remove_small; BH',
                  pooling='all admitted QC-retained drugs within cell-line/plate/time strata')
    inputs = [cardiotox.DTOXS_COUNTS, cardiotox.DTOXS_DEG_TAR,
              cardiotox.DTOXS_PROGRAMS, cardiotox.DTOXS_DESIGN, cardiotox.EDGE_R_SCRIPT]
    return frame, contrasts, metrics, method, inputs


def _gse244740_contrasts(repo, members, workdir):
    frame, contrasts, metrics = cardiotox.gse244740_class_contrasts(repo, members, 'TOX')
    method = dict(name='ScreenSeq pooled Welch',
                  test='Welch on log2(CPM+1); source-full library; BH per contrast',
                  pooling='all admitted drug doses within plate vs same-plate DMSO')
    inputs = [cardiotox.GSE244740_COUNTS, cardiotox.GSE244740_SERIES]
    return frame, contrasts, metrics, method, inputs


_LEGACY_ADAPTERS = {
    'gse243674': (_gse243674_contrasts, _gse243674_membership),
    'dtoxs_ipsc': (_dtoxs_ipsc_contrasts, _dtoxs_ipsc_membership),
    'gse244740': (_gse244740_contrasts, _gse244740_membership),
}


def source_contrasts(repo, spec, workdir):
    """Compute a dataset once from the independently admitted member union."""
    repo, workdir = Path(repo), Path(workdir)
    source, drugs = spec['source'], spec['drugs']
    members = drugs.source_id.astype(str).tolist()
    adapter = _LEGACY_ADAPTERS.get(source)
    if adapter is None:
        from . import tox_sources
        frame, contrasts, metadata = tox_sources.source_contrasts(repo, source, drugs, workdir)
        return dict(source=source, tier=spec['tier'], item=spec['item'],
                    contrasts=frame, contrast_records=contrasts,
                    member_drugs=drugs.drug_name.astype(str).tolist(),
                    member_source_ids=members,
                    audit=metadata['audit'], method=metadata['method'],
                    inputs=metadata['inputs']), metadata.get('metrics', {})
    frame, contrasts, metrics, method, inputs = adapter[0](repo, members, workdir)
    audit = legacy_membership(repo, source, drugs)
    return dict(source=source, tier=spec['tier'], item=spec['item'],
                contrasts=frame, contrast_records=contrasts,
                member_drugs=drugs.drug_name.astype(str).tolist(),
                member_source_ids=members, audit=audit, method=method,
                inputs=inputs), metrics


def build_candidate(repo, output, workdir=None):
    """Recompute all independently admitted TOX sources across every declared tier."""
    repo, output = Path(repo), Path(output)
    workdir = Path(workdir) if workdir else output.parent / (output.name + '-work')
    workdir.mkdir(parents=True, exist_ok=True)
    admission_path = workdir / 'admission'
    tox_admission.build_ledger(repo, admission_path)
    admission = tox_admission.ToxAdmission(admission_path, repo).value
    plan = source_plan(admission)
    axis, resolver = cardiotox.load_resolver(repo)
    inputs = dict(admission['inputs'])
    implementation = ['tox.py', 'tox_sources.py', 'bulk_deg.py', 'cardiotox.py', 'deg_rank.py',
                      'disease_deg.py', 'disease_deg.R', 'r_env.py', 'release.py',
                      'tox_admission.py', 'cardiotox_deseq2.R', 'cardiotox_edgeR.R']
    for name in implementation:
        cardiotox.add_input(inputs, repo / 'stage2/01_bases/output/api/bases' / name, repo)
    for path in (cardiotox.AXIS, cardiotox.CROSSWALK, cardiotox.HGNC):
        cardiotox.add_input(inputs, repo / path, repo)
    records, metrics, union = [], {}, set()
    entered_tiers = []
    for tier, specs in itertools.groupby(plan, key=lambda s: s['tier']):
        entered_tiers.append(tier)
        for spec in specs:
            print('TOX source start: ' + spec['source'], flush=True)
            record, runtime = source_contrasts(repo, spec, workdir / spec['source'])
            validate_membership(record['audit'])
            record['drug_identities'] = spec['drugs'][['source_id', 'drug_name', 'parent_inchikey']].to_dict('records')
            # Actual measured source paths, not former class DEG artifacts.
            pinned = {}
            paths = record['inputs']
            for path in paths if isinstance(paths, list) else paths.keys():
                cardiotox.add_input(pinned, repo / path, repo)
            record['inputs'] = pinned
            inputs.update(pinned)
            table, _, stats = cardiotox.class_deg(record['contrasts'], resolver)
            if stats['informative']:
                union.update(int(i) for i in table.loc[table.deg, 'gene_index'])
            records.append(record)
            metrics[spec['source']] = runtime
            print(f"TOX source done: {spec['source']} n_deg={stats['n_deg']}", flush=True)
    catalog = assemble_candidate(output, records, axis, inputs, sha256(repo / cardiotox.AXIS),
                                 resolver=resolver, admission=admission,
                                 environment=dict(cardiotox.python_environment(),
                                                  r=cardiotox.r_environment()))
    # Runtime peaks never enter deterministic scientific payload digests.
    write_json(workdir / 'run_metrics.json', dict(sources=metrics, tiers=entered_tiers,
                                                 n_deg=catalog['n_deg']))
    return catalog


def verify_candidate(repo, output, workdir=None):
    """Replay every entered dataset from pinned measurements, without cache reuse."""
    output = Path(output)
    expected = ToxCandidate(output).digests()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='tox-replay-', dir=output.parent) as temporary:
        replay = Path(temporary) / 'candidate'
        build_candidate(repo, replay, workdir=workdir or Path(temporary) / 'work')
        if ToxCandidate(replay).digests() != expected:
            raise ReleaseError('TOX candidate source replay differs')
        if sha256(replay / 'tox/manifest.json') != sha256(output / 'tox/manifest.json'):
            raise ReleaseError('TOX candidate manifest replay differs')
    return expected


class ToxCandidate:
    """Verify and inspect candidate dataset membership without opening counts."""
    def __init__(self, root):
        self.root = Path(root) / 'tox'
        self.manifest = read_json(self.root / 'manifest.json')
        check_manifest(self.root, self.manifest)
        self.catalog = disease_deg.read_catalog(self.root / 'catalog.json')

    def audit(self, source):
        entry = next(r for r in self.catalog['datasets'] if r['dataset_id'] == source)
        name = entry['audit']
        if sha256(self.root / name) != self.manifest['outputs'][name]['sha256']:
            raise ReleaseError('TOX membership audit digest changed')
        return read_json(self.root / entry['audit'])

    def deg1000(self):
        path = self.root / 'tox.deg1000.csv'
        if sha256(path) != self.manifest['outputs'][path.name]['sha256']:
            raise ReleaseError('TOX DEG1000 digest changed')
        return pd.read_csv(path, float_precision='round_trip')

    def digests(self):
        return {name: record['sha256'] for name, record in self.manifest['outputs'].items()}
