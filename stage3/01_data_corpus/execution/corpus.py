#!/usr/bin/env python3
"""Offline frozen-evidence phenotype corpus build and verification."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
import statistics as standard_statistics
from pathlib import Path
import subprocess
import sys
import tempfile

import pyarrow as pa
import pyarrow.parquet as pq
import scipy
from scipy import stats as scipy_stats

CORE_CONDITION = ('condition_cell_type', 'condition_drug_time', 'condition_hypoxia_time',
                  'condition_reoxygen_time', 'condition_drug_schedule')
RULES = {'schema_version': 1, 'rule_version': 'fc-bh-v1', 'fc_low': .8, 'fc_high': 1.2,
         'adjusted_p_cutoff': .05, 'dose_unit': 'mol/L', 'time_unit': 'h',
         'family_fields': ['source', 'condition', 'phenotype', 'assay_method', 'batch'],
         'interpretation': 'Operational change rule; zero is not equivalence; no final-subset FDR guarantee.'}
ANNOTATION = pa.struct([
    pa.field('measurement_id', pa.string(), nullable=False),
    pa.field('comparison_id', pa.string(), nullable=False),
    pa.field('assay_method', pa.string(), nullable=False),
    pa.field('unit', pa.string(), nullable=False),
    pa.field('representation', pa.string(), nullable=False),
    pa.field('asset', pa.string(), nullable=False),
    pa.field('selector', pa.string(), nullable=False),
    pa.field('experimental_unit', pa.string()),
    pa.field('normalized_value', pa.float64()),
    pa.field('source_metadata', pa.string(), nullable=False),
])
SCHEMA = pa.schema([
    pa.field('sample_id', pa.string(), nullable=False),
    pa.field('condition_id', pa.string(), nullable=False),
    pa.field('condition', pa.string(), nullable=False),
    pa.field('condition_cell_type', pa.string(), nullable=False),
    pa.field('condition_drug_time', pa.float64()),
    pa.field('condition_hypoxia_time', pa.float64()),
    pa.field('condition_reoxygen_time', pa.float64()),
    pa.field('condition_drug_schedule', pa.string(), nullable=False),
    pa.field('phenotype', pa.string(), nullable=False),
    pa.field('source', pa.string(), nullable=False),
    pa.field('drug', pa.string(), nullable=False),
    pa.field('dose', pa.float64(), nullable=False),
    pa.field('dose_unit', pa.string(), nullable=False),
    pa.field('endpoint', pa.int8(), nullable=False),
    pa.field('readout', pa.list_(pa.field('element', pa.float64(), nullable=False)), nullable=False),
    pa.field('annotations', pa.list_(pa.field('element', ANNOTATION, nullable=False)), nullable=False),
])


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False, separators=(',', ':'))


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def file_hash(path):
    hasher = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            hasher.update(block)
    return hasher.hexdigest()


def root_path():
    for directory in (Path(__file__).resolve().parent, *Path(__file__).resolve().parents):
        if (directory / '.git').exists():
            return directory
    raise ValueError('Repository root with .git not found')


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def probability(value):
    return finite(value) and 0 <= value <= 1


def condition_record(record):
    condition = dict(record['condition'])
    if any(key not in condition for key in CORE_CONDITION):
        raise ValueError('Missing required condition column')
    if any(not key.startswith('condition_') for key in condition):
        raise ValueError('Condition keys must begin condition_')
    for key in CORE_CONDITION[1:4]:
        if condition[key] is not None and (not finite(condition[key]) or condition[key] < 0):
            raise ValueError('Invalid condition duration')
    for key in CORE_CONDITION[1:4]:
        if condition[key] is not None:
            condition[key] = float(condition[key])
    if not isinstance(condition['condition_cell_type'], str) or not condition['condition_cell_type']:
        raise ValueError('Invalid cell type')
    if not isinstance(condition['condition_drug_schedule'], str) or not condition['condition_drug_schedule']:
        raise ValueError('Invalid drug schedule')
    if any(value is None or value == 'unknown' for value in condition.values()):
        condition.setdefault('condition_unresolved_context', record['source'])
    canonical(condition)
    return condition


def validate_bundle(bundle, test_only):
    if bundle.get('schema_version') != 1 or not isinstance(bundle.get('sources'), list) or not isinstance(bundle.get('records'), list):
        raise ValueError('Invalid frozen-input bundle')
    sources = {}
    assets = {}
    root = root_path()
    registered = {}
    if not test_only:
        registry = json.loads((root/'pheno_data/registry.json').read_text())
        registered = {entry['dataset_id']:entry for entry in registry['items']}
    for source in bundle['sources']:
        name = source['source']
        if name in sources or not isinstance(name, str) or not name:
            raise ValueError('Duplicate or invalid source')
        if bool(source.get('test_only')) != bool(test_only):
            raise ValueError('Production/test-only source boundary violated')
        if not test_only:
            entry = registered.get(name)
            if not entry or source.get('item') != 'pheno_data/' + entry['item']:
                raise ValueError('Unregistered production dataset')
            if source.get('status') in ('candidate','unavailable'):
                if entry['status'] not in ('candidate','unavailable'):
                    raise ValueError('Active source cannot be silently marked unavailable')
            elif entry['status'] != 'active':
                raise ValueError('Production records require active archived source')
        sources[name] = source
        if source.get('status') in ('candidate','unavailable'):
            if source.get('assets') or not source.get('reason'):
                raise ValueError('Unavailable source must have reason and no verified assets')
            continue
        if not source.get('assets'):
            raise ValueError('Active source has no frozen assets')
        if not test_only:
            item_dir = root/source['item']
            provenance = json.loads((item_dir/'author_log.json').read_text())
            if provenance.get('dataset_id') != name:
                raise ValueError('Source provenance identity mismatch')
            manifest = item_dir/'manifest.sha256'
            expected = {}
            for line in manifest.read_text().splitlines():
                sha, relative = line.split('  ', 1)
                relative = relative.removeprefix('source/')
                expected[str(Path(source['item'])/'source'/relative)] = sha
            actual_assets = {asset['path']:asset['sha256'] for asset in source['assets']}
            if not expected or any(actual_assets.get(path) != sha for path,sha in expected.items()):
                raise ValueError('Source registry omits or changes original manifest assets')
            local_paths = {str(Path(source['item'])/entry['path']) for entry in provenance.get('local_additions',[])}
            if any(path not in expected and path not in local_paths for path in actual_assets):
                raise ValueError('Undeclared local curation asset')
        for asset in source['assets']:
            path = Path(asset['path'])
            absolute = path if path.is_absolute() else root / path
            if not test_only:
                item = Path(source['item'])
                if path.is_absolute() or '..' in path.parts or item.parts[0] != 'pheno_data' or not path.is_relative_to(item):
                    raise ValueError('Production source asset outside registered phenotype item')
            if not absolute.is_file() or file_hash(absolute) != asset['sha256']:
                raise ValueError(f'Frozen asset missing/hash mismatch: {path}')
            if asset['path'] in assets and assets[asset['path']] != asset['sha256']:
                raise ValueError('Conflicting frozen asset identities')
            assets[asset['path']] = asset['sha256']
    seen_comparisons, seen_measurements = set(), set()
    for record in bundle['records']:
        if record['source'] not in sources or sources[record['source']].get('status') in ('candidate','unavailable'):
            raise ValueError('Record references unknown/unavailable source')
        comparison_id = record['comparison_id']
        if not comparison_id or comparison_id in seen_comparisons:
            raise ValueError('Duplicate/empty comparison identity')
        seen_comparisons.add(comparison_id)
        condition_record(record)
        for field in ('phenotype', 'assay_method'):
            if not isinstance(record[field], str) or not record[field]:
                raise ValueError(f'Invalid {field}')
        if record['assay_method'] in ('MTT', 'CCK8') and record['phenotype'] != 'cell_viability':
            raise ValueError('MTT/CCK8 must map to cell_viability')
        if record['assay_method'] in ('CTG', 'Calcein') and record['phenotype'] == 'cell_viability':
            raise ValueError('Unsupported assay to cell_viability mapping')
        if not isinstance(record['measurements'], list):
            raise ValueError('Measurements must be a list')
        own_assets = {asset['path'] for asset in sources[record['source']]['assets']}
        for measurement in record['measurements']:
            identity = (record['source'], measurement['measurement_id'])
            if not measurement['measurement_id'] or identity in seen_measurements:
                raise ValueError('Duplicate measurement identity/representation')
            seen_measurements.add(identity)
            if measurement['asset'] not in own_assets or not measurement['selector']:
                raise ValueError('Broken measurement source evidence')
            if measurement['representation'] not in ('raw', 'normalized_only', 'summary'):
                raise ValueError('Invalid measurement representation')
            if measurement.get('value') is not None and not finite(measurement['value']):
                raise ValueError('Nonfinite measurement')
            if measurement.get('normalized_value') is not None and not finite(measurement['normalized_value']):
                raise ValueError('Nonfinite normalized measurement')
            if measurement['assay_method'] != record['assay_method'] or not measurement['unit']:
                raise ValueError('Measurement assay/unit mismatch')
    return sources


def family_identity(record):
    key = [record['source'], condition_record(record), record['phenotype'], record['assay_method'], record.get('batch')]
    return 'family-' + digest(key)[:20]


def recompute_test(stats):
    """Summarize technical values by independently documented experimental unit."""
    if stats.get('independent_units') is not True:
        return None
    design = stats.get('design')
    if design == 'independent_summary':
        treated, control = stats.get('treated_summary',{}), stats.get('control_summary',{})
        for summary in (treated,control):
            if (not finite(summary.get('mean')) or not finite(summary.get('sd')) or summary['sd'] < 0 or
                    not isinstance(summary.get('n'),int) or isinstance(summary['n'],bool) or summary['n'] < 2):
                return None
        if treated['sd'] == 0 and control['sd'] == 0:
            return None
        outcome = scipy_stats.ttest_ind_from_stats(treated['mean'],treated['sd'],treated['n'],
                    control['mean'],control['sd'],control['n'],equal_var=False,alternative='two-sided')
        return float(outcome.pvalue), 'recomputed_Welch_summary', treated['mean'], control['mean']
    if design not in ('independent','paired'):
        return None
    units = []
    for side in ('treated','control'):
        values,identities = stats.get(side+'_values'),stats.get(side+'_unit_ids')
        if (not isinstance(values,list) or not isinstance(identities,list) or len(values)!=len(identities) or
                not values or not all(finite(value) for value in values) or
                not all(isinstance(identity,str) and identity for identity in identities)):
            return None
        grouped=defaultdict(list)
        for identity,value in zip(identities,values):
            grouped[identity].append(value)
        means={identity:standard_statistics.mean(values) for identity,values in grouped.items()}
        if len(means)<2:
            return None
        units.append(means)
    treated,control=units
    if design=='paired':
        if set(treated)!=set(control):
            return None
        identities=sorted(treated)
        tvalues=[treated[identity] for identity in identities]
        cvalues=[control[identity] for identity in identities]
        if standard_statistics.variance([a-b for a,b in zip(tvalues,cvalues)])==0:
            return None
        outcome=scipy_stats.ttest_rel(tvalues,cvalues,alternative='two-sided')
        origin='recomputed_paired'
    else:
        if set(treated)&set(control):
            return None
        tvalues,cvalues=list(treated.values()),list(control.values())
        if standard_statistics.variance(tvalues)==0 and standard_statistics.variance(cvalues)==0:
            return None
        outcome=scipy_stats.ttest_ind(tvalues,cvalues,equal_var=False,alternative='two-sided')
        origin='recomputed_Welch'
    return float(outcome.pvalue),origin,standard_statistics.mean(tvalues),standard_statistics.mean(cvalues)


def test_result(record):
    stats = record.get('statistics', {})
    result = {'p_raw': None, 'p_adjusted': None, 'statistical_origin': 'unavailable',
              'test_missing_reason': 'no_applicable_numerical_test'}
    known = (stats.get('sidedness') == 'two-sided' and bool(stats.get('method')) and
             bool(stats.get('evidence')) and bool(record.get('batch')))
    if not known:
        return result
    if (probability(stats.get('p_adjusted')) and stats.get('correction') == 'BH' and
            stats.get('family_matches') is True and stats.get('family_complete') is True):
        return dict(result, p_adjusted=float(stats['p_adjusted']), statistical_origin='source_BH', test_missing_reason=None)
    if probability(stats.get('p_raw')) and stats.get('correction') in ('none','raw'):
        return dict(result, p_raw=float(stats['p_raw']), statistical_origin='source_raw_p',
                    test_missing_reason='family_adjustment_unavailable')
    recomputed = recompute_test(stats)
    if recomputed is not None and probability(recomputed[0]):
        return dict(result,p_raw=recomputed[0],statistical_origin=recomputed[1],
                    tested_treated_mean=recomputed[2],tested_control_mean=recomputed[3],
                    test_missing_reason='family_adjustment_unavailable')
    return result


def adjusted_results(bundle):
    results = {record['comparison_id']: test_result(record) for record in bundle['records']}
    by_family = defaultdict(list)
    for record in bundle['records']:
        by_family[family_identity(record)].append(record)
    for records in by_family.values():
        actual_members = sorted(record['comparison_id'] for record in records)
        for record in records:
            result = results[record['comparison_id']]
            declared = record.get('statistics',{}).get('family_members')
            if declared is not None and (not isinstance(declared,list) or sorted(declared) != actual_members):
                raise ValueError('Declared BH family membership differs from frozen comparisons')
            if declared is None and result['p_adjusted'] is not None:
                result.update(p_adjusted=None, statistical_origin='unavailable', test_missing_reason='undeclared_BH_family')
        eligible = [record for record in records if probability(results[record['comparison_id']]['p_raw'])]
        complete = all(record.get('statistics',{}).get('family_complete') is True and
                       record.get('statistics',{}).get('family_matches') is True and
                       record.get('statistics',{}).get('family_members') == actual_members for record in eligible)
        if eligible and complete:
            adjusted = scipy_stats.false_discovery_control([results[record['comparison_id']]['p_raw'] for record in eligible], method='bh')
            for record,p_adj in zip(eligible,adjusted):
                result=results[record['comparison_id']]
                result.update(p_adjusted=float(p_adj), test_missing_reason=None)
    return results


def coverage_summary(comparisons, rows):
    """Count evidence and admission at one declared source/condition scope."""
    measurements = [(record,measurement) for record in comparisons for measurement in record['measurements']]
    units = {(record['source'],record.get('batch'),measurement['experimental_unit'])
             for record,measurement in measurements if measurement.get('experimental_unit')}
    unknown_units = sum(not measurement.get('experimental_unit') for _,measurement in measurements)
    independent = [(record,measurement) for record,measurement in measurements
                   if measurement.get('experimental_unit') and
                   (measurement.get('experimental_unit_independence') == 'independent' or
                    record.get('statistics',{}).get('independent_units') is True)]
    unknown_independence = len(measurements)-len(independent)
    independent_units = {(r['source'],r.get('batch'),m['experimental_unit']) for r,m in independent}
    if not measurements:
        units_status = 'no_measurements'
    elif unknown_units == len(measurements):
        units_status = 'unknown'
    elif unknown_units:
        units_status = 'partially_known'
    else:
        units_status = 'known_identifiers'
    identity_eligible = lambda record: bool(record.get('drug') and record.get('identity_evidence'))
    dose_eligible = lambda record: finite(record.get('dose')) and record['dose'] > 0 and bool(record.get('dose_evidence'))
    return {
        'assessment_status':'assessed' if comparisons else 'not_run',
        'comparisons':len(comparisons),
        'measurement_entries':len(measurements),
        'finite_source_values':sum(finite(m.get('value')) for _,m in measurements),
        'representation_counts':dict(Counter(m['representation'] for _,m in measurements)),
        'known_unit_count':len(units),
        'unknown_unit_measurements':unknown_units,
        'experimental_units':len(units) if not unknown_units else None,
        'experimental_units_status':units_status,
        'documented_independent_unit_count':len(independent_units),
        'independence_unknown_entries':unknown_independence,
        'independent_experimental_units':len(independent_units) if not unknown_independence else None,
        'independence_status':'not_run' if not measurements else ('unknown_or_incomplete' if unknown_independence else 'documented'),
        'supplied_normalized_representations':sum(m.get('normalized_value') is not None for _,m in measurements),
        'identity_eligible_comparisons':sum(identity_eligible(r) for r in comparisons),
        'dose_eligible_comparisons':sum(dose_eligible(r) for r in comparisons),
        'identity_and_dose_eligible_comparisons':sum(identity_eligible(r) and dose_eligible(r) for r in comparisons),
        'fc_computable_comparisons':sum(r['fc'] is not None for r in comparisons),
        'raw_p_available_comparisons':sum(r['p_raw'] is not None for r in comparisons),
        'adjusted_p_available_comparisons':sum(r['p_adjusted'] is not None for r in comparisons),
        'admitted_comparisons':sum(not r['exclusion_reasons'] for r in comparisons),
        'excluded_comparisons':sum(bool(r['exclusion_reasons']) for r in comparisons),
        'exclusion_reasons_overlapping':dict(Counter(reason for r in comparisons for reason in r['exclusion_reasons'])),
        'exclusion_count_semantics':'One comparison may have multiple reasons; reason counts overlap and must not be summed as mutually exclusive attrition.',
        'sample_rows':len(rows),
        'class_rows':{str(label):sum(row['endpoint']==label for row in rows) for label in (-1,0,1)},
        'drug_dose_union':len({(r['condition_id'],r['phenotype'],r['drug'],r['dose']) for r in rows}),
        'distinct_drugs':len({r['drug'] for r in rows}),
        'count_scope':'Drug-dose union preserves complete condition/phenotype, removes source/endpoint; unit counts identify known units and separately disclose unknown independence.'}


def derive(bundle, test_only):
    sources = validate_bundle(bundle, test_only)
    statistics_results = adjusted_results(bundle)
    comparisons, excluded, partitions = [], [], {}
    groups = {}
    families = {}
    for record in sorted(bundle['records'], key=lambda row: (row['source'], row['comparison_id'])):
        condition = condition_record(record)
        condition_id = 'condition-' + digest({'version': 1, 'condition': condition})[:20]
        partition_key = (condition_id, record['phenotype'])
        partitions[partition_key] = {'condition_id': condition_id, 'condition': condition, 'phenotype': record['phenotype']}
        stats = record.get('statistics', {})
        result = statistics_results[record['comparison_id']]
        family_id = family_identity(record)
        family = families.setdefault(family_id, {'family_id': family_id, 'source': record['source'],
                    'condition_id': condition_id, 'phenotype': record['phenotype'], 'assay_method': record['assay_method'],
                    'batch': record.get('batch'), 'members': [], 'eligible_members': [], 'method': 'BH'})
        family['members'].append(record['comparison_id'])
        if result['p_raw'] is not None or result['p_adjusted'] is not None:
            family['eligible_members'].append(record['comparison_id'])
        reasons = list(record.get('quality_reasons', []))
        if not record.get('drug') or not record.get('identity_evidence'):
            reasons.append('unresolved_drug_identity')
        if not finite(record.get('dose')) or record['dose'] <= 0 or not record.get('dose_evidence'):
            reasons.append('unresolved_canonical_dose')
        treated = result.get('tested_treated_mean', stats.get('treated_mean'))
        control = result.get('tested_control_mean', stats.get('control_mean'))
        fc = None
        if (stats.get('scale_compatible') is True and stats.get('control_matched') is True and
                record.get('control_evidence') and finite(treated) and treated >= 0 and finite(control) and control > 0):
            fc = treated / control
            if not finite(fc):
                fc = None
        if fc is None:
            reasons.append('unevaluable_fold_change')
        endpoint = None
        if fc is not None:
            if .8 < fc < 1.2:
                endpoint = 0
            elif result['p_adjusted'] is None:
                reasons.append('missing_required_adjusted_p')
            elif result['p_adjusted'] >= .05:
                endpoint = 0
            else:
                endpoint = 1 if fc >= 1.2 else -1
        valid_values = [m for m in record['measurements'] if finite(m.get('value'))]
        if not valid_values:
            reasons.append('no_valid_readout')
        reasons = sorted(set(reasons))
        comparison = dict(record, condition=condition, condition_id=condition_id, family_id=family_id,
                          fc=fc, endpoint=endpoint if not reasons else None, exclusion_reasons=reasons, **result)
        comparisons.append(comparison)
        if reasons:
            excluded.append({'comparison_id': record['comparison_id'], 'source': record['source'], 'reasons': reasons})
            continue
        key = (condition_id, record['phenotype'], record['source'], record['drug'], float(record['dose']), endpoint)
        if key not in groups:
            groups[key] = {'sample_id': 'sample-' + digest(key)[:24], 'condition_id': condition_id,
                'condition': canonical(condition), **{field: condition[field] for field in CORE_CONDITION},
                'phenotype': record['phenotype'], 'source': record['source'], 'drug': record['drug'],
                'dose': float(record['dose']), 'dose_unit': 'mol/L', 'endpoint': endpoint, 'readout': [], 'annotations': []}
        row = groups[key]
        for measurement in valid_values:
            row['readout'].append(float(measurement['value']))
            row['annotations'].append({field: measurement.get(field) for field in
                ('measurement_id', 'assay_method', 'unit', 'representation', 'asset', 'selector', 'experimental_unit', 'normalized_value')} |
                {'comparison_id': record['comparison_id'], 'source_metadata': canonical(measurement)})
    rows = sorted(groups.values(), key=lambda row: row['sample_id'])
    for partition in partitions.values():
        partition['rows'] = [row for row in rows if row['condition_id'] == partition['condition_id'] and row['phenotype'] == partition['phenotype']]
        partition['path'] = partition['condition_id'] + '-' + digest(partition['phenotype'])[:12] + '.parquet'
        partition['data_status'] = 'numerically_valid_small_corpus' if partition['rows'] else 'empty_after_QC'
    coverage = coverage_summary(comparisons,rows)
    coverage.update(source_count=len(sources),
        source_status={name:source.get('status','archived') for name,source in sources.items()},
        source_status_counts=dict(Counter(source.get('status','archived') for source in sources.values())),
        archived_source_count=sum(source.get('status','archived') in ('archived','active') for source in sources.values()),
        candidate_source_count=sum(source.get('status')=='candidate' for source in sources.values()),
        unavailable_source_count=sum(source.get('status')=='unavailable' for source in sources.values()),
        model_view_500_assessment='not_run',training_status='not_trained',
        scope='Drug-dose union is within complete condition/phenotype; source and endpoint removed.')
    coverage['partition_coverage']=[]
    for partition in sorted(partitions.values(),key=lambda p:p['path']):
        subset=[r for r in comparisons if r['condition_id']==partition['condition_id'] and r['phenotype']==partition['phenotype']]
        coverage['partition_coverage'].append(dict(coverage_summary(subset,partition['rows']),
                condition_id=partition['condition_id'],phenotype=partition['phenotype'],condition=partition['condition']))
    coverage['per_source']={}
    for name,source in sorted(sources.items()):
        subset=[r for r in comparisons if r['source']==name]
        selected=[r for r in rows if r['source']==name]
        summary=dict(coverage_summary(subset,selected),source=name,status=source.get('status','archived'))
        summary['partition_coverage']=[]
        for partition in sorted(partitions.values(),key=lambda p:p['path']):
            partition_comparisons=[r for r in subset if r['condition_id']==partition['condition_id'] and r['phenotype']==partition['phenotype']]
            if not partition_comparisons:
                continue
            partition_rows=[r for r in selected if r['condition_id']==partition['condition_id'] and r['phenotype']==partition['phenotype']]
            summary['partition_coverage'].append(dict(coverage_summary(partition_comparisons,partition_rows),
                  condition_id=partition['condition_id'],phenotype=partition['phenotype']))
        coverage['per_source'][name]=summary
    return {'partitions': sorted(partitions.values(), key=lambda p:p['path']), 'comparisons': comparisons,
            'families': sorted(families.values(), key=lambda f:f['family_id']), 'exclusions': excluded, 'coverage': coverage}


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + '\n')


def payload(bundle, directory, test_only, frozen_rules=None):
    derived = derive(bundle, test_only)
    rules = dict(RULES, software={'python':sys.version.split()[0], 'pyarrow':pa.__version__, 'scipy':scipy.__version__},
                 implementation_sha256=file_hash(__file__),
                 code_commit=subprocess.run(['git','log','-1','--format=%H','--',str(Path(__file__).resolve().relative_to(root_path()))],cwd=root_path(),capture_output=True,text=True,check=True).stdout.strip())
    if frozen_rules is not None:
        if any(frozen_rules.get(key) != value for key,value in RULES.items()):
            raise ValueError('Unsupported frozen rule version or thresholds')
        rules = frozen_rules
    write_json(directory/'input.json', bundle)
    write_json(directory/'rules.json', rules)
    write_json(directory/'schema.json', {'schema_version':1, 'arrow_schema':str(SCHEMA),
        'endpoint':[-1,0,1], 'readout':'nonempty finite float64 list aligned with annotations',
        'key':['condition_id','phenotype','source','drug','dose','endpoint']})
    write_json(directory/'source_registry.json',bundle['sources'])
    write_json(directory/'source_audits.json',bundle.get('source_audits',{}))
    for name in ('comparisons','families','exclusions','coverage'):
        write_json(directory/(name+'.json'), derived[name])
    partitions = []
    for partition in derived['partitions']:
        table = pa.Table.from_pylist(partition['rows'], schema=SCHEMA)
        pq.write_table(table,directory/partition['path'],compression='zstd',version='2.6')
        partitions.append({key:value for key,value in partition.items() if key != 'rows'} |
                          {'sample_rows':len(partition['rows'])})
    write_json(directory/'manifest.json',{'schema_version':1, 'input_sha256':digest(bundle), 'test_only':test_only,
        'rule_version':RULES['rule_version'], 'partitions':partitions, 'data_status':'numerically_valid' if derived['coverage']['sample_rows'] else 'empty_after_QC'})
    hashes = {path.name:file_hash(path) for path in sorted(directory.iterdir())}
    (directory/'checksums.sha256').write_text(''.join(f'{sha}  {name}\n' for name,sha in hashes.items()))
    return hashes


def build(input_path, output, test_only=False):
    bundle = json.loads(Path(input_path).read_text())
    output = Path(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.corpus-',dir=output.parent) as temp:
        directory=Path(temp)
        hashes=payload(bundle,directory,test_only)
        if output.exists():
            existing={path.name:file_hash(path) for path in output.iterdir() if path.is_file() and path.name!='checksums.sha256'}
            if existing!=hashes or (output/'checksums.sha256').read_bytes()!=(directory/'checksums.sha256').read_bytes():
                raise ValueError('Immutable release conflict: choose a new release identity')
        else:
            directory.rename(output)
    return verify(output,test_only)


def verify(release, test_only=False):
    release=Path(release)
    entries={}
    for line in (release/'checksums.sha256').read_text().splitlines():
        sha,name=line.split('  ',1)
        if name in entries or Path(name).name!=name or name=='checksums.sha256':
            raise ValueError('Invalid release checksum inventory')
        entries[name]=sha
    actual={path.name for path in release.iterdir()}
    if actual!=set(entries)|{'checksums.sha256'}:
        raise ValueError('Release file inventory mismatch')
    for name,sha in entries.items():
        if file_hash(release/name)!=sha:
            raise ValueError('Release hash mismatch: '+name)
    manifest=json.loads((release/'manifest.json').read_text())
    if manifest['test_only']!=test_only:
        raise ValueError('Production/test-only release boundary violated')
    seen=set()
    for partition in manifest['partitions']:
        table=pq.read_table(release/partition['path'])
        if not table.schema.equals(SCHEMA):
            raise ValueError('Invalid Parquet schema')
        for row in table.to_pylist():
            key=tuple(row[field] for field in ('condition_id','phenotype','source','drug','dose','endpoint'))
            if key in seen or row['endpoint'] not in (-1,0,1) or not row['readout'] or not all(finite(value) for value in row['readout']):
                raise ValueError('Invalid corpus key/endpoint/readout')
            seen.add(key)
            if len(row['annotations'])!=len(row['readout']) or row['condition_id']!=partition['condition_id'] or row['phenotype']!=partition['phenotype'] or json.loads(row['condition'])!=partition['condition']:
                raise ValueError('Broken partition/measurement alignment')
    bundle=json.loads((release/'input.json').read_text())
    with tempfile.TemporaryDirectory(prefix='arc-corpus-verify-') as temp:
        expected=payload(bundle,Path(temp),test_only,json.loads((release/'rules.json').read_text()))
        if expected!=entries:
            raise ValueError('Release differs from frozen evidence/rules/recomputed labels')
    return {'status':'PASS','sample_rows':len(seen),'partitions':len(manifest['partitions']), 'data_status':manifest['data_status']}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    build_parser=sub.add_parser('build')
    build_parser.add_argument('--input',required=True)
    build_parser.add_argument('--output',required=True)
    build_parser.add_argument('--test-only',action='store_true')
    verify_parser=sub.add_parser('verify')
    verify_parser.add_argument('--release',required=True)
    verify_parser.add_argument('--test-only',action='store_true')
    args=parser.parse_args()
    try:
        result=build(args.input,args.output,args.test_only) if args.command=='build' else verify(args.release,args.test_only)
        print(json.dumps(result,sort_keys=True))
    except (ValueError,KeyError,TypeError,OSError,pa.ArrowException) as error:
        print('FAIL: '+str(error),file=sys.stderr)
        return 1
    return 0

if __name__=='__main__':
    sys.exit(main())
