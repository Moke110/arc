#!/usr/bin/env python3
"""Extract the full frozen Sirenko 24 h beat-count workbook offline."""
import argparse
from collections import Counter
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import xlrd


def root_path():
    for parent in (Path(__file__).resolve().parent, *Path(__file__).resolve().parents):
        if (parent / '.git').exists():
            return parent
    raise ValueError('Repository root not found')


def number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def json_bytes(value):
    return (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n').encode()


def extract(root):
    item = root / 'pheno_data/Sirenko-2013'
    asset = 'pheno_data/Sirenko-2013/source/publisher-mmc3.xls'
    sheet = xlrd.open_workbook(root / asset).sheet_by_name('24hr_beat_count')
    if sheet.ncols != 22 or sheet.cell_value(0, 3) != 'Chemical_Name':
        raise ValueError('Unexpected source workbook schema')
    records, identities, anomalies = [], [], []
    coverage = Counter()
    for index in range(1, sheet.nrows):
        row = sheet.row_values(index)
        excel_row = index + 1
        identity = {
            'source_row': excel_row, 'plate_location': row[0], 'catalog_number': row[1],
            'original_cas': row[2], 'administered_form_name': row[3].strip(),
            'library_Control': row[4], 'canonical_parent': None,
            'identity_status': 'unresolved_parent_and_stereochemistry',
            'pharmacological_qualification': 'unknown',
            'qualification_reason': 'Library cardiotoxicity class does not establish therapeutic/tool use.',
            'identity_evidence': f'{asset}#24hr_beat_count!A{excel_row}:E{excel_row}',
        }
        identities.append(identity)
        for dose_index, dose in enumerate((.3, 3., 10., 30.)):
            raw_columns = (5 + 2*dose_index, 6 + 2*dose_index)
            normalized_columns = (14 + 2*dose_index, 15 + 2*dose_index)
            raw = [row[column] if number(row[column]) else None for column in raw_columns]
            normalized = [row[column] if number(row[column]) else None for column in normalized_columns]
            coverage['comparisons'] += 1
            coverage['raw_at_least_one'] += any(value is not None for value in raw)
            coverage['raw_complete_pairs'] += all(value is not None for value in raw)
            coverage['raw_finite_values'] += sum(value is not None for value in raw)
            coverage['normalized_finite_values'] += sum(value is not None for value in normalized)
            coverage['raw_blank_cells'] += sum(value is None for value in raw)
            coverage['raw_double_zero_comparisons'] += raw == [0, 0]
            measurements = []
            comparison_id = f'sirenko-24hr-beat-r{excel_row}-dose{dose:g}'
            for replicate, (raw_column, normalized_column) in enumerate(zip(raw_columns, normalized_columns)):
                selector = f'24hr_beat_count!{xlrd.colname(raw_column)}{excel_row}'
                normalized_selector = f'24hr_beat_count!{xlrd.colname(normalized_column)}{excel_row}'
                if raw[replicate] is None and normalized[replicate] == 0:
                    anomalies.append({'raw_selector': selector, 'normalized_selector': normalized_selector,
                                      'source_name': identity['administered_form_name'], 'dose_uM': dose})
                measurements.append({
                    'measurement_id': comparison_id + f'-r{replicate+1}',
                    'value': raw[replicate], 'normalized_value': normalized[replicate],
                    'assay_method': 'FLIPR_Calcium5', 'unit': 'beat_count_per_source_acquisition',
                    'representation': 'raw', 'asset': asset, 'selector': selector,
                    'normalized_selector': normalized_selector, 'normalized_unit': 'source_normalized_scale',
                    'experimental_unit': None, 'replicate_column': f'r{replicate+1}',
                    'replication_design': 'source_duplicate_measurements_independence_unknown',
                    'zero_status': 'arrest_or_technical_validity_unresolved' if raw[replicate] == 0 else None,
                    'original_identity': identity,
                })
            records.append({
                'comparison_id': comparison_id, 'source': 'sirenko_2013',
                'condition': {
                    'condition_cell_type': 'hiPSC-CM_iCell', 'condition_species': 'Homo sapiens',
                    'condition_drug_time': 24., 'condition_hypoxia_time': 0.,
                    'condition_reoxygen_time': 0., 'condition_drug_schedule': 'direct_exposure',
                    'condition_vehicle': 'DMSO', 'condition_vehicle_percent': .3 if dose == 30 else .1,
                    'condition_maturation': '8-10_days_post_plating',
                    'condition_acquisition_duration': 'up_to800frames_8fps_exact_window_unknown',
                    'condition_unresolved_context': 'sirenko_2013',
                },
                'phenotype': 'beat_count', 'assay_method': 'FLIPR_Calcium5',
                'batch': None,
                'drug': None, 'dose': float(Decimal(str(dose))*Decimal('1e-6')), 'original_dose': dose,
                'original_dose_unit': 'uM', 'dose_basis': 'source_administered_form',
                'identity_evidence': identity['identity_evidence'],
                'dose_evidence': 'article.xml P18 four24h concentrations; source dose headers; uM*1e-6 mol/L; parent basis unresolved',
                'control_evidence': 'article.xml P7/P9/P18 plate-specific DMSO wells; source workbook publishes no vehicle means/variance; Control is a library flag',
                'measurements': measurements,
                'quality_reasons': ['pharmacological_qualification_unknown', 'parent_dose_basis_unresolved'] +
                    (['zero_measurement_validity_unresolved'] if 0 in raw else []),
                'statistics': {'treated_mean': None, 'control_mean': None,
                    'scale_compatible': False, 'control_matched': False,
                    'p_raw': None, 'p_adjusted': None, 'family_complete': False,
                    'evidence': 'No numeric dosewise test, control moments, exact normalization formula or independent-unit design; no guessed100 denominator.'},
            })
    audit = {
        'dataset_id': 'sirenko_2013', 'compound_labels': len(identities),
        'unique_original_names': len({row['administered_form_name'] for row in identities}),
        'unique_original_cas_values': len({str(row['original_cas']) for row in identities}),
        **dict(coverage), 'raw_missing_normalized_zero': anomalies,
        'canonical_drug_count': None, 'pharmacological_qualification_counts': {'unknown': len(identities)},
        'identity_status': 'Source labels retained; full parent/use curation remains unresolved.',
        'scientific_status': 'All comparisons lack validated FC/control; exclusions preserve all source values.',
        'model_view_500_assessment': 'not_run', 'training_status': 'not_trained',
    }
    sources = []
    for line in (item/'manifest.sha256').read_text().splitlines():
        sha, relative = line.split('  ', 1)
        sources.append({'path': str(Path('pheno_data/Sirenko-2013/source')/relative.removeprefix('source/')), 'sha256': sha})
    provenance = json.loads((item/'author_log.json').read_text())
    for addition in provenance.get('local_additions', []):
        path = item/addition['path']
        if hashlib.sha256(path.read_bytes()).hexdigest() != addition['sha256']:
            raise ValueError('Frozen local curation hash mismatch')
        sources.append({'path': path.relative_to(root).as_posix(), 'sha256': addition['sha256']})
    return {'schema_version': 1, 'sources': [{'source': 'sirenko_2013', 'item': 'pheno_data/Sirenko-2013',
            'test_only': False, 'status': 'archived', 'assets': sources}], 'records': records,
            'source_audits': {'sirenko_2013': audit}}, identities


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--freeze-curation', action='store_true',
                        help='Declare immutable v1 identity and coverage extractions in the source item.')
    args = parser.parse_args()
    root = root_path()
    bundle, identities = extract(root)
    if args.freeze_curation:
        item = root/'pheno_data/Sirenko-2013'
        log_path = item/'author_log.json'
        log = json.loads(log_path.read_text())
        for name, data in [('identities-v1.json', identities),
                           ('source-audit-v1.json', bundle['source_audits']['sirenko_2013'])]:
            path = item/'local'/name
            payload = json_bytes(data)
            if path.exists() and path.read_bytes() != payload:
                raise ValueError('Immutable curation conflict: choose a new version')
            path.parent.mkdir(exist_ok=True)
            path.write_bytes(payload)
            relative = path.relative_to(item).as_posix()
            if not any(entry['path'] == relative for entry in log['local_additions']):
                log['local_additions'].append({'path': relative, 'bytes': len(payload),
                    'sha256': hashlib.sha256(payload).hexdigest(),
                    'origin': 'Offline extraction of frozen Supplemental Table3 by stage3/01_data_corpus/execution/sirenko_adapter.py; unresolved parents and usage explicitly retained',
                    'version': 'v1'})
        log_path.write_bytes(json_bytes(log))
        bundle, identities = extract(root)
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(json_bytes(bundle))
    print(json.dumps(bundle['source_audits'], sort_keys=True))


if __name__ == '__main__':
    main()
