#!/usr/bin/env python3
"""Audit frozen fold 0 Tahoe exposure and shared-parent responses."""
import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

sys.dont_write_bytecode = True
REPO = next(path for path in Path(__file__).resolve().parents if (path / '.git').exists())
STEP05 = Path('stage2/05_pdiff_transfer')
sys.path.insert(0, str(REPO / STEP05 / 'execution/script'))

import numpy as np

from cold_tahoe import FITTED, PATHS, ROLES, STEP04, frozen_inputs, pin, verify_files
from artifacts import read_json, rows, sha256, write_csv, write_json

BRACKETS = {0.2: (0.05, 0.5), 1.0: (0.5, 5.0)}
TEN_ONLY = 10.0


@dataclass(frozen=True)
class Release:
    root: Path
    hashes: dict

    def path(self, key):
        verify_files(self.root, {key: PATHS[key]}, self.hashes)
        return self.root / PATHS[key]


def case_tables(release, roles):
    drugs = {row['drug_key']: row for row in rows(release.path('drugs'))}
    return {
        role: [{**case, **drugs[case['drug_key']]}
               for case in rows(release.path(f'{role}_cases'))]
        for role in roles
    }


def by_parent(cases):
    result = defaultdict(list)
    for case in cases:
        result[case['parent_connectivity']].append(case)
    return result


def exposure_rows(tenv, historical, cold):
    target = {role: by_parent(tenv[role]) for role in ROLES}
    forbidden = target['validation'].keys() | target['test'].keys()
    for role in FITTED:
        if cold[role] != [case for case in historical[role]
                          if case['parent_connectivity'] not in forbidden]:
            raise ValueError(f'globally cold Tahoe {role} membership or order differs from exclusion rule')
    result = []
    for source, tables in (('historical', historical), ('globally_cold', cold)):
        for source_role, cases in tables.items():
            source_parents = by_parent(cases)
            for target_role, target_parents in target.items():
                shared = sorted(source_parents.keys() & target_parents.keys())
                source_cases = [case for parent in shared for case in source_parents[parent]]
                target_cases = [case for parent in shared for case in target_parents[parent]]
                result.append(dict(
                    source_release=source, tahoe_role=source_role, tenv_role=target_role,
                    tahoe_parent_count=len(source_parents), tenv_parent_count=len(target_parents),
                    shared_parent_count=len(shared), shared_parents='|'.join(shared),
                    tahoe_shared_profiles=len(source_cases),
                    tahoe_shared_groups=len({case['group_id'] for case in source_cases}),
                    tenv_shared_wells=len(target_cases),
                    tenv_shared_groups=len({case['group_id'] for case in target_cases}),
                ))
    return result


def group_profiles(release, role, cases, shared, n_genes):
    """Average only the shared-parent Δ rows, preserving frozen row offsets."""
    members = defaultdict(list)
    for index, case in enumerate(cases):
        parent = case['parent_connectivity']
        if parent in shared:
            key = (parent, case['drug_key'], case['parent_smiles'], case.get('cell_line_id', ''),
                   float(case['dose_um']))
            members[key].append(index)
    matrix = np.load(release.path(f'{role}_expression'), mmap_mode='r', allow_pickle=False)
    if matrix.shape != (len(cases), n_genes) or matrix.dtype != np.float32:
        raise ValueError(f'invalid {role} Δ array for {release.root}')
    profiles = defaultdict(dict)
    for (*context, dose), indices in sorted(members.items()):
        total = np.zeros(n_genes, dtype=np.float64)
        for index in indices:
            total += matrix[index]
        mean = (total / len(indices)).astype(np.float32)
        if not np.isfinite(mean).all():
            raise ValueError(f'nonfinite mean Δ in {role} for {context}, {dose}')
        profiles[context[0]].setdefault((*context[1:], role), {})[dose] = (mean, len(indices))
    return profiles


def amplitude(values):
    return float(np.sqrt(np.mean(np.square(values, dtype=np.float64))))


def comparison(tenv_profile, tahoe_profile, *, tenv_only=False):
    tenv_delta, tenv_support = tenv_profile or (None, 0)
    tahoe_delta, tahoe_support = tahoe_profile or (None, 0)
    joint_count, concordance = None, None
    if tenv_only:
        status = 'tenv_only' if tenv_support else 'missing_tenv_dose'
    elif not tenv_support or not tahoe_support:
        status = ('missing_both_doses' if not tenv_support and not tahoe_support else
                  'missing_tenv_dose' if not tenv_support else 'missing_tahoe_dose')
    else:
        joint = (tenv_delta != 0) & (tahoe_delta != 0)
        joint_count = int(np.count_nonzero(joint))
        status = 'no_joint_nonzero_genes'
        if joint_count:
            status = 'observed'
            concordance = float(np.mean(np.signbit(tenv_delta[joint]) == np.signbit(tahoe_delta[joint])))
    return dict(
        tenv_wells=tenv_support, tahoe_profiles=None if tenv_only else tahoe_support,
        tenv_delta_rms=amplitude(tenv_delta) if tenv_support else None,
        tahoe_delta_rms=amplitude(tahoe_delta) if tahoe_support else None,
        joint_nonzero_genes=joint_count, sign_concordance=concordance, status=status,
    )


def delta_panel(tenv, cold, tenv_train, cold_cases, shared, n_genes):
    target_profiles = group_profiles(tenv, 'train', tenv_train, shared, n_genes)
    source_profiles = {role: group_profiles(cold, role, cold_cases[role], shared, n_genes)
                       for role in FITTED}
    for parent in sorted(shared):
        tahoe_contexts = {context: doses for role in FITTED
                          for context, doses in source_profiles[role][parent].items()}
        for (drug, smiles, _, _), tenv_doses in target_profiles[parent].items():
            identity = dict(parent_connectivity=parent, tenv_drug_key=drug, tenv_parent_smiles=smiles)
            for (tahoe_drug, tahoe_smiles, line, role), tahoe_doses in tahoe_contexts.items():
                for tenv_dose, endpoints in BRACKETS.items():
                    for label, tahoe_dose in zip(('lower', 'upper'), endpoints, strict=True):
                        yield dict(
                            **identity, tenv_dose_um=tenv_dose, bracket=label, tahoe_dose_um=tahoe_dose,
                            tahoe_role=role, tahoe_cell_line_id=line, tahoe_drug_key=tahoe_drug,
                            tahoe_parent_smiles=tahoe_smiles, exact_parent_smiles=smiles == tahoe_smiles,
                            **comparison(tenv_doses.get(tenv_dose), tahoe_doses.get(tahoe_dose)),
                        )
            yield dict(
                **identity, tenv_dose_um=TEN_ONLY, bracket='tenv_only', tahoe_dose_um=None,
                tahoe_role=None, tahoe_cell_line_id=None, tahoe_drug_key=None, tahoe_parent_smiles=None,
                exact_parent_smiles=any(smiles == context[1] for context in tahoe_contexts),
                **comparison(tenv_doses.get(TEN_ONLY), None, tenv_only=True),
            )


def write_audit(output, tenv, historical, cold, n_genes, *, frozen=None, runtime=REPO):
    """Publish an audit in a fresh directory from read-only input releases."""
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    releases = {'tenv': tenv, 'historical_tahoe': historical, 'globally_cold_tahoe': cold}
    genes = {sha256(release.path('genes')) for release in releases.values()}
    if len(genes) != 1 or len(rows(tenv.root / PATHS['genes'])) != n_genes:
        raise ValueError('source and target must use the same ordered gene axis')

    tenv_cases = case_tables(tenv, ROLES)
    historical_cases = case_tables(historical, FITTED)
    cold_cases = case_tables(cold, FITTED)
    exposure = exposure_rows(tenv_cases, historical_cases, cold_cases)
    shared = by_parent(tenv_cases['train']).keys() & {
        case['parent_connectivity'] for cases in cold_cases.values() for case in cases
    }
    panel = list(delta_panel(tenv, cold, tenv_cases['train'], cold_cases, shared, n_genes))
    exact = [row for row in panel if row['exact_parent_smiles']]
    exact_parents = {row['parent_connectivity'] for row in exact}

    if frozen:
        observed = {
            (row['source_release'], row['tahoe_role'], row['tenv_role']): row['shared_parent_count']
            for row in exposure
        }
        expected = {
            (source, role, target): count
            for (source, role), counts in {
                ('historical', 'train'): (6, 1, 2),
                ('historical', 'validation'): (3, 0, 0),
                ('globally_cold', 'train'): (6, 0, 0),
                ('globally_cold', 'validation'): (3, 0, 0),
            }.items()
            for target, count in zip(ROLES, counts, strict=True)
        }
        if (observed, len(shared), len(exact_parents)) != (expected, 9, 7):
            raise ValueError('audit coverage differs from the accepted fold 0 design')

    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent, prefix='.source-audit-') as temporary:
        pending = Path(temporary)
        write_csv(pending / 'exposure.csv', exposure)
        write_csv(pending / 'delta_panel.csv', panel)
        write_csv(pending / 'exact_smiles_panel.csv', exact, fields=list(panel[0]))
        write_json(pending / 'lineage.json', {
            'schema': 'step06-source-audit-v1', 'fold': 0, 'descriptive_only': True,
            'n_genes': n_genes, 'gene_axis_sha256': genes.pop(),
            'shared_training_parents': sorted(shared),
            'exact_parent_smiles_training_parents': sorted(exact_parents),
            'source': frozen,
            'inputs': {
                name: {'root': str(release.root.relative_to(runtime) if release.root.is_relative_to(runtime)
                                   else release.root),
                       'files': {PATHS[key]: release.hashes[key] for key in input_keys}}
                for (name, release), input_keys in zip(releases.items(), (
                    ['genes', 'drugs', 'train_cases', 'validation_cases', 'test_cases', 'train_expression'],
                    ['genes', 'drugs', 'train_cases', 'validation_cases'],
                    ['genes', 'drugs', 'train_cases', 'validation_cases', 'train_expression',
                     'validation_expression'],
                ), strict=True)
            },
            'script_sha256': sha256(Path(__file__)),
            'outputs': {path.name: sha256(path) for path in sorted(pending.iterdir())},
        })
        pending.rename(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--runtime-root', type=Path, default=REPO)
    args = parser.parse_args()
    runtime = args.runtime_root.resolve()
    step05 = runtime / STEP05 / 'execution'
    step04 = runtime / STEP04 / 'execution'
    split = step05 / 'output/cold-tenv-ticket76/split.json'
    historical_root, historical, target, _, source = frozen_inputs(
        runtime, step04 / 'output/tahoe-formal-ticket73', split)
    tenv = read_json(split)
    if tenv['n_genes'] != historical['n_genes'] or any(
            tenv['predictor_files'][key] != str(Path(tenv['input_root']) / path)
            for key, path in PATHS.items()):
        raise ValueError('unexpected frozen tenv file contract')

    cold_root = step05 / 'output/cold-tahoe-ticket79'
    cold_hash = read_json(step05 / 'cold_tahoe_record/qualification.json')['local_artifact_sha256'][
        'cold-tahoe-ticket79/release.json']
    cold = read_json(pin(cold_root, 'release.json', cold_hash))
    if (cold['tahoe_input_identity'] != source['tahoe_input_identity']
            or cold['tenv_split_sha256'] != source['tenv_split_sha256']
            or set(cold['forbidden_parents']) != target['validation'] | target['test']):
        raise ValueError('globally cold release has different predecessor inputs or exclusions')

    write_audit(
        args.output, Release(Path(tenv['input_root']), tenv['file_sha256']),
        Release(historical_root, historical['file_sha256']), Release(cold_root, cold['file_sha256']),
        tenv['n_genes'], frozen={
            'tenv_split_sha256': source['tenv_split_sha256'],
            'tenv_split_identity': source['tenv_split_identity'],
            'historical_tahoe_plan_sha256': source['tahoe_plan_sha256'],
            'historical_tahoe_request_sha256': source['tahoe_request_sha256'],
            'historical_tahoe_input_identity': source['tahoe_input_identity'],
            'globally_cold_release_sha256': cold_hash,
        }, runtime=runtime,
    )


if __name__ == '__main__':
    main()
