"""DEG-independent small-molecule TOX admission (ADR 0072).

This ledger is an input to analysis, never inferred from a DEG table or ATC/MOA.
Names are resolved against the pinned structural annotation; evidence transfers
between sources only for the same full parent InChIKey. Missing or conflicting
labels stay pending. The old mechanism-class API remains readable until release.
"""
from pathlib import Path
import re
import pandas as pd
from .release import ReleaseError, read_json, sha256, write_json

EVIDENCE = 'stage2/01_bases/execution/tox_drug_evidence.json'
SOURCES = 'stage2/01_bases/execution/tox_source_admission.json'
ANNOTATION = 'pert_data/Annotation/local/cm_compound_annotation.csv'
NON_SMALL_MOLECULES = {'bevacizmab', 'bevacizumab', 'cetuximab', 'rituximab',
                       'trastuzumab', 'endothelin', 'igf1', 'tnfalpha', 'carfilzomib'}


def name_key(value):
    return re.sub('[^a-z0-9]', '', str(value).casefold())


def assess_drugs(records, annotation, evidence):
    """Return one auditable admission decision per source code, before DEG work."""
    rows = []
    for record in records:
        row = dict(record)
        name = name_key(row['drug_name'])
        matches = annotation[annotation.lookup_name.map(name_key).eq(name)]
        if not matches.empty:
            own = matches[matches.dataset_inchikey.fillna('').ne('')]
            if not own.empty:
                matches = own
        keys = set(matches.parent_inchikey.fillna('')) - {''}
        row.update(parent_inchikey='', small_molecule=False, evidence_label='pending',
                   evidence=[], disposition='excluded', reason='unresolved chemical identity',
                   atc='', moa='', uncertainty='no independent positive evidence established')
        if len(keys) == 1:
            key = next(iter(keys))
            row['parent_inchikey'] = key
            row['atc'] = '; '.join(sorted(set(matches.atc.fillna('')) - {''}))
            row['moa'] = '; '.join(sorted(set(matches.moa.fillna('')) - {''}))
            row['small_molecule'] = name not in NON_SMALL_MOLECULES
            facts = [dict(e) for e in evidence if e.get('parent_inchikey') == key]
            for fact in facts:
                if not all(fact.get(k) for k in ('url', 'locator', 'retrieved_at', 'rationale', 'uncertainty')):
                    raise ReleaseError('incomplete independent TOX evidence')
                if fact['label'] not in ('positive', 'negative', 'pending'):
                    raise ReleaseError('invalid TOX evidence label')
            labels = {e['label'] for e in facts} - {'pending'}
            row['evidence'] = facts
            row['evidence_label'] = next(iter(labels)) if len(labels) == 1 else 'pending'
            row['uncertainty'] = '; '.join(sorted({e['uncertainty'] for e in facts})) or row['uncertainty']
            if not row['small_molecule']:
                row['reason'] = 'antibody or peptide; outside small-molecule TOX definition'
            elif row['evidence_label'] == 'positive':
                row.update(disposition='admitted', reason='identity-resolved drug with independent positive cardiotoxicity evidence')
            elif row['evidence_label'] == 'negative':
                row['reason'] = 'explicit negative reference drug; never pooled as TOX'
            else:
                row.update(disposition='pending', reason='conflicting evidence' if len(labels) > 1 else 'insufficient independent cardiotoxicity evidence')
        rows.append(row)
    return pd.DataFrame(rows).sort_values(['source', 'source_id'], kind='stable').reset_index(drop=True)


def build_ledger(repo, output):
    """Build a pinned admission artifact without loading expression or DEG data."""
    from . import cardiotox
    repo, output = Path(repo), Path(output)
    evidence = read_json(repo / EVIDENCE)
    sources = read_json(repo / SOURCES)
    annotation = pd.read_csv(repo / ANNOTATION, dtype=str, keep_default_na=False)
    records = []
    names = cardiotox.dtoxs_series_names(repo / cardiotox.DTOXS_SERIES)
    groups = cardiotox.dtoxs_drug_groups(repo / cardiotox.DTOXS_PARAMS)
    records.extend(dict(source='dtoxs_ipsc', source_id=code, drug_name=names[code])
                   for code in groups.Drug.astype(str))
    panel = cardiotox.gse244740_panel(repo / cardiotox.GSE244740_SERIES)
    records.extend(dict(source='gse244740', source_id=name, drug_name=name)
                   for name in sorted(set(panel.compound) - {'DMSO', 'NoTreatment'}))
    for source in sources['sources']:
        records.extend(dict(source=source['source'], **drug) for drug in source.get('drugs', []))
    frame = assess_drugs(records, annotation, evidence['drugs'])
    inputs = [EVIDENCE, SOURCES, ANNOTATION, cardiotox.DTOXS_SERIES,
              cardiotox.DTOXS_PARAMS, cardiotox.GSE244740_SERIES]
    inputs += evidence['pinned_sources']
    inputs += [p for source in sources['sources'] for p in source.get('inputs', [])]
    value = dict(schema='stage2-tox-admission-v1', disease='TOX',
                 definition='transcriptional response to independently supported cardiotoxic small-molecule drugs; no toxicity assertion at each exposure',
                 exposure_policy='all experiment-quality doses/times with matched vehicle; no mechanism or minimum-drug restriction',
                 inputs={p: {'sha256': sha256(repo / p)} for p in sorted(set(inputs))},
                 counts=frame.disposition.value_counts().sort_index().to_dict(),
                 sources=sources['sources'], discovery=sources['discovery'],
                 drugs=frame.to_dict(orient='records'))
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / 'admission.json', value)
    return value


class ToxAdmission:
    """Digest-checkable consumer seam for the frozen pre-analysis decision ledger."""
    def __init__(self, path, repo=None):
        self.value = read_json(Path(path) / 'admission.json')
        if self.value['schema'] != 'stage2-tox-admission-v1':
            raise ReleaseError('unknown TOX admission schema')
        if repo is not None:
            for path, spec in self.value['inputs'].items():
                if sha256(Path(repo) / path) != spec['sha256']:
                    raise ReleaseError(f'TOX admission input changed: {path}')

    def drugs(self, source=None, admitted_only=False):
        frame = pd.DataFrame(self.value['drugs'])
        if source is not None:
            frame = frame[frame.source.eq(source)]
        if admitted_only:
            frame = frame[frame.disposition.eq('admitted')]
        return frame.reset_index(drop=True)
