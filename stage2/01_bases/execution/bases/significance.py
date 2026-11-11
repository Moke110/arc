"""Evaluator-only parent-dose truth; training corpora never import this module."""
import numpy as np


def ranked_indices(padj, log2_fold_change):
    """Rank measured significant indices by adjusted P, effect magnitude, index."""
    padj, effect = np.asarray(padj), np.asarray(log2_fold_change)
    if padj.ndim != 1 or effect.shape != padj.shape:
        raise ValueError('invalid DESeq2 statistic shapes')
    finite = np.isfinite(padj)
    if np.any(np.isinf(padj)) or np.any((padj[finite] < 0) | (padj[finite] > 1)):
        raise ValueError('invalid adjusted p-values')
    selected = padj < .05
    if np.any(selected & (~np.isfinite(effect) | (effect == 0))):
        raise ValueError('significant gene has invalid direction')
    lists = []
    for direction in (1, -1):
        indices = np.flatnonzero(selected & (np.sign(effect) == direction))
        order = np.lexsort((indices, -np.abs(effect[indices]), padj[indices]))
        lists.append(indices[order].tolist())
    return tuple(lists)

import hashlib
import json
import pandas as pd


def membership(values):
    return hashlib.sha256(json.dumps(sorted(values), separators=(',', ':')).encode()).hexdigest()


def group_key(parent, dose):
    dose = float(dose)
    if not str(parent) or not np.isfinite(dose) or dose <= 0:
        raise ValueError('invalid parent-dose key')
    return hashlib.sha256(f'{parent}|{dose:.12g}'.encode()).hexdigest()[:24]


def contrast_plan(treated, native):
    """Bind each parent-dose contrast to all retained native wells and exact DMSO controls."""
    if (treated.empty or treated.sample_id.duplicated().any()
            or treated.source_profile_id.duplicated().any() or native.index.has_duplicates):
        raise ValueError('duplicate or missing treated identities')
    rows = treated.source_row.to_numpy(dtype=int)
    if (np.any(rows < 0) or np.any(rows >= len(native)) or len(set(rows)) != len(rows)
            or list(native.index[rows]) != treated.source_profile_id.tolist()):
        raise ValueError('invalid native treated row joins')
    raw = native.iloc[rows]
    if (raw.is_control.any() or raw.is_blank.any()
            or ('drug_key' in treated and list(raw.drug_key) != treated.drug_key.tolist())
            or not np.array_equal(raw.dose_um, treated.dose_um)):
        raise ValueError('native treatment differs from training corpus')
    treated = treated.copy()
    treated['group_id'] = [group_key(p, d) for p, d in zip(treated.drug_parent, treated.dose_um)]
    records, samples = [], []
    for key, part in treated.groupby('group_id', sort=True):
        native_treated = native.iloc[part.source_row.to_numpy(dtype=int)]
        if len(part[['drug_parent', 'dose_um']].drop_duplicates()) != 1:
            raise ValueError('parent-dose key collision')
        for field in ('cell_type', 'cell_line'):
            if field in native_treated and native_treated[field].nunique(dropna=False) != 1:
                raise ValueError(f'incompatible parent-dose {field}: {key}')
        conditions = native_treated[['time_h', 'vehicle']].drop_duplicates()
        if (len(conditions) != 1 or conditions.isna().any().any()
                or str(conditions.iloc[0].vehicle).strip().upper() != 'DMSO'):
            raise ValueError(f'incompatible parent-dose conditions: {key}')
        controls = set()
        for row in part.itertuples(index=False):
            indices = json.loads(row.matched_control_rows)
            names = json.loads(row.matched_control_ids)
            if (not indices or len(set(indices)) != len(indices)
                    or any(not isinstance(i, int) or not 0 <= i < len(native) for i in indices)):
                raise ValueError('invalid DMSO control indices')
            match = native.iloc[indices]
            plate = native.iloc[row.source_row].plate
            if (list(match.index) != names or not match.is_control.all() or match.is_blank.any()
                    or not match.compound_name.str.strip().str.upper().eq('DMSO').all()
                    or not match.vehicle.str.strip().str.upper().eq('DMSO').all()
                    or not match.plate.eq(plate).all()
                    or not match.time_h.eq(native.iloc[row.source_row].time_h).all()):
                raise ValueError('controls are not exact same-plate DMSO')
            expected = native.index[(native.is_control & ~native.is_blank
                & native.compound_name.str.strip().str.upper().eq('DMSO')
                & native.vehicle.str.strip().str.upper().eq('DMSO') & native.plate.eq(plate))]
            if set(names) != set(expected):
                raise ValueError('DMSO control membership is incomplete')
            controls.update(indices)
        controls = sorted(controls)
        plates = sorted(native_treated.plate.astype(str).unique())
        selected = pd.concat([
            part[['sample_id', 'source_profile_id', 'source_row']].assign(condition='treated'),
            pd.DataFrame({'sample_id': ['' for _ in controls],
                'source_profile_id': native.index[controls], 'source_row': controls,
                'condition': 'control'}),
        ], ignore_index=True).sort_values('source_row')
        if selected.source_row.duplicated().any():
            raise ValueError('duplicated contrast sample')
        selected['plate'] = native.iloc[selected.source_row.to_numpy()].plate.to_numpy()
        selected['group_id'] = key
        records.append({'group_id': key, 'drug_parent': str(part.drug_parent.iloc[0]),
            'dose_um': float(part.dose_um.iloc[0]), 'n_treated': len(part),
            'n_controls': len(controls), 'n_plates': len(plates),
            'design': '~ plate + condition' if len(plates) > 1 else '~ condition',
            'time_h': float(conditions.iloc[0].time_h), 'vehicle': str(conditions.iloc[0].vehicle),
            'sample_ids_sha256': membership(part.sample_id.tolist()),
            'source_ids_sha256': membership(part.source_profile_id.tolist()),
            'control_ids_sha256': membership(native.index[controls].tolist())})
        samples.append(selected)
    return pd.DataFrame(records), pd.concat(samples, ignore_index=True)

from pathlib import Path

SCHEMA = 'stage2-parent-dose-truth-v1'
FIELDS = ('pvalue', 'padj', 'log2_fold_change')


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def checked_file(root, name, digest):
    path = Path(root) / name
    if path.resolve() != path.absolute() or not path.is_file() or sha256(path) != digest:
        raise ValueError(f'changed or unsafe truth payload: {path}')
    return path


def bh_adjust(pvalues):
    """Independent complete-family BH replay, including untestable genes in n."""
    values = np.asarray(pvalues)
    finite = np.isfinite(values)
    if np.any(np.isinf(values)) or np.any((values[finite] < 0) | (values[finite] > 1)):
        raise ValueError('invalid raw p-values')
    result = np.full(values.shape, np.nan)
    positions = np.flatnonzero(finite)
    order = positions[np.argsort(values[positions], kind='stable')]
    scaled = values[order] * len(values) / np.arange(1, len(order) + 1)
    result[order] = np.minimum(1, np.minimum.accumulate(scaled[::-1])[::-1])
    return result


def publish_statistics(destination, genes, groups, samples, statistics, provenance, *, fixture=False):
    """Seal DESeq2 outputs into separately keyed, ranked evaluator truth files.

    Statistics must be from raw-count axis-only fits. This publishing seam
    independently replays BH and list ordering; it never estimates significance
    from predicted deltas or expression-normalized training arrays.
    """
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    if (list(genes.gene_index) != list(range(len(genes))) or genes.hgnc_id.duplicated().any()
            or (not fixture and len(genes) != 18164)):
        raise ValueError('invalid indexed significance gene family')
    if (groups.group_id.duplicated().any() or groups[['drug_parent', 'dose_um']].duplicated().any()
            or any(group_key(r.drug_parent, r.dose_um) != r.group_id for r in groups.itertuples())):
        raise ValueError('duplicate or invalid parent-dose keys')
    treated = samples[samples.condition.eq('treated')]
    if (treated.sample_id.duplicated().any() or treated.source_profile_id.duplicated().any()
            or set(treated.group_id) != set(groups.group_id)):
        raise ValueError('invalid released sample memberships')
    if any(np.shape(statistics[field]) != (len(groups), len(genes)) for field in FIELDS):
        raise ValueError('invalid full-axis statistic shape')
    destination.mkdir(parents=True)
    (destination / 'groups').mkdir()
    for name in FIELDS:
        output = np.lib.format.open_memmap(destination / f'{name}.npy', mode='w+',
                                          dtype='<f8', shape=(len(groups), len(genes)))
        for i in range(len(groups)):
            output[i] = statistics[name][i]
        output.flush()
        del output
    for i, row in enumerate(groups.itertuples(index=False)):
        part = treated[treated.group_id.eq(row.group_id)]
        if membership(part.sample_id.tolist()) != row.sample_ids_sha256:
            raise ValueError('released treated membership differs')
        pvalue = np.asarray(statistics['pvalue'][i])
        padj = np.asarray(statistics['padj'][i])
        if not np.allclose(padj, bh_adjust(pvalue), rtol=1e-12, atol=1e-15, equal_nan=True):
            raise ValueError('adjusted p-values differ from full-axis BH')
        eligible = bool(np.isfinite(pvalue).any())
        if row.status not in ('eligible', 'untestable') or eligible != (row.status == 'eligible'):
            raise ValueError('invalid test status')
        if not eligible and not str(row.reason):
            raise ValueError('untestable contrast lacks reason')
        up, down = ranked_indices(padj, statistics['log2_fold_change'][i])
        write_json(destination / 'groups' / f'{row.group_id}.json', {
            'group_id': row.group_id, 'drug_parent': row.drug_parent, 'dose_um': float(row.dose_um),
            'status': row.status, 'reason': str(row.reason), 'up_indices': up, 'down_indices': down,
            'sample_ids_sha256': row.sample_ids_sha256,
        })
    genes.to_csv(destination / 'genes.csv', index=False)
    groups.to_csv(destination / 'groups.csv', index=False)
    samples.to_csv(destination / 'contrast_samples.csv', index=False)
    files = {str(p.relative_to(destination)): sha256(p) for p in sorted(destination.rglob('*')) if p.is_file()}
    write_json(destination / 'manifest.json', {
        'schema': SCHEMA, 'dataset': 'tenv', 'fixture': fixture, 'gene_count': len(genes),
        'group_count': len(groups), 'fdr': .05, 'complete': True,
        'method': 'raw-count DESeq2 two-sided Wald; within-contrast BH over exactly indexed axis',
        'independent_filtering': False, 'prediction_significance': False,
        'provenance': provenance, 'files': files,
    })
    return destination


class PartitionTruth:
    """Evaluator capability bound to exact sample membership of one partition.

    Only the requested role's per-group ranked lists are opened. Ordinary
    training uses bases.corpora, which never imports this evaluator module.
    Assignments can use stable sample_id or source_profile_id/case_id, allowing
    frozen historical cases to remain byte-identical.
    """
    def __init__(self, release, assignments, role, *, allow_fixture=False):
        self.root = Path(release).absolute()
        self.manifest = json.loads((self.root / 'manifest.json').read_text())
        manifest = self.manifest
        if (manifest.get('schema') != SCHEMA or manifest.get('dataset') != 'tenv'
                or manifest.get('fdr') != .05 or not manifest.get('complete')
                or (manifest.get('fixture') and not allow_fixture)
                or (not manifest.get('fixture') and manifest.get('gene_count') != 18164)):
            raise ValueError('invalid evaluator truth release')
        def table(name):
            return pd.read_csv(checked_file(self.root, name, manifest['files'][name]), keep_default_na=False)
        self.groups = table('groups.csv').set_index('group_id', drop=False)
        samples = table('contrast_samples.csv')
        treated = samples[samples.condition.eq('treated')]
        self.genes = table('genes.csv')
        if (self.groups.index.has_duplicates or treated.sample_id.duplicated().any()
                or self.groups[['drug_parent', 'dose_um']].duplicated().any()
                or list(self.genes.gene_index) != list(range(manifest['gene_count']))):
            raise ValueError('invalid truth release keys')
        assignments = assignments.copy()
        if 'sample_id' not in assignments:
            column = 'source_profile_id' if 'source_profile_id' in assignments else 'case_id'
            mapping = treated.set_index('source_profile_id').sample_id
            assignments['sample_id'] = assignments[column].map(mapping)
        if (assignments.empty or assignments.sample_id.isna().any()
                or assignments.sample_id.duplicated().any()
                or not set(assignments.role).issubset({'train', 'validation', 'test'})
                or role not in {'train', 'validation', 'test'}):
            raise ValueError('invalid evaluator partition assignments')
        mapping = treated.set_index('sample_id').group_id
        assignments['group_id'] = assignments.sample_id.map(mapping)
        if assignments.group_id.isna().any():
            raise ValueError('partition contains unknown samples')
        for key, part in assignments.groupby('group_id'):
            if part.role.nunique() != 1:
                raise ValueError(f'contrast crosses partitions: {key}')
            if membership(part.sample_id.tolist()) != self.groups.loc[key].sample_ids_sha256:
                raise ValueError(f'partition has incomplete contrast: {key}')
        self.role = role
        self.memberships = {key: sorted(part.sample_id.tolist()) for key, part in
                            assignments[assignments.role.eq(role)].groupby('group_id')}
        self._source_ids = treated.set_index('source_profile_id').sample_id.to_dict()
        self._labels = {}
        for key in self.memberships:
            name = f'groups/{key}.json'
            value = json.loads(checked_file(self.root, name, manifest['files'][name]).read_text())
            row = self.groups.loc[key]
            up, down = value['up_indices'], value['down_indices']
            indices = [*up, *down]
            if (value['group_id'] != key or value['drug_parent'] != row.drug_parent
                    or float(value['dose_um']) != float(row.dose_um)
                    or group_key(row.drug_parent, row.dose_um) != key
                    or value['sample_ids_sha256'] != row.sample_ids_sha256
                    or value['status'] != row.status
                    or value['status'] not in ('eligible', 'untestable')
                    or len(set(indices)) != len(indices)
                    or any(type(i) is not int or not 0 <= i < manifest['gene_count'] for i in indices)
                    or (value['status'] == 'untestable' and (indices or not value['reason']))):
                raise ValueError('invalid partition truth list')
            self._labels[key] = value

    def ranked(self, parent, dose):
        key = group_key(parent, dose)
        if key not in self._labels:
            raise ValueError('truth group is outside the evaluator partition')
        value = self._labels[key]
        return list(value['up_indices']), list(value['down_indices'])

    def truth(self, parent, dose, sample_ids):
        key = group_key(parent, dose)
        if key not in self._labels:
            raise ValueError('truth group is outside the evaluator partition')
        resolved = [self._source_ids.get(s, s) for s in sample_ids]
        if sorted(resolved) != self.memberships[key]:
            raise ValueError('truth contrast sample membership differs from evaluator partition')
        value = self._labels[key]
        if value['status'] == 'untestable':
            return None
        truth = np.zeros(self.manifest['gene_count'], dtype=np.int8)
        truth[value['up_indices']] = 1
        truth[value['down_indices']] = -1
        return truth
