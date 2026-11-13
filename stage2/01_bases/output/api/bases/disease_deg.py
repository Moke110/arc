"""Disease-DEG publication and reader for the Stage 2 gene axis (ADR 0066).

Tickets #187/#188/#189 publish, under ``output/disease-deg/``, one compact
Control-vs-Disease table and one within-dataset rank per source dataset; the
fusion ticket (#190) consumes those level-1 inputs. This module is the single
public boundary for that product family:

* ``align_source_table`` — attach the frozen ``gene_index`` to an item-local DEG
  table (symbol first, then source identifier/Ensembl) through
  :class:`bases.deg_rank.SymbolResolver`; unresolved rows are reported, never
  guessed.
* ``build_dataset_products`` — collapse an item-local per-contrast DEG table with
  the ticket-01 contract (:func:`bases.deg_rank.dataset_deg_rank`), carry the
  winning contrast's ``stat``/``p_value``/``base_mean`` and the frozen
  ``gene_symbol``, and return the compact table, the full-axis rank frame and the
  statistics.
* ``genotype_split`` — the genotype-first donor split the Reichart source
  declares (positive donors leave ``nHF``; the single NCCM donor is not used).
* ``contrast_selection`` — the one-joint-design arm collapse: one fit produces
  several contrast rows, each published dataset selects its declared labels.
* ``DiseaseDeg`` — the reader every later ticket imports. It returns per
  category, and per dataset inside it, the compact table, the rank frame and the
  declared eligibility / shared-control record, all digest-checked against the
  category release manifest.

Layout (one product directory per disease category)::

    output/disease-deg/<category>/
        catalog.json                 # category -> dataset -> paths + eligibility
        manifest.json                # candidate-then-replace release manifest
        <dataset>.deg.csv            # compact table, one row per measured gene
        <dataset>.rank.csv           # frozen-axis rank frame (18,164 rows)

``.deg.csv`` columns are ``DEG_COLUMNS``; ``.rank.csv`` columns are
``RANK_COLUMNS``; both are asserted by :func:`validate_table` /
:func:`validate_rank`. Ranks are within-dataset only and must never be read as
effect sizes across datasets.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from .deg_rank import (
    AXIS_SIZE,
    SCHEMA as RANK_SCHEMA,
    SymbolResolver,
    align_table,
    dataset_deg_rank,
    deg_effect,
)

CATALOG_SCHEMA = 'stage2-disease-deg-catalog-v1'
RELEASE_SCHEMA = 'stage2-disease-deg-release-v1'
LOCAL_SCHEMA = 'stage2-item-deg-v1'

LOCAL_DEG_COLUMNS = ('source_gene_id', 'gene_symbol', 'contrast', 'log2fc', 'stat',
                     'p_value', 'padj', 'base_mean')
DEG_COLUMNS = ('gene_index', 'gene_symbol', 'contrast', 'log2fc', 'stat', 'p_value',
               'padj', 'base_mean', 'direction', 'significant', 'n_contrasts',
               'n_contrasts_supporting', 'effect', 'rank', 'rank_norm')
RANK_COLUMNS = ('gene_index', 'rank', 'rank_norm')

GENOTYPE_GROUPS = ('dcm_genotype_negative', 'dcm_genotype_positive', 'arvc',
                   'normal', 'nccm')


def read_axis(path):
    """Read the frozen 18,164-gene identity table (``gene-identity/genes.csv``)."""
    axis = pd.read_csv(path, dtype={'symbol': str, 'hgnc_id': str,
                                    'ensembl_id': str, 'match_type': str})
    required = {'gene_index', 'symbol', 'hgnc_id', 'ensembl_id'}
    if not required <= set(axis.columns):
        raise ValueError(f'axis table requires {sorted(required)}')
    if len(axis) != AXIS_SIZE:
        raise ValueError(f'axis must hold {AXIS_SIZE} genes, found {len(axis)}')
    if axis.gene_index.tolist() != list(range(AXIS_SIZE)):
        raise ValueError('axis gene_index is not the frozen 0..18163 order')
    return axis


def build_resolver(axis, crosswalk=()):
    """The ticket-01 resolver over the frozen axis and pinned crosswalk evidence.

    The pinned HGNC complete set is deliberately not an input here: it is a
    gitignored external source, so using it would make the published tables
    depend on whether that tree happens to be materialised. Canonical symbol,
    Ensembl id and the tracked crosswalk's alias/previous rows are enough for
    this ticket's sources, and every loss is reported.
    """
    return SymbolResolver(axis, crosswalk=crosswalk)


def read_crosswalk(path):
    if path is None or not Path(path).is_file():
        return pd.DataFrame(columns=['dataset_id', 'source_gene_id', 'hgnc_id',
                                     'match_type'])
    return pd.read_csv(path, dtype=str)


def read_local_deg(path):
    """Read one item-local full DEG table (``local/deg/*_deg.csv``)."""
    frame = pd.read_csv(path)
    missing = set(LOCAL_DEG_COLUMNS) - set(frame.columns)
    if missing:
        raise ValueError(f'{path}: local DEG table lacks {sorted(missing)}')
    if frame.empty:
        raise ValueError(f'{path}: local DEG table is empty')
    return frame


def align_source_table(rows, resolver):
    """Attach ``gene_index`` to item-local DEG rows, symbol first then identifier.

    A row's symbol resolves through the canonical/alias/previous/Ensembl
    precedence; when the symbol is empty or unresolved, the source identifier is
    tried (an Ensembl id resolves through the axis or crosswalk). The returned
    table keeps only resolved rows and carries ``alignment``; the report counts
    every category, so an unresolved row is visible instead of silently dropped.
    """
    rows = pd.DataFrame(rows)
    if rows.empty:
        raise ValueError('an item DEG table needs at least one row')
    symbols = rows['gene_symbol'].astype(str).to_numpy(dtype=object)
    identifiers = rows['source_gene_id'].astype(str).to_numpy(dtype=object)
    indices, categories = resolver.resolve_many(symbols.tolist())
    fallback_index, fallback_category = resolver.resolve_many(identifiers.tolist())
    resolved = np.asarray(indices, dtype=np.int64)
    category = np.asarray(categories, dtype=object)
    fallback_resolved = np.asarray(fallback_index, dtype=np.int64) >= 0
    take = (resolved < 0) & fallback_resolved & (category != 'ambiguous')
    resolved[take] = np.asarray(fallback_index, dtype=np.int64)[take]
    category[take] = np.asarray(fallback_category, dtype=object)[take]
    aligned = rows.assign(gene_index=resolved, alignment=category)
    first_pass = resolver.loss_summary(indices, categories)
    final = resolver.loss_summary(resolved, category)
    report = {
        'n_source_rows': int(len(rows)),
        'n_resolved': int((resolved >= 0).sum()),
        'n_axis_genes_resolved': int(len(set(resolved[resolved >= 0].tolist()))),
        'alignment_categories': {name: int(sum(1 for value in category if value == name))
                                 for name in sorted(set(category.tolist()))},
        # A symbol that is ambiguous under the crosswalk stays ambiguous: the
        # identifier fallback never silently re-labels it.
        'n_recovered_by_identifier': int(take.sum()),
        'n_ambiguous_kept': int(sum(1 for value in category if value == 'ambiguous')),
        'symbol_first_pass': first_pass,
    }
    table = aligned[aligned.gene_index >= 0].copy()
    table['gene_index'] = table['gene_index'].astype(np.int64)
    return table.reset_index(drop=True), report


def winning_contrast_rows(aligned):
    """One row per (gene_index, contrast): the maximum-effect duplicate wins.

    This is the same rule ``deg_rank.per_gene_deg`` applies across contrasts; it
    is applied here within one contrast so two source identifiers that alias to
    one axis gene cannot emit two different statistics for the same cell of the
    published table.
    """
    aligned = pd.DataFrame(aligned)
    if aligned.empty:
        raise ValueError('no aligned rows')
    work = aligned.assign(effect=deg_effect(aligned.log2fc, aligned.padj))
    work = work.sort_values(['gene_index', 'contrast', 'effect'],
                            ascending=[True, True, False], kind='stable')
    return work.drop_duplicates(['gene_index', 'contrast'], keep='first')


def build_dataset_products(local_deg, resolver, axis_symbols, dataset_id,
                           axis_size=AXIS_SIZE):
    """Build one dataset's compact table, full-axis rank and statistics.

    ``local_deg`` is an item-local per-contrast table (``LOCAL_DEG_COLUMNS``);
    ``resolver``/``axis_symbols`` bind the frozen axis. Returns
    ``(table, rank, stats, alignment_report)`` where ``table`` follows
    ``DEG_COLUMNS`` sorted by rank and ``rank`` follows ``RANK_COLUMNS`` over the
    whole frozen axis (unmeasured genes at ``n_measured + 1``).
    """
    aligned, report = align_source_table(local_deg, resolver)
    stats_rows = winning_contrast_rows(aligned)
    table, rank, stats = dataset_deg_rank(
        aligned[['gene_index', 'log2fc', 'padj', 'contrast']], axis_size=axis_size)
    keep = stats_rows[['gene_index', 'contrast', 'stat', 'p_value', 'base_mean']]
    if keep.duplicated(['gene_index', 'contrast']).any():
        raise ValueError(f'{dataset_id}: duplicate winning contrast rows')
    table = table.merge(keep, on=['gene_index', 'contrast'], how='left',
                        validate='one_to_one')
    symbols = pd.Series(np.asarray(axis_symbols, dtype=object),
                        index=np.arange(axis_size, dtype=np.int64))
    table['gene_symbol'] = table.gene_index.map(symbols)
    table = table[list(DEG_COLUMNS)]
    table = table.sort_values(['rank', 'gene_index'], kind='stable').reset_index(drop=True)
    stats = dict(stats)
    stats.update({'dataset_id': dataset_id, 'alignment': report,
                  'n_contrasts_published': int(table.contrast.nunique())})
    return table, rank, stats, report


def validate_table(table, axis_size=AXIS_SIZE):
    """Validate a published compact DEG table against the ticket-01 contract."""
    table = pd.DataFrame(table)
    if tuple(table.columns) != DEG_COLUMNS:
        raise ValueError(f'DEG table columns differ from the contract: {list(table.columns)}')
    if table.empty:
        raise ValueError('DEG table is empty')
    if not table.gene_index.is_unique:
        raise ValueError('DEG table repeats a gene_index')
    if table.gene_index.min() < 0 or table.gene_index.max() >= axis_size:
        raise ValueError('DEG table gene_index outside the frozen axis')
    if not table['rank'].tolist() == list(range(1, len(table) + 1)):
        raise ValueError('DEG table is not sorted by the within-dataset rank')
    if not np.allclose(table.rank_norm.to_numpy(), table['rank'].to_numpy() / len(table)):
        raise ValueError('DEG table rank_norm is not rank / n_measured')
    order = np.lexsort((table.gene_index.to_numpy(dtype=np.int64),
                        -table.effect.to_numpy(dtype=np.float64)))
    if table.gene_index.to_numpy(dtype=np.int64).tolist() != \
            table.gene_index.to_numpy(dtype=np.int64)[order].tolist():
        raise ValueError('DEG table is not ordered by effect with gene_index tie-break')
    return table


def validate_rank(rank, axis_size=AXIS_SIZE):
    """Validate a published full-axis rank frame."""
    rank = pd.DataFrame(rank)
    if tuple(rank.columns) != RANK_COLUMNS:
        raise ValueError(f'rank frame columns differ from the contract: {list(rank.columns)}')
    if len(rank) != axis_size:
        raise ValueError(f'rank frame must hold {axis_size} genes, found {len(rank)}')
    if rank.gene_index.tolist() != list(range(axis_size)):
        raise ValueError('rank frame is not the frozen 0..18163 order')
    measured = rank['rank'] <= rank['rank'].max()
    if (rank['rank'] < 1).any():
        raise ValueError('rank frame has a rank below 1')
    n_measured = int((rank['rank'] != rank['rank'].max()).sum())
    return rank


def genotype_split(genotypes):
    """Split Reichart donors by ``Primary.Genetic.Diagnosis`` (ADR 0066).

    ``genotypes`` is a frame with ``donor_id`` and ``genotype``. The split is
    genotype-first, never the study's disease label: a donor whose genotype
    declares a pathogenic variant is *not* part of ``nHF`` regardless of the
    clinical label, the single NCCM donor is not used at all, and normal donors
    supply the ``nHF`` controls. Returns a frame with ``donor_id``,
    ``genotype``, ``genotype_status`` and ``nHF_role`` (``case``, ``control`` or
    ``excluded``/``export_only``), plus a counts summary.
    """
    frame = pd.DataFrame(genotypes)
    if not {'donor_id', 'genotype'} <= set(frame.columns):
        raise ValueError('genotype table requires donor_id and genotype')
    if frame.donor_id.duplicated().any():
        raise ValueError('genotype table repeats a donor')
    roles, statuses = [], []
    for genotype in frame.genotype.astype(str):
        key = genotype.strip().lower()
        if key in ('normal', 'donor', 'healthy', 'non-dcm', 'control'):
            statuses.append('negative')
            roles.append('control')
        elif 'arvc' in key:
            statuses.append('pathogenic')
            roles.append('export_only')
        elif 'nccm' in key:
            statuses.append('pathogenic')
            roles.append('excluded')
        elif key.startswith('dcm') or 'dcm' in key:
            if any(token in key for token in ('pathogenic', 'positive', 'variant',
                                              'mutation', 'carrier')):
                statuses.append('pathogenic')
                roles.append('export_only')
            else:
                statuses.append('negative')
                roles.append('case')
        else:
            raise ValueError(f'unrecognised genotype label: {genotype!r}')
    result = frame.assign(genotype_status=statuses, nHF_role=roles)
    summary = result.groupby(['genotype_status', 'nHF_role'], sort=True).size()
    summary = {f'{status}:{role}': int(count)
               for (status, role), count in summary.items()}
    return result, summary


def contrast_selection(rows, dataset_id, contrasts):
    """Select one dataset's declared contrasts out of one joint-design DEG table.

    A shared-control source is fitted once (``~ arm`` / ``~ drug``); each
    published dataset is a declared selection of contrast labels from that single
    fit. Selecting an undeclared label or producing an empty selection is an
    error, so the joint design is recorded, not re-fitted per dataset.
    """
    rows = pd.DataFrame(rows)
    if 'contrast' not in rows:
        raise ValueError('a joint-design DEG table requires a contrast column')
    wanted = list(contrasts)
    absent = sorted(set(wanted) - set(rows.contrast.astype(str)))
    if absent:
        raise ValueError(f'{dataset_id}: contrasts not produced by the joint design: {absent}')
    selected = rows[rows.contrast.astype(str).isin(wanted)].copy()
    if selected.empty:
        raise ValueError(f'{dataset_id}: contrast selection is empty')
    return selected


# The single live DEG worker (ADR 0066 D9). The `arc` conda env's R wrapper is
# mis-patched on this workstation, so the R entry point is invoked as
# R_HOME=<env>/lib/R <env>/lib/R/bin/exec/R (the documented pattern for the
# repository's pinned R envs). The path is configuration, not science: item
# scripts pass it and record only the env name/version, never a host path.
DEFAULT_R_HOME = '/home/chang/miniconda3/envs/arc/lib/R'
DEG_WORKER = Path(__file__).with_name('disease_deg.R')
DEG_SUMMARY_PREFIX = 'DEG_SUMMARY '


def r_entry_point(r_home=None):
    """Return the R executable of the pinned `arc` env (see module docstring)."""
    home = Path(r_home or os.environ.get('ARC_R_HOME', DEFAULT_R_HOME))
    executable = home / 'bin' / 'exec' / 'R'
    if not executable.is_file():
        raise FileNotFoundError(f'pinned R entry point not found: {executable}')
    return executable, home


def run_deg_worker(matrix, samples, design, out_csv, r_home=None, timeout=7200):
    """Fit one joint design through ``disease_deg.R`` and return its summary.

    ``matrix``/``samples``/``design``/``out_csv`` are paths; the design is the
    JSON-ready dict documented in ``disease_deg.R``. The R worker writes the
    long per-contrast table and prints its fit description as the last
    ``DEG_SUMMARY <json>`` line, which is returned (never parsed for science).
    """
    import json
    import subprocess
    import tempfile

    executable, home = r_entry_point(r_home)
    with tempfile.TemporaryDirectory() as scratch:
        design_path = Path(scratch) / 'design.json'
        design_path.write_text(json.dumps(design, indent=2, sort_keys=True) + '\n')
        command = [str(executable), '--no-save', '--slave', '-f', str(DEG_WORKER),
                   '--args', '--matrix', str(matrix), '--samples', str(samples),
                   '--design', str(design_path), '--out', str(out_csv)]
        env = dict(os.environ)
        env['R_HOME'] = str(home)
        env['OMP_NUM_THREADS'] = env.get('OMP_NUM_THREADS', '2')
        env['MKL_NUM_THREADS'] = env.get('MKL_NUM_THREADS', '2')
        env['OPENBLAS_NUM_THREADS'] = env.get('OPENBLAS_NUM_THREADS', '2')
        completed = subprocess.run(command, capture_output=True, text=True,
                                   env=env, timeout=timeout)
    if completed.returncode != 0:
        raise RuntimeError('disease_deg.R failed:\n' + completed.stderr[-4000:])
    for line in reversed(completed.stdout.splitlines()):
        if line.startswith(DEG_SUMMARY_PREFIX):
            return json.loads(line[len(DEG_SUMMARY_PREFIX):])
    raise RuntimeError('disease_deg.R produced no DEG_SUMMARY line')


def category_root(output_root, directory):
    """The accepted product directory of one disease category (\`\`nhf\`\`, \`\`icm\`\`, ...)."""
    return Path(output_root) / directory


def catalog_path(output_root, category):
    return category_root(output_root, category) / 'catalog.json'


def manifest_path(output_root, category):
    return category_root(output_root, category) / 'manifest.json'


def catalog_digest_bytes(catalog):
    """Canonical catalog bytes without the self-referential ``catalog_sha256``.

    The rule (shared with the #189 producer) is: serialise the catalog without
    the field using ``json.dumps(..., indent=2, sort_keys=True) + "\\n"`` and
    sha256 those bytes; that digest is stored in ``catalog_sha256`` and repeated
    in the category manifest. Verification pops the field, re-serialises the same
    way and compares, so the digest is a function of the catalog content only.
    """
    body = dict(catalog)
    body.pop('catalog_sha256', None)
    return (json.dumps(body, indent=2, sort_keys=True) + '\n').encode()


def catalog_digest(catalog):
    return hashlib.sha256(catalog_digest_bytes(catalog)).hexdigest()


def write_catalog(path, catalog):
    """Write a catalog with its canonical ``catalog_sha256``; return the digest."""
    digest = catalog_digest(catalog)
    body = dict(catalog)
    body['catalog_sha256'] = digest
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(body, indent=2, sort_keys=True) + '\n')
    return digest


def read_catalog(path):
    """Read a catalog and verify its embedded ``catalog_sha256``."""
    path = Path(path)
    catalog = json.loads(path.read_text())
    declared = catalog.get('catalog_sha256')
    if not isinstance(declared, str) or len(declared) != 64:
        raise ValueError(f'{path}: catalog has no canonical catalog_sha256')
    if declared != catalog_digest(catalog):
        raise ValueError(f'{path}: catalog_sha256 does not match the catalog content')
    return catalog


def build_catalog(category, directory, tier, evidence_type, axis, axis_sha256,
                  source_records, rejected=(), axis_size=AXIS_SIZE):
    """Assemble a category catalog from the declared source/eligibility records.

    Canonical #190 shape: ``datasets`` is a **list of full record objects** (no
    separate id list / record map), and the content digest is carried as
    ``catalog_sha256``.
    """
    seen = set()
    records = []
    for record in source_records:
        dataset_id = record['dataset_id']
        if dataset_id in seen:
            raise ValueError(f'{category}: duplicate dataset id {dataset_id}')
        seen.add(dataset_id)
        records.append(record)
    catalog = {
        'schema': CATALOG_SCHEMA,
        'category': category,
        'directory': directory,
        'tier': int(tier),
        'evidence_type': evidence_type,
        'axis_size': int(axis_size),
        'axis': {'path': 'stage2/01_bases/output/gene-identity/genes.csv',
                 'sha256': axis_sha256, 'size': int(axis_size)},
        'rank_contract': RANK_SCHEMA,
        'datasets': records,
        'rejected': list(rejected),
    }
    return catalog


def assemble_manifest(category, recipe, inputs, outputs, counts, catalog_sha256):
    """Assemble the per-category release manifest (ADR 0062 shape)."""
    manifest = {
        'schema': RELEASE_SCHEMA,
        'category': category,
        'recipe': recipe,
        'axis': {'path': 'stage2/01_bases/output/gene-identity/genes.csv'},
        'inputs': dict(sorted(inputs.items())),
        'outputs': dict(sorted(outputs.items())),
        'counts': dict(sorted(counts.items())),
        'catalog_sha256': catalog_sha256,
    }
    return manifest


def _sha256(path):
    import hashlib
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


class DiseaseDeg:
    """Reader for the published per-category disease-DEG products.

    Every read is digest-checked against the category ``manifest.json``, so a
    later ticket can trust that the table it fuses is the accepted one. The
    reader is importable without the external item trees: it only reads
    ``output/disease-deg/``.
    """

    def __init__(self, root):
        self.root = Path(root)

    # ---- structure -----------------------------------------------------
    def categories(self):
        """Product directory names (``af``, ``icm``, ``nhf``), sorted."""
        found = []
        for path in sorted(self.root.iterdir()):
            if path.is_dir() and (path / 'catalog.json').is_file():
                found.append(path.name)
        return found

    def directory(self, key):
        """Resolve a logical category (``nHF``) or a product directory (``nhf``)."""
        if (self.root / key / 'catalog.json').is_file():
            return key
        for name in self.categories():
            catalog = json.loads((self.root / name / 'catalog.json').read_text())
            if catalog.get('category') == key:
                return name
        raise KeyError(f'unknown disease category {key!r}')

    def catalog(self, category):
        path = catalog_path(self.root, self.directory(category))
        if not path.is_file():
            raise FileNotFoundError(f'no disease-DEG catalog for category {category!r}')
        catalog = read_catalog(path)
        if catalog.get('schema') != CATALOG_SCHEMA:
            raise ValueError(f'{path}: unexpected catalog schema')
        if catalog.get('category') != category and catalog.get('directory') != category:
            raise ValueError(f'{path}: catalog category mismatch')
        return catalog

    def manifest(self, category):
        path = manifest_path(self.root, self.directory(category))
        if not path.is_file():
            raise FileNotFoundError(f'no disease-DEG manifest for category {category!r}')
        manifest = json.loads(path.read_text())
        if manifest.get('schema') != RELEASE_SCHEMA:
            raise ValueError(f'{path}: unexpected release schema')
        return manifest

    def datasets(self, category):
        """Dataset ids of a category, in catalog order."""
        return [record['dataset_id'] for record in self.catalog(category)['datasets']]

    def record(self, category, dataset_id):
        catalog = self.catalog(category)
        for record in catalog['datasets']:
            if record['dataset_id'] == dataset_id:
                return record
        raise KeyError(f'{category}: unknown dataset {dataset_id!r}')

    def eligibility(self, category, dataset_id):
        """The declared eligibility/shared-control fields of one dataset."""
        record = self.record(category, dataset_id)
        keys = ('category', 'source_item', 'species', 'model_system', 'evidence_type',
                'tier', 'genotype_status', 'time_window_basis', 'exploratory',
                'shared_control_group', 'n_case', 'n_control', 'method')
        return {key: record[key] for key in keys if key in record}

    # ---- payloads ------------------------------------------------------
    def _verified_path(self, category, name):
        manifest = self.manifest(category)
        declared = manifest['outputs'].get(name)
        if not isinstance(declared, dict):
            raise ValueError(f'{category}: {name} is not a declared output')
        path = category_root(self.root, self.directory(category)) / name
        if not path.is_file():
            raise FileNotFoundError(f'{category}: missing declared output {name}')
        if _sha256(path) != declared['sha256']:
            raise ValueError(f'{category}: {name} digest differs from the manifest')
        return path

    def table(self, category, dataset_id):
        record = self.record(category, dataset_id)
        table = pd.read_csv(self._verified_path(category, record['deg']),
                            float_precision='round_trip')
        return validate_table(table)

    def rank(self, category, dataset_id):
        record = self.record(category, dataset_id)
        rank = pd.read_csv(self._verified_path(category, record['rank']),
                           dtype={'gene_index': np.int64},
                           float_precision='round_trip')
        return validate_rank(rank)

    def digest(self, category, name):
        """Declared sha256 of one accepted output."""
        return self.manifest(category)['outputs'][name]['sha256']
