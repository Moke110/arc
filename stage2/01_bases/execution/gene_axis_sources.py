#!/usr/bin/env python3
"""``pert_data/`` cardiotoxicity adapters for the gene-axis evidence product.

Two already-owned external items supply anthracycline full-transcriptome
evidence; neither is re-downloaded and neither is modified:

* ``pert_data/DToxS-iPSC/`` (GSE217424 SuperSeries, GSE217421 iPSC-CM): the
  upstream per-comparison ``Calc.CTRL-<CODE>.tsv`` tables already carry
  ``logFC``, ``PValue`` and ``FDR``. The five cardiotoxic positives present are
  DAU/DOX/EPI/IDA/TRA; there is no MTX table.
* ``pert_data/GSE244740/``: ScreenSeq hiPSC-CM raw counts plus the series
  matrix. Doxorubicin, Idarubicin and Mitoxantrone wells are compared with the
  DMSO vehicle.

Both adapters return a tidy per-contrast table
``contrast,gene_symbol,log2fc,pvalue,padj,n_control,n_treated``. Collapsing the
contrasts into one source (per gene **maximum** effect score) happens in
``build_gene_axis.py`` through ``bases.gene_axis.per_gene_deg`` so the rule is
shared with every other DEG source.

GSE244740 has no ship-with DEG table, so it is computed here with a declared,
deterministic method: library-size normalisation over **all** source genes of
the selected wells, ``log2(CPM + 1)``, a two-sample Welch t-test per gene and
Benjamini-Hochberg adjustment within each drug contrast. No moderation, no
imputation and no gene filter; the sample set is fixed by the series matrix.
"""
from pathlib import Path
import csv
import gzip
import io
import tarfile

import numpy as np
import pandas as pd
from scipy import stats

DTOXS_CODES = ('DAU', 'DOX', 'EPI', 'IDA', 'TRA')
GSE244740_DRUGS = ('Doxorubicin', 'Idarubicin', 'Mitoxantrone')
GSE244740_VEHICLE = 'DMSO'
GSE244740_TIME = '24h'
GFRAME_COLUMNS = ['contrast', 'gene_symbol', 'log2fc', 'pvalue', 'padj',
                  'n_control', 'n_treated']


def _write_cache(frame, cache_dir, name):
    if cache_dir is None:
        return None
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / name
    frame.to_csv(path, index=False, float_format='%.10g', lineterminator='\n')
    return path


def dtoxs_contrasts(tar_path, codes=DTOXS_CODES, cache_dir=None):
    """Read the DToxS per-comparison DEG tables for the cardiotoxic positives.

    Returns ``(frame, provenance)``. ``contrast`` is the upstream table stem so
    the cell line, hour and plate that produced each contrast stay visible.
    """
    tar_path = Path(tar_path)
    wanted = {f'CTRL-{code}.' for code in codes}
    frames = []
    wells = set()
    provenance = {'tar': tar_path.name, 'codes': list(codes), 'tables': []}
    with tarfile.open(tar_path, 'r|gz') as archive:
        for member in archive:
            name = Path(member.name).name
            if not name.startswith('Human-') or 'Calc.CTRL-' not in name:
                continue
            if not any(token in name for token in wanted) or not name.endswith('.tsv'):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ValueError(f'unreadable DToxS member: {member.name}')
            payload = handle.read()
            # 32 data fields: gene symbol, N raw control + M raw treated, the
            # same wells normalised, logFC, logCPM, PValue, FDR, Regulation.
            table = pd.read_csv(io.BytesIO(payload), sep='\t', header=None,
                                skiprows=1, dtype=str, keep_default_na=False)
            width = table.shape[1]
            if width < 6:
                raise ValueError(f'malformed DToxS table: {member.name}')
            table = table.iloc[:, [0, width - 5, width - 4, width - 3, width - 2]]
            table.columns = ['gene_symbol', 'log2fc', 'logcpm', 'pvalue', 'padj']
            header = io.BytesIO(payload).readline().decode().rstrip('\n').split('\t')
            raw = [column for column in header if not column.endswith('.Norm')
                   and column not in ('logFC', 'logCPM', 'PValue', 'FDR', 'Regulation')]
            wells.update(raw)
            treated = sum(1 for column in raw if column.split('.')[0] in codes)
            control = len(raw) - treated
            contrast = Path(member.name).stem
            table = table.assign(contrast=contrast,
                                 n_control=control, n_treated=treated)
            for column in ('log2fc', 'pvalue', 'padj'):
                table[column] = pd.to_numeric(table[column], errors='coerce')
            frames.append(table[GFRAME_COLUMNS])
            provenance['tables'].append({'member': member.name, 'rows': len(table),
                                         'n_control': control, 'n_treated': treated})
    if not frames:
        raise ValueError('no DToxS cardiotoxic table found in the archive')
    frame = pd.concat(frames, ignore_index=True)
    frame = frame.sort_values(['contrast', 'gene_symbol'], kind='stable').reset_index(drop=True)
    provenance['rows'] = len(frame)
    provenance['n_distinct_wells'] = len(wells)
    provenance['contrasts'] = sorted(frame.contrast.unique().tolist())
    _write_cache(frame, cache_dir, 'dtoxs_contrasts.csv')
    return frame, provenance


def _plate(well):
    """Plate prefix of a ScreenSeq well id (``S4881Nr8_A1`` -> ``S4881Nr8``)."""
    return str(well).split('_')[0]


def _series_samples(series_path):
    """well id -> {treatment, time} from the GSE244740 series matrix.

    The series matrix lists characteristics before the description, so the rows
    are collected by key first and paired by column position afterwards.
    """
    rows = {}
    with gzip.open(series_path, 'rt') as handle:
        for line in handle:
            if not line.startswith('!Sample_'):
                continue
            fields = [value.strip().strip('"') for value in line.rstrip('\n').split('\t')]
            key, values = fields[0], fields[1:]
            if key == '!Sample_description':
                rows['well'] = values
            elif key == '!Sample_characteristics_ch1' and values:
                if values[0].startswith('treatment:'):
                    rows['treatment'] = [value.split(':', 1)[1].strip() for value in values]
                elif values[0].startswith('time:'):
                    rows['time'] = [value.split(':', 1)[1].strip() for value in values]
    if 'well' not in rows:
        raise ValueError('GSE244740 series matrix has no sample description row')
    samples = {}
    for position, well in enumerate(rows['well']):
        samples[well] = {'well': well,
                         'treatment': rows.get('treatment', [''] * len(rows['well']))[position],
                         'time': rows.get('time', [''] * len(rows['well']))[position]}
    return samples


def _benjamini_hochberg(pvalues):
    pvalues = np.asarray(pvalues, dtype=np.float64)
    n = len(pvalues)
    order = np.argsort(pvalues, kind='stable')
    ranked = pvalues[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty(n, dtype=np.float64)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def _welch(counts, treated_mask, control_mask, library):
    """Per-gene Welch t-test and log2 fold change on log2(CPM+1) values.

    ``library`` is the per-well total of **all** source genes (not only the
    axis-mapped ones), as declared in the method.
    """
    linear = counts / library[None, :] * 1e6
    values = np.log2(linear + 1.0)
    treated, control = values[:, treated_mask], values[:, control_mask]
    n_t, n_c = treated.shape[1], control.shape[1]
    mean_t, mean_c = treated.mean(axis=1), control.mean(axis=1)
    var_t = treated.var(axis=1, ddof=1) if n_t > 1 else np.zeros_like(mean_t)
    var_c = control.var(axis=1, ddof=1) if n_c > 1 else np.zeros_like(mean_c)
    log2fc = mean_t - mean_c
    se_t, se_c = var_t / n_t, var_c / n_c
    total = se_t + se_c
    with np.errstate(invalid='ignore', divide='ignore'):
        statistic = np.where(total > 0, log2fc / np.sqrt(total), 0.0)
        degrees = np.where(total > 0,
                           total ** 2 / (se_t ** 2 / max(n_t - 1, 1) + se_c ** 2 / max(n_c - 1, 1)),
                           np.inf)
        pvalue = np.where(total > 0, 2.0 * stats.t.sf(np.abs(statistic), degrees), 1.0)
    return log2fc, pvalue, n_c, n_t


def gse244740_contrasts(counts_path, series_path, ensembl_symbol, drugs=GSE244740_DRUGS,
                        vehicle=GSE244740_VEHICLE, time_point=GSE244740_TIME, cache_dir=None):
    """Derive per-drug DEG contrasts for the GSE244740 ScreenSeq panel.

    ``ensembl_symbol`` maps a version-stripped Ensembl gene id to the frozen-axis
    symbol; genes without a mapping are dropped from the contrast table but their
    counts still enter the library size (the normalisation is over all source
    genes). Returns ``(frame, provenance)``.
    """
    counts_path, series_path = Path(counts_path), Path(series_path)
    samples = _series_samples(series_path)
    # Treatment labels look like ``Doxorubicin, Con5`` / ``DMSO, Con0``: the
    # first field is the compound, the rest is the dose index. All doses of a
    # compound are pooled into that compound's contrast (declared). Controls are
    # the **same-plate** DMSO wells: the three anthracyclines sit on one plate
    # (S4881Nr8 here) whose DMSO wells are the only batch-matched vehicle.
    treated = {}
    for well, info in samples.items():
        if info.get('time') != time_point:
            continue
        treatment = info.get('treatment', '').split(',')[0].strip()
        if treatment in drugs:
            treated[well] = treatment
    if not treated:
        raise ValueError('GSE244740 series matrix yields no selected drug wells')
    plates = sorted({_plate(well) for well in treated})
    controls = {}
    for well, info in samples.items():
        if info.get('time') != time_point or _plate(well) not in plates:
            continue
        if info.get('treatment', '').split(',')[0].strip() == vehicle:
            controls[well] = vehicle
    if not controls:
        raise ValueError('GSE244740 has no same-plate vehicle well')
    selected = {**controls, **treated}
    wells = sorted(selected)
    labels = np.asarray([selected[well] for well in wells])

    header = pd.read_csv(counts_path, sep='\t', nrows=0, compression='gzip')
    ensembl_column = header.columns[0]
    usecols = [ensembl_column] + wells
    chunks = []
    library = np.zeros(len(wells), dtype=np.float64)
    n_unmapped = 0
    reader = pd.read_csv(counts_path, sep='\t', usecols=usecols, chunksize=4000,
                         compression='gzip', dtype={ensembl_column: str})
    for chunk in reader:
        identifiers = chunk[ensembl_column].astype(str).str.split('.').str[0]
        values = chunk[wells].to_numpy(dtype=np.float32)
        library += values.sum(axis=0, dtype=np.float64)
        symbols = identifiers.map(ensembl_symbol)
        mapped = symbols.notna().to_numpy()
        n_unmapped += int((~mapped).sum())
        if mapped.any():
            chunks.append(pd.DataFrame(values[mapped], columns=wells).assign(
                gene_symbol=symbols[mapped].to_numpy()))
    if not chunks:
        raise ValueError('GSE244740 has no axis-mapped gene')
    matrix = pd.concat(chunks, ignore_index=True)
    genes = matrix.pop('gene_symbol').to_numpy()
    values = matrix[wells].to_numpy(dtype=np.float64)
    if (library <= 0).any():
        raise ValueError('GSE244740 has an empty library')

    frames = []
    provenance = {'counts': counts_path.name, 'series': series_path.name,
                  'n_wells': int(len(wells)), 'n_genes_read': int(len(genes)),
                  'plates': plates, 'n_vehicle_wells': int((labels == vehicle).sum()),
                  'n_source_genes_unmapped': int(n_unmapped),
                  'vehicle': vehicle, 'time': time_point, 'contrasts': []}
    for drug in drugs:
        treated_mask = labels == drug
        control_mask = labels == vehicle
        if not treated_mask.any():
            continue
        n_t, n_c = int(treated_mask.sum()), int(control_mask.sum())
        if n_t < 2 or n_c < 2:
            raise ValueError('GSE244740 contrast needs at least two wells per arm')
        log2fc, pvalue, n_c, n_t = _welch(values, treated_mask, control_mask, library)
        padj = _benjamini_hochberg(pvalue)
        frames.append(pd.DataFrame({
            'contrast': drug, 'gene_symbol': genes, 'log2fc': log2fc,
            'pvalue': pvalue, 'padj': padj, 'n_control': n_c, 'n_treated': n_t}))
        provenance['contrasts'].append({'drug': drug, 'n_treated': int(n_t),
                                        'n_control': int(n_c)})
    frame = pd.concat(frames, ignore_index=True)
    frame = frame.sort_values(['contrast', 'gene_symbol'], kind='stable').reset_index(drop=True)
    provenance['rows'] = len(frame)
    _write_cache(frame, cache_dir, 'gse244740_contrasts.csv')
    return frame, provenance
