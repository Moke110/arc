#!/usr/bin/env python3
"""Write the human-facing gene-axis rule document from the released products.

The document is a reading aid, not a release product. It opens with the ADR 0071
**v2** build rule (DEG definition, tiered source collection and its stop test,
``DEG2000``, ``f2000`` and the tenv-HVG fallback), so the released numbers cannot
be read under the retired ``pN`` rule; it then carries the source dataset table
(one row per published dataset plus the tenv-HVG fusion source), the per-disease
table (source-dataset count, DEG union, ``DEG2000``), the coverage table and the
released evidence figure.

Its inputs are the published coverage CSV (``output/gene-axis/
source_topn_coverage.csv``), the released PNG, the declared published
``output/disease-deg/<category>/catalog.json`` files, the published
``output/disease-fusion/catalog.json`` (the per-disease level-1 record), the
frozen axis identity table the catalogs declare (``axis.path``), the release
``manifest.json`` (whose ``release.algorithm`` marker refuses a pre-v2 release)
and the tracked ``gene_axis_source_metadata`` mapping. Only the rule text lives here as a
constant, and it mirrors ``bases.gene_axis.FILL_RULE``.

    python stage2/01_bases/execution/make_gene_axis_docx.py \
        --coverage stage2/01_bases/output/gene-axis/source_topn_coverage.csv \
        --figure stage2/01_bases/figure/gene_axis_rank_sum_histogram.png \
        --out stage2/01_bases/output/gene-axis/gene_axis_fill_rule.docx

Every path defaults to the repo layout and can be overridden: ``--catalog
<category>=<path>`` (repeatable), ``--fusion-catalog PATH``, ``--manifest PATH``,
``--metadata PATH`` (imported from that file dynamically) and ``--axis PATH``.
The DOCX container is re-zipped with fixed entry timestamps, so two runs from the
same inputs and the same ``python-docx`` produce identical bytes.
"""
import argparse
import importlib.util
import io
import json
import zipfile
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Cm, Pt, RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

import public_api  # noqa: F401  (binds the single public implementation in output/api)
import gene_axis_source_metadata as source_metadata
from gene_axis_source_metadata import TENV_TREATED_WELLS, TENV_VEHICLE_WELLS
from bases import gene_axis

TITLE = 'Stage 2 基因轴取用规则与来源覆盖'
#: The published coverage table's long-format contract, imported from the module
#: that owns it so the document cannot drift from the release.
COVERAGE_COLUMNS = gene_axis.COVERAGE_COLUMNS
#: The columns rendered in the document's coverage table (``note`` is long prose
#: that the section's own paragraph carries instead of repeating per row).
COVERAGE_TABLE_COLUMNS = ('N', 'group_label', 'numerator', 'denominator', 'rate')
SOURCE_TABLE_COLUMNS = ('category', 'accession', 'sample type', 'groups',
                        'sample sizes', 'sequencing technology', 'CM DEG count')
SOURCE_TABLE_CAPTION = '表 1　已发布来源数据集信息'
SOURCE_TABLE_HEADING = '2 来源数据集信息'
DISEASE_TABLE_COLUMNS = ('disease', 'source datasets', 'tier path',
                         'DEG union', 'DEG2000', 'eligible', 'note')
DISEASE_TABLE_CAPTION = '表 2　每个疾病的来源数据集数、DEG 并集与 DEG2000'
DISEASE_TABLE_HEADING = '3 每个疾病的来源数与 DEG 资格'
#: A ``DEG2000`` this small is named gene by gene rather than described as a
#: count, so a single-gene result cannot look like a rounding of a larger one.
NAMED_DEG2000_LIMIT = 5

RULE_HEADING = '1 取用规则 ADR 0072'

DEG_DEFINITION = (
    'DEG 定义：一个基因在某数据集里是 DEG，当它在至少一个已声明对比中满足 '
    'padj < 0.05 且 |log2FC| ≥ 2（FC ≥ 4）。任一对比通过即算（pass-any-contrast），'
    '不是“最佳对比”塌缩出的单一胜负集；padj 为 NA 视为不显著；效应方向被记录，'
    '但不参与筛选。本文件的 CM DEG 数就是这一口径下的基因数，'
    '它与 padj<0.05 的 padj-only 显著基因数不是同一个量。'
)

RULE_STEPS = (
    '第一步（来源内秩）：每个来源只在自己的 DEG 内排序，秩取 1..n_DEG，按效应 '
    '|log2FC|·(−log10 padj) 降序、gene_index 破平；非 DEG 一律记 n_DEG+1；'
    'rank_norm = rank / n_DEG。分母是 DEG 数，绝不是测得基因数。'
    '旧的 pN 与“测得基因秩”随 ADR 0071 退役，本文件不再使用。',
    '第二步（每病融合，level 1）：disease_rank_sum = Σ（该病各来源的 rank_norm）；'
    '疾病秩是该和值的升序位置，gene_index 破平。DEG2000_d = 该病 DEG 并集中疾病秩 '
    '≤ 2,000 的基因；并集小于 2,000 时 DEG2000_d 就是整个并集。'
    'DEG 并集的大小不受约束，也不做任何填充。',
    '第三步（分层收集）：按 tier 升序收取来源。tier 1 = 人房/室组织且 CM 特异'
    '（单细胞/单核 CM pseudobulk，或 CM 核分选 bulk）；tier 2 = hiPSC/hESC-CM、'
    '类器官或工程化心脏组织（bulk），其中临床类别（nHF/ICM/AF）在本层只收非遗传模型；'
    'tier 3 = 人房/室 bulk 全组织（非 CM 特异）。每收完一层，用解析到冻结 18,164 轴的'
    '累计 DEG 并集与 2,000 比较：达到即停，不再进入更低层；未声明的层从阶梯里缺席'
    '（AF 没有 tier 2，TOX 从 tier 2 进入）。',
    '第四步（疾病资格）：每病完成来源分析、level 1 和 DEG2000 后，冻结轴 DEG 并集至少 100 个才进入跨病融合。AF 的一个基因及上游证据保留，但本次不参与 rank_sum、f2000 或疾病驱动选轴。TOX 每个数据集先合并独立证据支持的药物暴露与匹配溶媒对比，一份数据集提供一个来源秩；机制仅作注释。',
    '第五步（轴填充，level 2）：对 N ∈ {2,000, 5,000}，前半 0.5N 个基因取 '
    'tenv-HVG 秩 ≤ N/2。后半 0.5N 个基因取自剩余候选中 f2000 ≥ 1 的基因，按 '
    'f2000 降序 → rank_sum 升序 → gene_index 升序。'
    'f2000 = 合格疾病里 DEG2000_d 含该基因的疾病数（tenv 不是疾病项，不参与）；'
    'rank_sum = Σ_合格疾病d rank_d / n_d，基因在该病 DEG2000_d 之外时 rank_d 记 n_d+1，'
    '绝不省略。',
    '第六步（不足时的补齐）：若 f2000 ≥ 1 的剩余候选不足 0.5N，缺口用尚未入选的 '
    'tenv-HVG 基因按 tenv-HVG 秩补齐；f2000 = 0 的基因绝不进入疾病那半。',
    '合并两半后按 gene_index 升序发布；轴是冻结 18,164 轴的有序子集，保留原始 '
    'gene_index 整数，绝不重编号为 0..N-1。',
)

F2000_DEFINITION = (
    '合格疾病中 DEG2000_d 含该基因的疾病数量。本次五个声明疾病有四个合格，取值 0..4。'
    'tenv-HVG 不计入；AF 保留证据但不计入。f2000 = 0 的基因只能由 tenv-HVG 分配或补齐进入轴。'
)
CLARIFICATIONS = (
    ('pN 已退役', 'f2000 只计合格疾病的 DEG2000，不计 tenv-HVG。'),
    ('排序方向', '先 f2000 降序，再 rank_sum 升序，最后 gene_index 升序。'),
    ('小集合与资格不同', 'DEG2000 不填充。少于 2000 不代表失败；少于 100 个轴内 DEG 才不参与跨病融合。'),
    ('覆盖率分母', '合格疾病各用自己的 DEG2000 大小；tenv-HVG 用 N。排除疾病保留空值和原因，不能当作零覆盖或成功覆盖。'),
)
SOURCE_TABLE_FOOTNOTE = (
    'CM DEG 数来自各 catalog 的 n_deg，是 pass-any-contrast DEG 数，不是 padj-only 显著基因数。'
    '样本量口径：非 TOX 来源为已发布病例和对照数；TOX 是已声明对比的 pooled wells，'
    '有多个对比时显示范围，n_case_basis 与 n_control_basis 记录计数口径，不能作为全系列独立供体总数。TOX 表中药物仅指独立心脏毒性证据准入药物；'
    '保留匹配溶媒、供体阻断及原有细胞系、板和时间分层。没有要求每个剂量和时间都出现毒性。'
)
DISEASE_TABLE_FOOTNOTE = (
    'DEG union 是全部来源解析到冻结 18164 轴的去重 DEG 并集。DEG2000 是该并集中排名前 2000 的基因，'
    '不足时不填充。eligible 来自疾病资格记录；只有合格疾病才贡献 level 2，AF 的一个基因保留不参与。'
    '各数值来自 disease-fusion/catalog.json 并与 disease-deg catalogs 交叉校验。'
)

#: The tenv-HVG fusion row is not a case/control contrast. Its sample sizes are
#: the published corpus populations pinned in ``output/view-2k/manifest.json``
#: ``sources.tenv`` (treated = 4,206, vehicle = 272; also
#: ``log/raw_counts_provenance.json``: 4,510 native wells kept). The two numbers
#: are declared once, in ``gene_axis_source_metadata``, and imported here.
TENV_HVG_ROW = {
    'category': 'tenv-HVG (fusion source)',
    'sample type': 'hiPSC-CM 96-well plate lysates (TempO-Seq corpus)',
    'groups': 'not a contrast: ranked by the published per-well dispersion rank',
    'sample sizes': f'{TENV_TREATED_WELLS:,} / {TENV_VEHICLE_WELLS:,} (treated / vehicle wells)',
    'CM DEG count': 'not a DEG source',
}

REPRODUCE = (
    'python stage2/01_bases/execution/build_gene_axis.py --verify',
    'python stage2/01_bases/figure/script/render_gene_axis_rank_sum_histogram.py '
    '--csv stage2/01_bases/figure/data/gene_axis_rank_sum_histogram.csv '
    '--png stage2/01_bases/figure/gene_axis_rank_sum_histogram.png',
    'python stage2/01_bases/execution/make_gene_axis_docx.py '
    '--coverage stage2/01_bases/output/gene-axis/source_topn_coverage.csv '
    '--figure stage2/01_bases/figure/gene_axis_rank_sum_histogram.png '
    '--out stage2/01_bases/output/gene-axis/gene_axis_fill_rule.docx',
)

FIXED_ZIP_TIME = (1980, 1, 1, 0, 0, 0)


def normalize_docx(path):
    """Re-zip a DOCX with fixed entry timestamps so its bytes replay deterministically."""
    path = Path(path)
    buffer = io.BytesIO()
    with zipfile.ZipFile(path) as source, zipfile.ZipFile(buffer, 'w') as target:
        for info in source.infolist():
            entry = zipfile.ZipInfo(info.filename, date_time=FIXED_ZIP_TIME)
            entry.compress_type = info.compress_type
            entry.external_attr = info.external_attr
            entry.create_system = info.create_system
            target.writestr(entry, source.read(info.filename))
    path.write_bytes(buffer.getvalue())
    return path


def load_metadata_module(path):
    """Import a ``gene_axis_source_metadata``-shaped module from an arbitrary file path.

    This is the ``--metadata`` CLI input: loading by path (rather than by the
    package name) lets the generator be pointed at a modified copy, which is how
    the wiring test proves the argument is actually read.
    """
    path = Path(path)
    spec = importlib.util.spec_from_file_location('gene_axis_source_metadata_input', path)
    if spec is None or spec.loader is None:
        raise ValueError(f'{path}: cannot be imported as a Python module')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_coverage(path):
    """Read and validate the released coverage table (long format, one row per group and axis).

    The v2 contract carries the numerator and the denominator explicitly, so a
    reader never has to guess whether a row is divided by ``N`` (tenv-HVG) or by
    the disease's own ``|DEG2000_d|``.
    """
    table = pd.read_csv(path)
    if list(table.columns) != list(COVERAGE_COLUMNS):
        raise ValueError(f'{path}: coverage columns differ from the v2 contract '
                         f'{list(COVERAGE_COLUMNS)}')
    if len(table) == 0:
        raise ValueError(f'{path}: coverage table is empty')
    keys = list(zip(table['N'], table['group']))
    if len(set(keys)) != len(keys):
        raise ValueError(f'{path}: coverage table repeats an (N, group) row')
    for row in table.itertuples(index=False):
        if pd.isna(row.numerator) and pd.isna(row.denominator) and pd.isna(row.rate):
            if not str(row.note).strip():
                raise ValueError('excluded coverage requires a reason')
            continue
        if not 0 <= row.numerator <= row.denominator:
            raise ValueError(f'{path}: N={row.N} group={row.group!r} has numerator '
                             f'{row.numerator} outside 0..{row.denominator}')
    return table


def _accession_cell(source):
    """Render the accession column: the accession plus its series note if any."""
    accession = source['accession']
    note = source['accession_note']
    if accession is None:
        return note
    return f'{accession} ({note})' if note else accession


def _shared_control_vehicle(shared_control_group):
    """Derive the vehicle phrase from a cardiotox catalog's ``shared_control_group``.

    ``"GSE243674:same-time vehicle wells (6 donors ...)"`` becomes
    ``"same-time vehicle wells"``: the catalog publishes the wording, so keeping
    no second hard-coded copy is what lets table wording drift fail a test.
    """
    _, separator, tail = shared_control_group.partition(':')
    if not separator:
        raise ValueError(f'shared_control_group has no ":" prefix: {shared_control_group!r}')
    vehicle = tail.strip().split(' (', 1)[0]
    if not vehicle:
        raise ValueError(f'shared_control_group has no vehicle phrase: {shared_control_group!r}')
    return vehicle


def _groups(dataset):
    """Render the groups column: the catalog's arm, or the cardiotox class rule."""
    arm = dataset.get('arm_or_class')
    if arm is not None:
        return arm
    if dataset.get('category') == 'TOX':
        return f"TOX pooled {len(dataset['member_drugs'])} admitted drugs vs matched vehicle"
    vehicle = _shared_control_vehicle(dataset['shared_control_group'])
    drugs = ', '.join(dataset['member_drugs'])
    return f"{dataset['class_label']} ({dataset['class_id']}; {drugs}) vs {vehicle}"


def _tenv_hvg_row(metadata):
    source = metadata.SOURCES[metadata.TENV_HVG]
    row = dict(TENV_HVG_ROW)
    row['accession'] = _accession_cell(source)
    row['sequencing technology'] = source['technology']
    return row


def read_source_table_rows(catalogs, metadata):
    """Read the declared catalogs into the source table's 34 data rows.

    ``catalogs`` maps category -> ``catalog.json`` path and ``metadata`` is the
    verified mapping module. Rows follow the released
    ``gene_axis.RELEASE_CATEGORIES`` order and each catalog's ``datasets[]``
    file order; the tenv-HVG fusion row is appended last. Row counts, sample
    type, groups, sample sizes and the v2 CM DEG counts (``n_deg``, never the
    padj-only ``n_significant``) come from the catalogs; only
    accession and sequencing technology come from the mapping. The metadata
    coverage is validated before any lookup, so a missing or extra dataset fails
    the build with an explicit error.
    """
    datasets = []
    dataset_ids = []
    for category in gene_axis.RELEASE_CATEGORIES:
        path = Path(catalogs[category])
        catalog = json.loads(path.read_text(encoding='utf-8'))
        if catalog.get('category') != gene_axis.CATEGORY_LABELS[category]:
            raise ValueError(f'{path}: catalog category {catalog.get("category")!r} '
                             f'differs from {gene_axis.CATEGORY_LABELS[category]!r}')
        for dataset in catalog['datasets']:
            dataset_ids.append(dataset['dataset_id'])
            datasets.append((catalog, dataset))
    metadata.validate(dataset_ids + [metadata.TENV_HVG])
    rows = []
    for catalog, dataset in datasets:
        source = metadata.SOURCES[dataset['dataset_id']]
        rows.append({
            'category': catalog['category'],
            'accession': _accession_cell(source),
            'sample type': dataset.get('model_system', 'human stem cell cardiomyocyte bulk'),
            'groups': _groups(dataset),
            'sample sizes': f"{dataset['n_case']} / {dataset['n_control']}",
            'sequencing technology': source['technology'],
            'CM DEG count': f"{dataset['n_deg']:,}",
        })
    rows.append(_tenv_hvg_row(metadata))
    return rows


def read_axis_symbols(path):
    """Read the frozen axis identity table the catalogs declare (``axis.path``)."""
    table = pd.read_csv(path, usecols=['gene_index', 'symbol'])
    if len(table) != len(set(table['gene_index'])):
        raise ValueError(f'{path}: duplicate gene_index in the frozen axis identity table')
    return dict(zip(table['gene_index'].astype(int), table['symbol'].astype(str)))


def read_manifest(path):
    """Read the release manifest and refuse anything that is not the v2 release.

    ``release.algorithm`` is the self-describing marker #207 added: a manifest
    without it predates the ADR 0071 rule, so rendering its numbers under the v2
    rule would describe them falsely. The returned mapping is the published
    tenv-HVG fallback slot count per axis size.
    """
    manifest = json.loads(Path(path).read_text(encoding='utf-8'))
    algorithm = (manifest.get('release') or {}).get('algorithm')
    if algorithm != gene_axis.ALGORITHM:
        raise ValueError(f'{path}: release.algorithm is {algorithm!r}, not '
                         f'{gene_axis.ALGORITHM!r}; this document only renders the '
                         f'v2 release')
    slots = (manifest['release'].get('fallback_second_half_slots') or {})
    return {int(size): int(count) for size, count in slots.items()}


def fallback_sentence(slots):
    """State the published tenv-HVG fallback count instead of asserting zero."""
    if not slots:
        return '本次发布的补齐数量未记录在 manifest 中。'
    if not any(slots.values()):
        return '本次发布两个轴的缺口都是 0（manifest 的 fallback_second_half_slots）。'
    return ('本次发布的补齐数量：'
            + '、'.join(f'N={size} 用 {count} 个 tenv-HVG 基因补齐'
                        for size, count in sorted(slots.items()))
            + '（manifest 的 fallback_second_half_slots）。')


def read_disease_table_rows(catalogs, fusion, fusion_path, symbols):
    """Read the published per-disease record into the second table's declared rows.

    ``fusion`` is ``output/disease-fusion/catalog.json`` -- the published level-1
    record the gene axis itself consumes -- and ``symbols`` maps ``gene_index`` to
    the frozen-axis symbol so a very small ``DEG2000`` can be named instead of
    only counted. Each row is cross-checked against that disease's
    ``output/disease-deg/<category>/catalog.json``: the source count, the DEG
    union and ``DEG2000`` must agree, so a build from mismatched products fails
    instead of publishing two contradictory tables.
    """
    by_label = {record['category']: record for record in fusion['diseases']}
    if set(by_label) != {gene_axis.CATEGORY_LABELS[category]
                         for category in gene_axis.RELEASE_CATEGORIES}:
        raise ValueError('the fusion catalog diseases differ from '
                         f'{list(gene_axis.CATEGORY_LABELS.values())}')
    n_reached = sum(1 for record in fusion['diseases']
                    if record['deg2000'] == gene_axis.DEG2000_LIMIT
                    and record['n_deg'] > record['deg2000'])
    rows = []
    for category in gene_axis.RELEASE_CATEGORIES:
        label = gene_axis.CATEGORY_LABELS[category]
        catalog = json.loads(Path(catalogs[category]).read_text(encoding='utf-8'))
        record = by_label[label]
        n_sources = len(catalog['datasets'])
        tiers = sorted({entry['tier'] for entry in catalog['tier_sources']})
        for published, recomputed in (
            (record['declared_source_count'], n_sources),
            (record['n_deg'], catalog['n_deg']),
            (record['deg2000'], catalog['deg2000']),
        ):
            if published != recomputed:
                raise ValueError(f'{label}: disease-fusion catalog says {published}, '
                                 f'disease-deg catalog says {recomputed}')
        deg2000_genes = _deg2000_symbols(fusion_path, record, symbols)
        rows.append({
            'disease': label,
            'source datasets': str(n_sources),
            'tier path': ' → '.join(str(tier) for tier in tiers),
            'DEG union': f"{record['n_deg']:,}",
            'DEG2000': f"{record['deg2000']:,}",
            'eligible': 'yes' if record.get('eligible', True) else 'no',
            'note': _disease_note(record, tiers, n_reached, deg2000_genes) + ('' if record.get('eligible', True) else '；排除跨病融合：轴内 DEG 并集少于 100'),
        })
    return rows


def _deg2000_symbols(fusion_path, record, symbols):
    """Resolve the disease's published ``DEG2000`` gene indices to symbols."""
    path = Path(fusion_path).parent / record['deg2000_file']
    table = pd.read_csv(path)
    if len(table) != record['deg2000']:
        raise ValueError(f'{path}: {len(table)} rows but the catalog says '
                         f'{record["deg2000"]} DEG2000 genes')
    return [symbols[index] for index in table['gene_index'].astype(int)]


def _disease_note(record, tiers, n_reached, deg2000_genes):
    """One honest sentence per disease: the cut, or the shortfall and its reason."""
    limit = gene_axis.DEG2000_LIMIT
    n_deg, deg2000 = record['n_deg'], record['deg2000']
    if deg2000 == limit:
        if n_deg == deg2000:
            return f'恰好达到 {limit:,}：DEG 并集 = DEG2000 = {limit:,} 个，无需截断'
        return (f'达到 {limit:,} 截断：DEG 并集 {n_deg:,}，DEG2000 取前 {limit:,} 个'
                + ('（本次唯一达到截断的类别）' if n_reached == 1 else ''))
    if n_deg != deg2000:
        return (f'DEG 并集 {n_deg:,}，DEG2000 {deg2000:,}；'
                f'未达到 {limit:,}，差 {limit - deg2000:,}')
    note = (f'未达到 {limit:,}：停在 tier {max(tiers)}，DEG 并集 {n_deg:,}，'
            f'差 {limit - deg2000:,}；声明的层级已穷尽、没有更低层可收，'
            f'DEG2000 即整个并集（不填充）')
    if 0 < deg2000 <= NAMED_DEG2000_LIMIT:
        note += ('。DEG2000 就是 '
                 + '、'.join(deg2000_genes)
                 + (' 这 1 个基因' if deg2000 == 1 else f' 这 {deg2000} 个基因'))
    excluded = record.get('excluded_sources') or []
    if excluded:
        note += (f'；{len(excluded)} 个来源 n_DEG = 0（'
                 + '、'.join(source['dataset_id'] for source in excluded)
                 + '），按 ADR 0068 不计入融合')
    return note


def _heading(document, text, level=1):
    return document.add_heading(text, level=level)


def _body(document, text, bold_lead=None):
    paragraph = document.add_paragraph()
    if bold_lead:
        run = paragraph.add_run(bold_lead)
        run.bold = True
    paragraph.add_run(text)
    return paragraph


def _grid(document, columns, rows, font_size):
    """Add one ``Table Grid`` table with a bold header and plain data rows."""
    grid = document.add_table(rows=1 + len(rows), cols=len(columns))
    grid.style = 'Table Grid'
    grid.autofit = False
    widths = ([3.5, 4.2, 4.8, 5.0, 2.0, 5.2, 2.0] if tuple(columns) == SOURCE_TABLE_COLUMNS
              else [3.5, 1.8, 1.8, 1.8, 1.8, 1.8, 14.2] if tuple(columns) == DISEASE_TABLE_COLUMNS
              else None)
    if widths:
        for column, width in zip(grid.columns, widths):
            column.width = Cm(width)
        for row in grid.rows:
            for cell, width in zip(row.cells, widths):
                cell.width = Cm(width)
    for column, name in enumerate(columns):
        cell = grid.cell(0, column)
        cell.text = ''
        run = cell.paragraphs[0].add_run(name)
        run.bold = True
        run.font.size = Pt(font_size)
    for row_index, values in enumerate(rows, start=1):
        for column, value in enumerate(values):
            cell = grid.cell(row_index, column)
            cell.text = ''
            run = cell.paragraphs[0].add_run(value)
            run.font.size = Pt(font_size)
    header = OxmlElement('w:tblHeader')
    header.set(qn('w:val'), '1')
    grid.rows[0]._tr.get_or_add_trPr().append(header)
    for row in grid.rows:
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                paragraph.paragraph_format.space_after = Pt(2)
                paragraph.paragraph_format.space_before = Pt(2)
                paragraph.paragraph_format.line_spacing = 1.0
        avoid_split = OxmlElement('w:cantSplit')
        row._tr.get_or_add_trPr().append(avoid_split)
        avoid_split.set(qn('w:val'), '1')
    return grid


def _percent(value):
    return f'{value:.2f}%'


def coverage_note(disease_rows):
    """The coverage table's reading rule, with the denominators read off table 2."""
    denominators = '、'.join(f"{row['disease']} {row['DEG2000']}" for row in disease_rows if row['eligible'] == 'yes')
    return (
        '每行 = 一个组（tenv-HVG 或一个疾病）在一个轴 N 上的覆盖率。tenv-HVG 行沿用 v1 口径：'
        '分子 = 其 top-N 集合中落进 G_N 的基因数，分母 = N。疾病行是 v2 口径：'
        f'分子 = |DEG2000_d ∩ G_N|，分母 = |DEG2000_d|（{denominators}），'
        'AF 的一个 DEG2000 基因保留，但因轴内 DEG 并集少于 100 被排除；其覆盖率三项留空。各行百分比不能跨组比较。'
    )


def _check_coverage_denominators(table, disease_rows):
    """The coverage table's disease denominators must be each disease's ``DEG2000``.

    tenv-HVG keeps ``N`` (its top-``N`` set), every disease divides by its own
    ``|DEG2000_d|``. Reading a disease row against ``N`` is the v1 misreading the
    long-format contract exists to prevent, so the document refuses to render a
    table whose two halves disagree.
    """
    published = {row['disease']: int(row['DEG2000'].replace(',', ''))
                 for row in disease_rows}
    for record in table.itertuples(index=False):
        if record.group == source_metadata.TENV_HVG:
            expected = record.N
        else:
            label = gene_axis.CATEGORY_LABELS.get(record.group)
            if label is None:
                raise ValueError(f'coverage table has an unknown group {record.group!r}')
            expected = published[label]
        excluded = record.group != source_metadata.TENV_HVG and next(row for row in disease_rows if row['disease'] == label)['eligible'] == 'no'
        if excluded:
            if not all(pd.isna(value) for value in (record.numerator, record.denominator, record.rate)):
                raise ValueError('excluded disease coverage must be blank')
            continue
        if record.denominator != expected:
            raise ValueError(f'coverage N={record.N} group={record.group!r} divides by '
                             f'{record.denominator}, expected {expected}')


def build(coverage_path, figure_path, output_path, catalogs=None, metadata=None,
          fusion_catalog=None, axis_path=None, manifest_path=None):
    """Render the document; every path defaults to the released repo layout."""
    defaults = default_paths()
    metadata = source_metadata if metadata is None else metadata
    if catalogs is None:
        catalogs = defaults['catalogs']
    fusion_catalog = defaults['fusion'] if fusion_catalog is None else Path(fusion_catalog)
    axis_path = defaults['axis'] if axis_path is None else Path(axis_path)
    manifest_path = defaults['manifest'] if manifest_path is None else Path(manifest_path)
    fallback_slots = read_manifest(manifest_path)
    table = read_coverage(coverage_path)
    source_rows = read_source_table_rows(catalogs, metadata)
    fusion = json.loads(Path(fusion_catalog).read_text(encoding='utf-8'))
    symbols = read_axis_symbols(axis_path)
    disease_rows = read_disease_table_rows(catalogs, fusion, fusion_catalog, symbols)
    _check_coverage_denominators(table, disease_rows)

    document = Document()
    section = document.sections[0]
    section.page_width = Cm(29.7)
    section.page_height = Cm(21)
    section.top_margin = Cm(1.5)
    section.bottom_margin = Cm(1.5)
    section.left_margin = Cm(1.5)
    section.right_margin = Cm(1.5)
    style = document.styles['Normal']
    style.font.size = Pt(10)
    style.paragraph_format.space_after = Pt(6)
    style.paragraph_format.line_spacing = 1.05
    style.font.name = 'Noto Sans CJK SC'
    style.element.rPr.rFonts.set(qn('w:eastAsia'), 'Noto Sans CJK SC')
    for name in ('Title', 'Heading 1', 'Heading 2'):
        document.styles[name].font.color.rgb = RGBColor(0, 0, 0)
        document.styles[name].font.name = 'Noto Sans CJK SC'

    compact = document.styles['No Spacing']
    compact.font.name = 'Noto Sans CJK SC'
    compact.font.size = Pt(8)
    compact.paragraph_format.line_spacing = 1.0

    title = document.add_heading(TITLE, level=0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for element in (title._p, document.styles['Title'].element):
        for borders in list(element.iter(qn('w:pBdr'))):
            borders.getparent().remove(borders)
    _body(document,
          '依据 ADR 0071 的 DEG 规则及 ADR 0072 的统一疾病资格与 TOX 定义。五个声明疾病中，'
          'nHF、ICM、genetic_cardiomyopathy 和 TOX 合格；AF 保留证据但排除跨病融合。'
          '本文件说明基因轴选择规则及各疾病实际证据，供读取和使用 G_2000 与 G_5000 时核对。'
          '规则文字是常量，并与 bases.gene_axis.FILL_RULE 保持一致；'
          '表格与图来自已发布产品，可由文末命令重放。')

    _heading(document, RULE_HEADING)
    _body(document, DEG_DEFINITION)
    steps = list(RULE_STEPS)
    steps[5] = steps[5] + fallback_sentence(fallback_slots)
    for index, step in enumerate(steps, start=1):
        document.add_paragraph(f'{index}) {step}')
    _body(document, F2000_DEFINITION, bold_lead='f2000 定义：')

    document.add_page_break()
    _heading(document, SOURCE_TABLE_HEADING)
    # Bound each source-table page explicitly: LO's automatic repeated-header
    # layout can paint the first line of a long row under the header.
    for offset in range(0, len(source_rows), 5):
        if offset:
            document.add_page_break()
        caption = document.add_paragraph(SOURCE_TABLE_CAPTION + (f' 续表 {offset // 5 + 1}' if offset else ''))
        caption.alignment = WD_ALIGN_PARAGRAPH.CENTER
        _grid(document, SOURCE_TABLE_COLUMNS,
              [[record[name] for name in SOURCE_TABLE_COLUMNS] for record in source_rows[offset:offset + 5]],
              font_size=8.5)
    _body(document, SOURCE_TABLE_FOOTNOTE)

    document.add_page_break()
    _heading(document, DISEASE_TABLE_HEADING)
    caption = document.add_paragraph(DISEASE_TABLE_CAPTION)
    caption.alignment = WD_ALIGN_PARAGRAPH.CENTER
    _grid(document, DISEASE_TABLE_COLUMNS,
          [[record[name] for name in DISEASE_TABLE_COLUMNS] for record in disease_rows],
          font_size=8.5)
    _body(document, DISEASE_TABLE_FOOTNOTE)

    _heading(document, '4 四处易错点')
    for lead, text in CLARIFICATIONS:
        _body(document, text, bold_lead=lead + ' ')

    _heading(document, '5 每个疾病 DEG2000 的轴覆盖率')
    _body(document, coverage_note(disease_rows))
    coverage_rows = [
        [str(int(record.N)), record.group_label, ('excluded' if pd.isna(record.numerator) else str(int(record.numerator))),
         ('—' if pd.isna(record.denominator) else str(int(record.denominator))), ('—' if pd.isna(record.rate) else _percent(record.rate))]
        for record in table.itertuples(index=False)
    ]
    _grid(document, COVERAGE_TABLE_COLUMNS, coverage_rows, font_size=8.5)

    document.add_page_break()
    _heading(document, '6 证据图')
    document.add_picture(str(figure_path), width=Cm(15.5))
    caption = document.add_paragraph(
        '图 1　每个子图对应一个基因轴，标题为 N = 2,000 / N = 5,000。灰柱是全部 '
        '18,164 个冻结轴候选基因的 rank_sum 分布；每根柱自顶端向下的彩色段表示该柱中'
        '被选进轴的基因数量（彩色总量 = N）。')
    caption.alignment = WD_ALIGN_PARAGRAPH.CENTER

    _heading(document, '7 复现命令与相关资产')
    for command in REPRODUCE:
        document.add_paragraph(command, style='No Spacing')
    assets = (
        ('candidate_rank_sum.csv',
         '主候选秩和表：18,164 行，含 tenv_hvg_rank、声明疾病秩及资格、rank_sum、f2000'),
        ('gene_axis_2000.csv / gene_axis_5000.csv', '两个已发布轴（gene_index, symbol）'),
        ('axis_composition.csv', '每轴的 axis_fill / DEG2000 / f2000 构成记录'),
        ('source_topn_coverage.csv', '本文件的覆盖率表（长表：每组 N、分子、分母、比率）'),
        ('manifest.json', '输入/输出摘要与 sha256（release.algorithm = f2000-deg2000-v2）'),
        ('disease-fusion/catalog.json', '每病 level-1 融合与 DEG2000 的发布记录（表 2 的来源）'),
        ('figure/data/gene_axis_rank_sum_histogram.csv', '图中全部绘制值'),
    )
    _grid(document, ('资产', '内容'), assets, font_size=8)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    document.save(output_path)
    normalize_docx(output_path)
    return output_path


def default_paths():
    root = next(parent for parent in Path(__file__).resolve().parents
                if (parent / '.git').exists())
    step = root / 'stage2/01_bases'
    return {
        'coverage': step / 'output/gene-axis/source_topn_coverage.csv',
        'figure': step / 'figure/gene_axis_rank_sum_histogram.png',
        'fusion': step / 'output/disease-fusion/catalog.json',
        'axis': step / 'output/gene-identity/genes.csv',
        'manifest': step / 'output/gene-axis/manifest.json',
        'catalogs': {category: step / 'output/disease-deg' / category / 'catalog.json'
                     for category in gene_axis.RELEASE_CATEGORIES},
        'out': step / 'output/gene-axis/gene_axis_fill_rule.docx',
    }


def main():
    defaults = default_paths()
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--coverage', type=Path, default=defaults['coverage'])
    parser.add_argument('--figure', type=Path, default=defaults['figure'])
    parser.add_argument('--catalog', action='append', default=None,
                        metavar='CATEGORY=PATH',
                        help='override one released catalog path (repeatable; default: '
                             f'the declared repo catalogs {sorted(defaults["catalogs"])})')
    parser.add_argument('--fusion-catalog', type=Path, default=defaults['fusion'],
                        help='override the published disease-fusion catalog (table 2)')
    parser.add_argument('--axis', type=Path, default=defaults['axis'],
                        help='override the frozen axis identity table (gene_index -> symbol)')
    parser.add_argument('--manifest', type=Path, default=defaults['manifest'],
                        help='override the release manifest (guards on release.algorithm)')
    parser.add_argument('--metadata', type=Path, default=None,
                        metavar='PATH',
                        help='override the tracked metadata-module path (imported from '
                             f'that file; default: {Path(source_metadata.__file__).resolve()})')
    parser.add_argument('--out', type=Path, default=defaults['out'])
    args = parser.parse_args()
    catalogs = dict(defaults['catalogs'])
    for override in args.catalog or ():
        category, separator, value = override.partition('=')
        if not separator or category not in catalogs:
            parser.error(f'--catalog expects CATEGORY=PATH with CATEGORY in '
                         f'{sorted(catalogs)}, got {override!r}')
        catalogs[category] = Path(value)
    metadata = (source_metadata if args.metadata is None
                else load_metadata_module(args.metadata))
    path = build(args.coverage, args.figure, args.out, catalogs, metadata=metadata,
                 fusion_catalog=args.fusion_catalog, axis_path=args.axis,
                 manifest_path=args.manifest)
    print(f'wrote {path} from {args.coverage}, {args.figure}, {len(catalogs)} catalogs, '
          f'{args.fusion_catalog}, {args.axis}, {args.manifest} and {metadata.__file__}')


if __name__ == '__main__':
    main()
