#!/usr/bin/env python3
"""Build the portrait A4 ADR0073 gene-axis rule document from one release.

Reads published disease/catalog, complete tier-first fusion, four-axis coverage
and evidence figure inputs. Source tables are separate per declared disease;
tenv-HVG is an explicitly separate non-disease source. No DEG analysis or axis
selection is recomputed here. Defaults write stage2/01_bases/gene_axis_fill_rule.docx.
All input paths may be overridden for an isolated coherent candidate. DOCX ZIP
entry times are fixed so the same inputs/toolchain reproduce identical bytes.
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
                        'sample sizes', 'sequencing technology', 'axis DEG count')
SOURCE_TABLE_CAPTION = '表 1　已发布来源数据集信息'
SOURCE_TABLE_HEADING = '2 来源数据集信息'
DISEASE_TABLE_COLUMNS = ('disease', 'source datasets', 'tier path',
                         'DEG union', 'DEG1000', 'eligible', 'note')
DISEASE_TABLE_CAPTION = '表 2　每个疾病的来源数据集数、DEG 并集与 DEG1000'
DISEASE_TABLE_HEADING = '3 每个疾病的来源数与 DEG 资格'
#: A ``DEG1000`` this small is named gene by gene rather than described as a
#: count, so a single-gene result cannot look like a rounding of a larger one.
NAMED_DEG1000_LIMIT = 5

RULE_HEADING = '1 取用规则 ADR 0073'

DEG_DEFINITION = (
    'DEG 定义：一个基因在某数据集里是 DEG，当它在至少一个已声明对比中满足 '
    'padj < 0.05 且 |log2FC| ≥ 2（FC ≥ 4）。任一对比通过即算（pass-any-contrast），'
    '不是“最佳对比”塌缩出的单一胜负集；padj 为 NA 视为不显著；效应方向被记录，'
    '但不参与筛选。本文件的 轴内 DEG 数就是这一口径下的基因数，'
    '它与 padj<0.05 的 padj-only 显著基因数不是同一个量。'
)

RULE_STEPS = (
    '第一步（来源内秩）：每个来源只在自己的 DEG 内排序，秩取 1..n_DEG，按效应 '
    '|log2FC|·(−log10 padj) 降序、gene_index 破平；非 DEG 一律记 n_DEG+1；'
    'rank_norm = rank / n_DEG。分母是 DEG 数，不是测得基因数。',
    '第二步（完整来源收集）：遍历每病全部合格的声明层级，不因累计 DEG 数达到 '
    '1,000 或 2,000 而停止。tier 1 = 人房室组织 CM 特异分析；tier 2 = 人干细胞来源 '
    'CM、类器官或工程心脏组织，临床类别仅收非遗传模型；tier 3 = 人房室全组织 bulk。'
    '未声明的层不补假来源；保留来源质量、独立数据集、药物证据与去重要求。',
    '第三步（每病融合）：disease_rank_sum = Σ来源 rank_norm，数值不加 tier 权重。'
    'source_tier_highest 是该基因实际 DEG 阳性的合格信息来源 tier 的最小值；'
    '仅测到该基因不算支持。完整 DEG 并集按 source_tier_highest 升序、'
    'disease_rank_sum 升序、gene_index 升序取得连续 1..n_d 秩。DEG1000_d '
    '精确取前 min(1,000, n_d) 个，不填充；并集及未入选基因的来源证据仍完整保留。',
    '第四步（疾病资格）：完整冻结轴 DEG 并集 n_d ≥ 100 才进入跨病融合。'
    '未达门槛的疾病保留上游证据与排除理由；AF 同样由数据判定，不作永久排除。'
    'TOX 每个独立数据集将已准入药物与匹配溶媒比较后提供一个来源秩，机制不增加票数。',
    '第五步（跨病分数与选轴）：N ∈ {2,000, 3,000, 4,000, 5,000}。'
    'f1000 是 DEG1000_d 含该基因的合格疾病数；tenv-HVG 与排除疾病均不贡献。'
    'rank_sum = Σ合格疾病 rank_d / n_d，分母 n_d 为完整并集大小。'
    'DEG1000 外的基因统一用 rank_d = n_d+1，包括真实 DEG 的第 1,001 位。'
    '前 N/2 取 tenv-HVG，后半从未选 f1000 ≥ 1 基因按 f1000 降序、rank_sum 升序、'
    'gene_index 升序选取；DEG1000 不限制 HVG 半轴或完整候选表。',
    '第六步（补齐与发布）：疾病候选不足 N/2 时，按 HVG 秩用未选 tenv-HVG 基因补齐。'
    '每档记录实际 fallback；合并后按冻结 gene_index 升序发布，不重新编号。'
    '四轴精确嵌套并保留原始 gene_index 与 symbol；仅发布 G_5000 的 tenv 和 tahoe_bulk views。 '
    '每行先取 5,000 个原始计数再计算库总量，CP10K 后取自然 log1p；每个溶媒独立归一化后求匹配对照均值，delta = X − control（ADR 0077）。',
)

F1000_DEFINITION = (
    '合格疾病中 DEG1000_d 含该基因的疾病数量。资格和数量以本次 catalog 为准；'
    'tenv-HVG 不计入，f1000 = 0 的基因只可由 HVG 分配或补齐进入轴。'
)
CLARIFICATIONS = (
    ('pN 已退役', 'f1000 只计合格疾病的 DEG1000，不计 tenv-HVG。'),
    ('排序方向', '疾病内先 tier 升序，再数值 rank sum；跨病先 f1000 降序，再 rank_sum 升序。'),
    ('小集合与资格不同', 'DEG1000 不填充；只有完整轴内并集少于 100 才排除跨病贡献。'),
    ('覆盖率分母', '合格疾病用自己的 DEG1000 大小；tenv-HVG 用 N。排除行为空并保留理由。'),
)
SOURCE_TABLE_FOOTNOTE = (
    '轴内 DEG 数来自各 catalog 的 n_deg，是 pass-any-contrast DEG 数，不是 padj-only 显著基因数。'
    '样本量口径：非 TOX 来源为已发布病例和对照数；TOX 是已声明对比的 pooled wells，'
    '有多个对比时显示范围，n_case_basis 与 n_control_basis 记录计数口径，不能作为全系列独立供体总数。TOX 表中药物仅指独立心脏毒性证据准入药物；'
    '保留匹配溶媒、供体阻断及原有细胞系、板和时间分层。没有要求每个剂量和时间都出现毒性。'
)
DISEASE_TABLE_FOOTNOTE = (
    'DEG union 是全部来源解析到冻结 18164 轴的去重 DEG 并集。DEG1000 是 tier 优先排名前 1000 的基因，'
    '不足时不填充。eligible 来自完整并集的 100 基因门槛；排除类别保留证据但不贡献跨病分数。'
    '各数值来自 disease-fusion/catalog.json 并与 disease-deg catalogs 交叉校验。'
)

#: The tenv-HVG fusion row is not a case/control contrast. Its sample sizes are
#: the published corpus populations pinned in ``output/view-5k/manifest.json``
#: ``sources.tenv`` (treated = 4,206, vehicle = 272; also
#: ``log/raw_counts_provenance.json``: 4,510 native wells kept). The two numbers
#: are declared once, in ``gene_axis_source_metadata``, and imported here.
TENV_HVG_ROW = {
    'category': 'tenv-HVG (fusion source)',
    'sample type': 'hiPSC-CM 96-well plate lysates (TempO-Seq corpus)',
    'groups': 'not a contrast: ranked by the published per-well dispersion rank',
    'sample sizes': f'{TENV_TREATED_WELLS:,} / {TENV_VEHICLE_WELLS:,} (treated / vehicle wells)',
    'axis DEG count': 'not a DEG source',
}

REPRODUCE = (
    'python stage2/01_bases/execution/build_gene_axis.py verify',
    'python stage2/01_bases/figure/script/render_gene_axis_rank_sum_histogram.py '
    '--csv stage2/01_bases/figure/data/gene_axis_rank_sum_histogram.csv '
    '--png stage2/01_bases/figure/gene_axis_rank_sum_histogram.png',
    'python stage2/01_bases/execution/make_gene_axis_docx.py '
    '--coverage stage2/01_bases/output/gene-axis/source_topn_coverage.csv '
    '--figure stage2/01_bases/figure/gene_axis_rank_sum_histogram.png '
    '--out stage2/01_bases/gene_axis_fill_rule.docx',
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
    the disease's own ``|DEG1000_d|``.
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
    """Read the declared catalogs into the complete source table.

    ``catalogs`` maps category -> ``catalog.json`` path and ``metadata`` is the
    verified mapping module. Rows follow the released
    ``gene_axis.RELEASE_CATEGORIES`` order and each catalog's ``datasets[]``
    file order; the tenv-HVG fusion row is appended last. Row counts, sample
    type, groups, sample sizes and the axis-resolved DEG counts (``n_deg``, never the
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
            'axis DEG count': f"{dataset['n_deg']:,}",
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
    """Read the release manifest and refuse anything that is not the ADR 0073 release.

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
                         f'ADR 0073 release')
    slots = (manifest['release'].get('fallback_second_half_slots') or {})
    parsed = {int(size): int(count) for size, count in slots.items()}
    if set(parsed) != set(gene_axis.AXIS_SIZES) or any(count < 0 or count > size // 2 for size, count in parsed.items()):
        raise ValueError('release manifest must declare bounded fallback for all four axes')
    return parsed


def fallback_sentence(slots):
    """State the published tenv-HVG fallback count instead of asserting zero."""
    if not slots:
        raise ValueError('release manifest must declare actual fallback counts')
    if not any(slots.values()):
        return '本次发布四个轴的缺口都是 0（manifest 的 fallback_second_half_slots）。'
    return ('本次发布的补齐数量：'
            + '、'.join(f'N={size} 用 {count} 个 tenv-HVG 基因补齐'
                        for size, count in sorted(slots.items()))
            + '（manifest 的 fallback_second_half_slots）。')


def read_disease_table_rows(catalogs, fusion, fusion_path, symbols):
    """Read the published per-disease record into the second table's declared rows.

    ``fusion`` is ``output/disease-fusion/catalog.json`` -- the published level-1
    record the gene axis itself consumes -- and ``symbols`` maps ``gene_index`` to
    the frozen-axis symbol so a very small ``DEG1000`` can be named instead of
    only counted. Each row is cross-checked against that disease's
    ``output/disease-deg/<category>/catalog.json``: the source count, the DEG
    union and ``DEG1000`` must agree, so a build from mismatched products fails
    instead of publishing two contradictory tables.
    """
    by_label = {record['category']: record for record in fusion['diseases']}
    if set(by_label) != {gene_axis.CATEGORY_LABELS[category]
                         for category in gene_axis.RELEASE_CATEGORIES}:
        raise ValueError('the fusion catalog diseases differ from '
                         f'{list(gene_axis.CATEGORY_LABELS.values())}')
    n_reached = sum(1 for record in fusion['diseases']
                    if record['deg1000'] == gene_axis.DEG1000_LIMIT
                    and record['n_deg'] > record['deg1000'])
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
            (record['deg1000'], catalog['deg1000']),
        ):
            if published != recomputed:
                raise ValueError(f'{label}: disease-fusion catalog says {published}, '
                                 f'disease-deg catalog says {recomputed}')
        deg1000_genes = _deg1000_symbols(fusion_path, record, symbols)
        rows.append({
            'disease': label,
            'source datasets': str(n_sources),
            'tier path': ' → '.join(str(tier) for tier in tiers),
            'DEG union': f"{record['n_deg']:,}",
            'DEG1000': f"{record['deg1000']:,}",
            'eligible': 'yes' if record.get('eligible', True) else 'no',
            'note': _disease_note(record, tiers, n_reached, deg1000_genes) + ('' if record.get('eligible', True) else '；排除跨病融合：轴内 DEG 并集少于 100'),
        })
    return rows


def _deg1000_symbols(fusion_path, record, symbols):
    """Resolve the disease's published ``DEG1000`` gene indices to symbols."""
    path = Path(fusion_path).parent / record['deg1000_file']
    table = pd.read_csv(path)
    if len(table) != record['deg1000']:
        raise ValueError(f'{path}: {len(table)} rows but the catalog says '
                         f'{record["deg1000"]} DEG1000 genes')
    return [symbols[index] for index in table['gene_index'].astype(int)]


def _disease_note(record, tiers, n_reached, deg1000_genes):
    """One honest sentence per disease: the cut, or the shortfall and its reason."""
    limit = gene_axis.DEG1000_LIMIT
    n_deg, deg1000 = record['n_deg'], record['deg1000']
    if deg1000 == limit:
        if n_deg == deg1000:
            return f'恰好达到 {limit:,}：DEG 并集 = DEG1000 = {limit:,} 个，无需截断'
        return (f'达到 {limit:,} 截断：DEG 并集 {n_deg:,}，DEG1000 取前 {limit:,} 个'
                + ('（本次唯一达到截断的类别）' if n_reached == 1 else ''))
    if n_deg != deg1000:
        return (f'DEG 并集 {n_deg:,}，DEG1000 {deg1000:,}；'
                f'未达到 {limit:,}，差 {limit - deg1000:,}')
    note = (f'未达到 {limit:,}：DEG 并集 {n_deg:,}，差 {limit - deg1000:,}；'
            '已遍历全部合格声明层级，DEG1000 即整个并集（不填充）')
    if 0 < deg1000 <= NAMED_DEG1000_LIMIT:
        note += ('。DEG1000 就是 '
                 + '、'.join(deg1000_genes)
                 + (' 这 1 个基因' if deg1000 == 1 else f' 这 {deg1000} 个基因'))
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
    widths = ([1.6, 2.1, 3.0, 3.8, 1.4, 3.8, 1.3] if tuple(columns) == SOURCE_TABLE_COLUMNS
              else [2.2, 1.8, 1.2, 1.6, 1.6, 1.6, 8.0] if tuple(columns) == DISEASE_TABLE_COLUMNS
              else [2.0, 5.4, 3.2, 3.2, 3.2] if tuple(columns) == COVERAGE_TABLE_COLUMNS
              else [6.0, 12.0])
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
    from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
    borders = OxmlElement('w:tblBorders')
    for edge in ('top', 'left', 'bottom', 'right', 'insideH', 'insideV'):
        element = OxmlElement('w:' + edge)
        for key, value in (('val', 'single'), ('sz', '4'), ('color', 'D9D9D9')):
            element.set(qn('w:' + key), value)
        borders.append(element)
    grid._tbl.tblPr.append(borders)
    for ri, row in enumerate(grid.rows):
        for ci, cell in enumerate(row.cells):
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            margins = OxmlElement('w:tcMar')
            for edge in ('top', 'left', 'bottom', 'right'):
                margin = OxmlElement('w:' + edge)
                margin.set(qn('w:w'), '35' if tuple(columns) == COVERAGE_TABLE_COLUMNS else '65')
                margin.set(qn('w:type'), 'dxa')
                margins.append(margin)
            cell._tc.get_or_add_tcPr().append(margins)
            if ri == 0:
                shading = OxmlElement('w:shd'); shading.set(qn('w:fill'), 'DCE6F1')
                cell._tc.get_or_add_tcPr().append(shading)
            if ri == 0 or ci in (4, 6) and tuple(columns) == SOURCE_TABLE_COLUMNS:
                cell.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER
    header = OxmlElement('w:tblHeader')
    header.set(qn('w:val'), '1')
    grid.rows[0]._tr.get_or_add_trPr().append(header)
    for row in grid.rows:
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                paragraph.paragraph_format.space_after = Pt(0.5 if tuple(columns) == COVERAGE_TABLE_COLUMNS else 2)
                paragraph.paragraph_format.space_before = Pt(0.5 if tuple(columns) == COVERAGE_TABLE_COLUMNS else 2)
                paragraph.paragraph_format.line_spacing = Pt(12) if tuple(columns) == COVERAGE_TABLE_COLUMNS else 1.0
        avoid_split = OxmlElement('w:cantSplit')
        row._tr.get_or_add_trPr().append(avoid_split)
        avoid_split.set(qn('w:val'), '1')
    return grid


def _percent(value):
    return f'{value:.2f}%'


def coverage_note(disease_rows):
    """The coverage table's reading rule, with the denominators read off table 2."""
    denominators = '、'.join(f"{row['disease']} {row['DEG1000']}" for row in disease_rows if row['eligible'] == 'yes')
    return (
        '每行 = 一个组（tenv-HVG 或一个疾病）在一个轴 N 上的覆盖率。tenv-HVG 行沿用 v1 口径：'
        '分子 = 其 top-N 集合中落进 G_N 的基因数，分母 = N。疾病行是 v2 口径：'
        f'分子 = |DEG1000_d ∩ G_N|，分母 = |DEG1000_d|（{denominators}），'
        '排除疾病的上游证据仍保留，覆盖率三项留空并附理由。各行百分比不能跨组比较。'
    )


def _check_coverage_denominators(table, disease_rows):
    """The coverage table's disease denominators must be each disease's ``DEG1000``.

    tenv-HVG keeps ``N`` (its top-``N`` set), every disease divides by its own
    ``|DEG1000_d|``. Reading a disease row against ``N`` is the v1 misreading the
    long-format contract exists to prevent, so the document refuses to render a
    table whose two halves disagree.
    """
    published = {row['disease']: int(row['DEG1000'].replace(',', ''))
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

    display_labels = gene_axis.CATEGORY_DISPLAY_LABELS
    for row in disease_rows:
        row['disease'] = display_labels.get(row['disease'], row['disease'])

    document = Document()
    section = document.sections[0]
    section.page_width = Cm(21)
    section.page_height = Cm(29.7)
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
    eligible = [row['disease'] for row in disease_rows if row['eligible'] == 'yes']
    excluded = [row['disease'] for row in disease_rows if row['eligible'] == 'no']
    _body(document,
          '依据 ADR 0073 完整收集、tier 优先 DEG1000 与 f1000 规则，说明四档嵌套基因轴及实际来源覆盖。'
          '使用 G_2000、G_3000、G_4000 和 G_5000 时，请核对疾病资格、完整并集分母及每档补齐量。'
          f'本次合格疾病为 {"、".join(eligible)}；排除疾病为 {"、".join(excluded) or "无"}。'
          '表格与图来自同一重建版本，规则文字与公开实现一致，可由文末命令重放。')

    _heading(document, RULE_HEADING)
    _body(document, DEG_DEFINITION)
    steps = list(RULE_STEPS)
    steps[5] = steps[5] + fallback_sentence(fallback_slots)
    for index, step in enumerate(steps, start=1):
        document.add_paragraph(f'{index}) {step}')
    _body(document, F1000_DEFINITION, bold_lead='f1000 定义：')

    _heading(document, SOURCE_TABLE_HEADING).paragraph_format.page_break_before = True
    for index, category in enumerate((*(gene_axis.CATEGORY_LABELS[c] for c in gene_axis.RELEASE_CATEGORIES), TENV_HVG_ROW['category']), start=1):
        records = [row for row in source_rows if row['category'] == category]
        if not records:
            raise ValueError(f'source table category {category!r} has no rows')
        caption = document.add_paragraph(f'表 1 {index}　{display_labels.get(category, category)}' + (' 非疾病来源' if category == TENV_HVG_ROW['category'] else ' 来源数据集'))
        caption.alignment = WD_ALIGN_PARAGRAPH.CENTER
        caption.paragraph_format.keep_with_next = True
        _grid(document, SOURCE_TABLE_COLUMNS,
              [[display_labels.get(category, category) if name == 'category' else record[name]
                for name in SOURCE_TABLE_COLUMNS] for record in records],
              font_size=8.5)
        if category != TENV_HVG_ROW['category']:
            catalog_category = next(key for key in gene_axis.RELEASE_CATEGORIES
                                    if gene_axis.CATEGORY_LABELS[key] == category)
            entries = json.loads(Path(catalogs[catalog_category]).read_text())['datasets']
            notes = {}
            for entry in entries:
                evidence = metadata.SOURCES[entry['dataset_id']]
                note = evidence.get('sample_size_note')
                if note:
                    notes.setdefault(note, []).append(evidence['accession'] or entry['dataset_id'])
            for note, accessions in notes.items():
                paragraph = _body(document, '样本量注 ' + '、'.join(dict.fromkeys(accessions)) + '：' + note)
                for run in paragraph.runs:
                    run.font.size = Pt(8.5)
    _body(document, SOURCE_TABLE_FOOTNOTE)

    _heading(document, DISEASE_TABLE_HEADING).paragraph_format.page_break_before = True
    caption = document.add_paragraph(DISEASE_TABLE_CAPTION)
    caption.alignment = WD_ALIGN_PARAGRAPH.CENTER
    _grid(document, DISEASE_TABLE_COLUMNS,
          [[record[name] for name in DISEASE_TABLE_COLUMNS] for record in disease_rows],
          font_size=8.5)
    _body(document, DISEASE_TABLE_FOOTNOTE)

    _heading(document, '4 四处易错点')
    for lead, text in CLARIFICATIONS:
        _body(document, text, bold_lead=lead + ' ')

    _heading(document, '5 四档轴的疾病覆盖率').paragraph_format.page_break_before = True
    _body(document, coverage_note(disease_rows))
    coverage_rows = [
        [str(int(record.N)), display_labels.get(record.group, record.group_label), ('excluded' if pd.isna(record.numerator) else str(int(record.numerator))),
         ('—' if pd.isna(record.denominator) else str(int(record.denominator))), ('—' if pd.isna(record.rate) else _percent(record.rate))]
        for record in table.itertuples(index=False)
    ]
    _grid(document, COVERAGE_TABLE_COLUMNS, coverage_rows, font_size=8.5)

    _heading(document, '6 证据图').paragraph_format.page_break_before = True
    document.add_picture(str(figure_path), width=Cm(15.5))
    caption = document.add_paragraph(
        '图 1　四个子图分别为 N = 2,000、3,000、4,000、5,000。灰柱是全部 '
        '18,164 个冻结轴候选基因的 rank_sum 分布；每根柱自顶端向下的彩色段表示该柱中'
        '被选进轴的基因数量（彩色总量 = N）。')
    caption.alignment = WD_ALIGN_PARAGRAPH.CENTER

    _heading(document, '7 复现命令与相关资产').paragraph_format.page_break_before = True
    for command in REPRODUCE:
        document.add_paragraph(command, style='No Spacing')
    assets = (
        ('candidate_rank_sum.csv',
         '主候选秩和表：18,164 行，含 tenv_hvg_rank、声明疾病秩及资格、rank_sum、f1000'),
        ('gene_axis_2000.csv 至 gene_axis_5000.csv', '四个已发布轴（gene_index, symbol）'),
        ('axis_composition.csv', '每轴的 axis_fill / DEG1000 / f1000 构成记录'),
        ('source_topn_coverage.csv', '本文件的覆盖率表（长表：每组 N、分子、分母、比率）'),
        ('manifest.json', '输入/输出摘要与 sha256（release.algorithm = tier-first-deg1000-f1000-four-axes-v3）'),
        ('disease-fusion/catalog.json', '每病 level-1 融合与 DEG1000 的发布记录（表 2 的来源）'),
        ('figure/data/gene_axis_rank_sum_histogram.csv', '图中全部绘制值'),
        ('execution/gene_axis_source_metadata.py', '逐来源平台、技术、样本量口径与已核验的原始记录链接'),
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
        'out': step / 'gene_axis_fill_rule.docx',
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
