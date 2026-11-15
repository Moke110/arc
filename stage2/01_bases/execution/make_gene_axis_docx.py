#!/usr/bin/env python3
"""Write the human-facing gene-axis rule document from the released products.

The document is a reading aid, not a release product: it restates the fill rule
with its common misreadings, carries the source top-N coverage table and embeds
the released evidence figure. Its inputs are the published coverage CSV
(``output/gene-axis/source_topn_coverage.csv``) and the released PNG; the rule
text is a constant here and mirrors ``bases.gene_axis.FILL_RULE`` / the README.

    python stage2/01_bases/execution/make_gene_axis_docx.py \
        --coverage stage2/01_bases/output/gene-axis/source_topn_coverage.csv \
        --figure stage2/01_bases/figure/gene_axis_rank_sum_histogram.png \
        --out stage2/01_bases/output/gene-axis/gene_axis_fill_rule.docx

The DOCX container is re-zipped with fixed entry timestamps, so two runs from the
same inputs and the same ``python-docx`` produce identical bytes.
"""
import argparse
import io
import zipfile
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Cm, Pt

TITLE = 'Stage 2 基因轴取用规则与来源覆盖（G_2000 / G_5000）'
COVERAGE_COLUMNS = ('N', '%tenv_hvg_gN', '%nhf_gN', '%icm_gN',
                    '%genetic_cardiomyopathy_gN', '%cardiotox_TOP2_gN',
                    '%cardiotox_other_gN')

RULE_STEPS = (
    '第一步：取 tenv_hvg_rank ≤ N/2 的基因，即 tenv-HVG 秩数值最低的 0.5N 个。'
    'tenv-HVG 秩是 1..18,164 的完整排列，所以这一步恰好取满 0.5N 个。',
    '第二步：剩余候选 = tenv_hvg_rank > N/2（必须排除第一步已选的基因）。'
    '把剩余候选按 pN 降序 → rank_sum 升序 → gene_index 升序排序，取前 0.5N 个。',
    '合并两半，按 gene_index 升序发布。轴是冻结 18,164 轴的有序子集，保留原始 '
    'gene_index 整数，绝不重编号为 0..N-1。',
)

PN_DEFINITION = (
    '在 6 个融合源（5 个疾病类别 + tenv-HVG）里，把该基因放进各自 top-N 集合的源'
    '的数量。N 与轴大小对应：G_2000 用 p2000，G_5000 用 p5000。'
)

CLARIFICATIONS = (
    ('pN 不是 n_sources_measured。',
     'n_sources_measured 是“测到该基因的源数”，两者只在 2/18,164 行（p2000）与 '
     '96/18,164 行（p5000）恰好相等；用 n_sources_measured 代替 pN 会把 G_2000 换掉 '
     '289 个基因、G_5000 换掉 425 个。'),
    ('rank_sum 越小越好（数值下降方向更好）。',
     'rank_sum 是各源 rank_norm 之和；按降序（数值越大越靠前）取会换掉 G_2000 的 300 个'
     '与 G_5000 的 1,004 个基因。'),
    ('剩余候选必须排除第一步已选。',
     '若对全部 18,164 个候选重排再取，会重复抽到已选基因，50/50 结构立即破坏。'),
    ('gene_index 升序破平是契约的一部分。',
     'rank_sum 由整数秩相加，重复很多：剩余候选中 N=2,000 有 2,038 个重复值覆盖 '
     '4,341 个基因，N=5,000 有 1,776 个覆盖 3,761 个。已发布的两个 N 其边界恰好落在'
     '唯一 (pN, rank_sum) 键上，所以现有两轴对该破平不敏感；换 N 或换排序实现时不是。'
     '发布表行序本就是 gene_index 升序，因此稳定排序恰好等价于显式破平。'),
    ('两轴不嵌套。',
     'pN 随 N 变化，G_2000 ∩ G_5000 = 1,972（Jaccard 0.3922）；98.60% 的 G_2000 落在 '
     'G_5000 内，两轴各自都含有对方没有的基因。'),
)

COVERAGE_NOTE = (
    '每格 = 该来源 top-N 集合中落进 G_N 的百分比。两个集合都恰有 N 个基因，所以这个'
    '计数同时是“该来源 top-N 被轴保留的份额”与“轴中来自该来源 top-N 的份额”。'
    '各来源 top-N 集合彼此重叠（pN 就是在数重叠），因此一行不等于 100%。'
    '参照：随机 N 基因集合覆盖轴的比例为 11.01%（G_2000）与 27.53%（G_5000）。'
)

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


def read_coverage(path):
    """Read and validate the released coverage table (one row per axis)."""
    table = pd.read_csv(path)
    if list(table.columns) != list(COVERAGE_COLUMNS):
        raise ValueError(f'{path}: coverage columns differ from the contract')
    if len(table) == 0:
        raise ValueError(f'{path}: coverage table is empty')
    return table


def _percent(value):
    return '—（排除）' if pd.isna(value) else f'{value:.2f}%'


def _heading(document, text, level=1):
    heading = document.add_heading(text, level=level)
    return heading


def _body(document, text, bold_lead=None):
    paragraph = document.add_paragraph()
    if bold_lead:
        run = paragraph.add_run(bold_lead)
        run.bold = True
    paragraph.add_run(text)
    return paragraph


def build(coverage_path, figure_path, output_path):
    table = read_coverage(coverage_path)
    document = Document()
    section = document.sections[0]
    section.left_margin = Cm(2.2)
    section.right_margin = Cm(2.2)
    style = document.styles['Normal']
    style.font.size = Pt(10.5)

    title = document.add_heading(TITLE, level=0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    _body(document,
          '依据 ADR 0066（发布契约）与 ADR 0068（融合修正）；产品目录 '
          'stage2/01_bases/output/gene-axis/。本文件是人类阅读文档，不是发布产品：'
          '规则文字是常量并与 bases.gene_axis.FILL_RULE、README 同步，表格与图来自'
          '已发布产品，可由文末命令重放。')

    _heading(document, '1. 取用规则（精确版本）')
    for index, step in enumerate(RULE_STEPS, start=1):
        document.add_paragraph(f'{index}) {step}')
    _body(document, PN_DEFINITION, bold_lead='pN 定义：')

    _heading(document, '2. 五处易错点')
    for lead, text in CLARIFICATIONS:
        _body(document, text, bold_lead=lead + ' ')

    _heading(document, '3. 每个来源 top-N 基因的覆盖率')
    _body(document, COVERAGE_NOTE)
    columns = list(table.columns)
    grid = document.add_table(rows=1 + len(table), cols=len(columns))
    grid.style = 'Table Grid'
    for column, name in enumerate(columns):
        cell = grid.cell(0, column)
        cell.text = ''
        run = cell.paragraphs[0].add_run(name)
        run.bold = True
        run.font.size = Pt(8.5)
    for row_index, record in enumerate(table.itertuples(index=False), start=1):
        values = [str(int(record.N))] + [_percent(value) for value in record[1:]]
        for column, value in enumerate(values):
            cell = grid.cell(row_index, column)
            cell.text = ''
            run = cell.paragraphs[0].add_run(value)
            run.font.size = Pt(8.5)

    _heading(document, '4. 证据图')
    document.add_picture(str(figure_path), width=Cm(15.5))
    caption = document.add_paragraph(
        '图 1　每个子图对应一个基因轴，标题为 N = 2,000 / N = 5,000。灰柱是全部 '
        '18,164 个冻结轴候选基因的 rank_sum 分布；每根柱自顶端向下的彩色段表示该柱中'
        '被选进轴的基因数量（彩色总量 = N）。')
    caption.alignment = WD_ALIGN_PARAGRAPH.CENTER

    _heading(document, '5. 复现命令与相关资产')
    for command in REPRODUCE:
        document.add_paragraph(command, style='No Spacing')
    assets = (
        ('candidate_rank_sum.csv', '主候选秩和表：18,164 行，含各源秩、rank_sum、pN'),
        ('gene_axis_2000.csv / gene_axis_5000.csv', '两个已发布轴（gene_index, symbol）'),
        ('axis_composition.csv', '各来源在轴内测到的基因数与 pN 分布'),
        ('source_topn_coverage.csv', '本文件的覆盖率表'),
        ('manifest.json', '输入/输出摘要与 sha256'),
        ('figure/data/gene_axis_rank_sum_histogram.csv', '图中全部绘制值'),
        ('figure/gene_axis_rank_sum_histogram.png', '证据图'),
    )
    assets_table = document.add_table(rows=1 + len(assets), cols=2)
    assets_table.style = 'Table Grid'
    for column, name in enumerate(('资产', '内容')):
        cell = assets_table.cell(0, column)
        cell.text = ''
        run = cell.paragraphs[0].add_run(name)
        run.bold = True
        run.font.size = Pt(9)
    for row_index, (name, detail) in enumerate(assets, start=1):
        for column, value in enumerate((name, detail)):
            cell = assets_table.cell(row_index, column)
            cell.text = ''
            run = cell.paragraphs[0].add_run(value)
            run.font.size = Pt(9)

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
        'out': step / 'output/gene-axis/gene_axis_fill_rule.docx',
    }


def main():
    defaults = default_paths()
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--coverage', type=Path, default=defaults['coverage'])
    parser.add_argument('--figure', type=Path, default=defaults['figure'])
    parser.add_argument('--out', type=Path, default=defaults['out'])
    args = parser.parse_args()
    path = build(args.coverage, args.figure, args.out)
    print(f'wrote {path} from {args.coverage} and {args.figure}')


if __name__ == '__main__':
    main()
