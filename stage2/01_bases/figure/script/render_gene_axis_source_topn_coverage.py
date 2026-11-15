#!/usr/bin/env python3
"""Render the released two-panel source-coverage figure from its plotted-values CSV.

This is a publication script, not scientific code: it reads only
``figure/data/gene_axis_source_topn_coverage.csv`` and writes the PNG. It never
opens ``output/`` products, matrices, checkpoints or truth, and it recomputes no
metric. The CSV carries every plotted value -- the two panel records (id, title,
y-axis label and the value column the panel draws) and, for each released axis
(the bar series) and fused source (the x group), the absolute number of that
source's top-``N`` genes the axis contains and its rate. Panel ``a`` is the count
and panel ``b`` the rate; the bars, colours and group layout are identical.

    python3 stage2/01_bases/figure/script/render_gene_axis_source_topn_coverage.py \
        --csv stage2/01_bases/figure/data/gene_axis_source_topn_coverage.csv \
        --png stage2/01_bases/figure/gene_axis_source_topn_coverage.png

Re-rendering from the CSV reproduces the committed PNG byte for byte.
"""
import argparse
import csv
import os
import sys
import tempfile
from pathlib import Path

# Keep matplotlib's cache out of a read-only home directory; the cache location
# does not affect the plotted bytes.
os.environ.setdefault('MPLCONFIGDIR',
                      str(Path(tempfile.gettempdir()) / 'arc-matplotlib'))

import matplotlib  # noqa: E402  (backend must be set before pyplot)

matplotlib.use('Agg')
from matplotlib import pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# The stage palette is a pure public implementation (ADR 0063): importing it
# reads no product, so the renderer still needs only its plotted-values CSV.
_API_ROOT = Path(__file__).resolve().parents[2] / 'output' / 'api'
if str(_API_ROOT) not in sys.path:
    sys.path.insert(0, str(_API_ROOT))
from bases.palette import FILL_SHADE, LINE_SHADE, color  # noqa: E402

RECORDS = ('panel', 'bar')
#: The numeric columns a panel may draw; every bar carries both.
VALUE_COLUMNS = ('count', 'rate')
#: One palette family per released axis, in bar-series order.
SERIES_FAMILIES = ('blue', 'orange')
#: Bars inside a group touch (no gap); groups keep a 0.2-wide gap.
GROUP_WIDTH = 0.8


def _number(row, name):
    try:
        return float((row[name] or '').strip())
    except ValueError as error:
        raise ValueError(f'a bar {name} is not a number') from error


def load(path):
    """Read and validate the complete plotted-values table.

    Returns ``(panels, series, groups, values)``: the panel records in CSV order
    as ``{'panel', 'title', 'y_label', 'value'}``, the bar series in CSV order as
    ``(name, label)``, the x groups in CSV order as ``(source, source_label)``
    and ``values[(series, source)] = {'count', 'rate'}``. The renderer refuses a
    table without a titled panel, a duplicated (series, source) cell, a missing
    cell, a value column a panel cannot draw, or a count that is not a
    non-negative integer / a rate outside 0..100.
    """
    with open(path, newline='') as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError('plotted-values CSV is empty')
    kinds = {row['record'] for row in rows}
    if not set(RECORDS) <= kinds:
        raise ValueError(f'plotted-values CSV is incomplete: {sorted(kinds)}')
    panels = []
    for row in rows:
        if row['record'] != 'panel':
            continue
        panel = {'panel': (row['panel'] or '').strip(),
                 'title': (row['panel_title'] or '').strip(),
                 'y_label': (row['y_label'] or '').strip(),
                 'value': (row['value'] or '').strip()}
        if not all(panel.values()):
            raise ValueError('a panel record is malformed')
        if panel['value'] not in VALUE_COLUMNS:
            raise ValueError(f"a panel draws an unknown column: {panel['value']}")
        if panel['panel'] in {item['panel'] for item in panels}:
            raise ValueError('panel ids must be unique')
        panels.append(panel)
    if not panels:
        raise ValueError('at least one panel record is required')
    series, groups, values = [], [], {}
    for row in rows:
        if row['record'] != 'bar':
            continue
        name = (row['series'] or '').strip()
        label = (row['series_label'] or '').strip()
        source = (row['source'] or '').strip()
        source_label = (row['source_label'] or '').strip()
        if not name or not label or not source or not source_label:
            raise ValueError('a bar record is malformed')
        count = _number(row, 'count')
        rate = _number(row, 'rate')
        if count < 0 or count != int(count):
            raise ValueError('a bar count must be a non-negative integer')
        if not 0.0 <= rate <= 100.0:
            raise ValueError('a bar rate must be in 0..100 percent')
        if name not in {item[0] for item in series}:
            series.append((name, label))
        if source not in {item[0] for item in groups}:
            groups.append((source, source_label))
        if (name, source) in values:
            raise ValueError('a (series, source) bar is repeated')
        values[(name, source)] = {'count': int(count), 'rate': rate}
    if not series or not groups:
        raise ValueError('the figure needs at least one series and one group')
    expected = {(name, source) for name, _ in series for source, _ in groups}
    if set(values) != expected:
        raise ValueError('the bar records do not cover every (series, source) cell')
    return panels, series, groups, values


def render(panels, series, groups, values, output):
    """Draw the grouped bars once per panel with fixed publication geometry."""
    output = Path(output)
    if output.suffix != '.png':
        raise ValueError('a .png output path is required')
    output.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, len(panels), figsize=(5.6 * len(panels), 4.2),
                                layout='constrained', squeeze=False)
    try:
        width = GROUP_WIDTH / len(series)
        centres = [float(index) for index in range(len(groups))]
        for axis, panel in zip(axes[0], panels):
            heights = [[values[(name, source)][panel['value']]
                        for source, _ in groups] for name, _ in series]
            top = max(max(row) for row in heights) or 1
            handles = []
            for index, (name, label) in enumerate(series):
                family = SERIES_FAMILIES[index % len(SERIES_FAMILIES)]
                fill = color(family, FILL_SHADE)
                edge = color(family, LINE_SHADE)
                left = [centre - GROUP_WIDTH / 2 + (index + 0.5) * width
                        for centre in centres]
                axis.bar(left, heights[index], width=width, color=fill,
                         edgecolor=edge, linewidth=0.6)
                handles.append(Patch(facecolor=fill, edgecolor=edge, label=label))
            axis.set_xticks(centres)
            axis.set_xticklabels([label for _, label in groups], rotation=30,
                                 ha='right', fontsize=8)
            axis.set(ylabel=panel['y_label'], title=panel['title'])
            axis.set_ylim(0, top * 1.12)
            axis.legend(handles=handles, frameon=False, fontsize=8,
                        loc='upper right')
        figure.savefig(output, dpi=160)
    finally:
        plt.close(figure)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--csv', type=Path, required=True,
                        help='figure/data plotted values')
    parser.add_argument('--png', type=Path, required=True, help='figure PNG to write')
    args = parser.parse_args()
    render(*load(args.csv), args.png)
    print(f'rendered {args.png} from {args.csv}')


if __name__ == '__main__':
    main()
