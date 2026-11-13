#!/usr/bin/env python3
"""Render the released gene-axis rank-sum evidence panel from its plotted-values CSV.

This is a publication script, not scientific code: it reads only
``figure/data/gene_axis_rank_sum_histogram.csv`` and writes the PNG. It never
opens ``output/`` products, matrices, checkpoints or truth, and it recomputes no
metric. The CSV carries every plotted value -- the rank-sum histogram bins over
all frozen-axis candidates, both released axes' ``rank_sum`` intervals with their
medians, and the explicit inflection diagnostic marker (or its ``no inflection``
note) -- so this renderer needs no other input.

    python3 stage2/01_bases/figure/script/render_gene_axis_rank_sum_histogram.py \
        --csv stage2/01_bases/figure/data/gene_axis_rank_sum_histogram.csv \
        --png stage2/01_bases/figure/gene_axis_rank_sum_histogram.png

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

# The stage palette is a pure public implementation (ADR 0063): importing it
# reads no product, so the renderer still needs only its plotted-values CSV.
_API_ROOT = Path(__file__).resolve().parents[2] / 'output' / 'api'
if str(_API_ROOT) not in sys.path:
    sys.path.insert(0, str(_API_ROOT))
from bases.palette import FILL_SHADE, LINE_SHADE, color  # noqa: E402

RECORDS = ('figure', 'panel', 'hist', 'interval', 'marker')
HISTOGRAM_FILL = color('gray', FILL_SHADE)
INTERVAL_COLORS = (color('blue', FILL_SHADE), color('orange', FILL_SHADE))
INTERVAL_LINES = (color('blue', LINE_SHADE), color('orange', LINE_SHADE))


def load(path):
    """Read and validate the complete plotted-values table."""
    with open(path, newline='') as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError('plotted-values CSV is empty')
    kinds = {row['record'] for row in rows}
    if not set(RECORDS) <= kinds:
        raise ValueError(f'plotted-values CSV is incomplete: {sorted(kinds)}')
    figures = [row for row in rows if row['record'] == 'figure']
    if len(figures) != 1 or not figures[0]['title']:
        raise ValueError('exactly one figure record with a title is required')
    panels = [row for row in rows if row['record'] == 'panel']
    if len(panels) != 1 or not panels[0]['panel_title']:
        raise ValueError('exactly one panel record with a title is required')
    histogram = []
    for row in rows:
        if row['record'] != 'hist':
            continue
        left, right, count = float(row['bin_left']), float(row['bin_right']), int(row['count'])
        if right <= left or count < 0:
            raise ValueError('histogram bin is malformed')
        histogram.append((left, right, count))
    if not histogram:
        raise ValueError('the rank-sum histogram values are missing')
    intervals = []
    for row in rows:
        if row['record'] != 'interval':
            continue
        low = float(row['interval_low'])
        high = float(row['interval_high'])
        median = float(row['median'])
        if not low <= median <= high:
            raise ValueError('an axis interval does not contain its median')
        intervals.append({'axis': row['axis'], 'low': low, 'high': high,
                          'median': median})
    if len(intervals) != 2 or len({item['axis'] for item in intervals}) != 2:
        raise ValueError('exactly two distinct axis intervals are required')
    markers = [row for row in rows if row['record'] == 'marker']
    if len(markers) != 1:
        raise ValueError('exactly one inflection marker record is required')
    marker_row = markers[0]
    inflection = str(marker_row['inflection']).strip().lower() == 'true'
    marker = {'inflection': inflection,
              'n': int(marker_row['inflection_n']) if marker_row['inflection_n'] else None,
              'note': marker_row['note']}
    if inflection and marker['n'] is None:
        raise ValueError('an accepted inflection needs its candidate N')
    return figures[0]['title'], panels[0]['panel_title'], histogram, intervals, marker


def render(title, panel_title, histogram, intervals, marker, output):
    """Draw the one-panel evidence figure with fixed publication geometry."""
    output = Path(output)
    if output.suffix != '.png':
        raise ValueError('a .png output path is required')
    output.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(9.6, 5.6), layout='constrained')
    try:
        left = [item[0] for item in histogram]
        width = [item[1] - item[0] for item in histogram]
        axis.bar(left, [item[2] for item in histogram], width=width,
                 align='edge', color=HISTOGRAM_FILL, edgecolor='none',
                 label='frozen-axis candidates')
        for item, fill, line in zip(intervals, INTERVAL_COLORS, INTERVAL_LINES):
            axis.axvspan(item['low'], item['high'], color=fill, alpha=0.30,
                         label=f"G_{item['axis']} rank_sum interval")
            axis.axvline(item['median'], color=line, linewidth=1.4, linestyle='--')
        axis.set(xlabel='rank_sum (ascending is better)',
                 ylabel='genes', title=panel_title)
        axis.grid(alpha=0.18, axis='y')
        axis.legend(frameon=False, fontsize=8, loc='upper right')
        axis.text(0.02, 0.95, marker['note'], transform=axis.transAxes,
                  va='top', fontsize=8.5)
        figure.suptitle(title, fontsize=11.5)
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
