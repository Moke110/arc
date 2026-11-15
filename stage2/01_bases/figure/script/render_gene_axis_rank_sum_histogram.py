#!/usr/bin/env python3
"""Render the released gene-axis rank-sum evidence figure from its plotted-values CSV.

This is a publication script, not scientific code: it reads only
``figure/data/gene_axis_rank_sum_histogram.csv`` and writes the PNG. It never
opens ``output/`` products, matrices, checkpoints or truth, and it recomputes no
metric. The CSV carries every plotted value -- the candidate histogram over all
frozen-axis candidates, one panel per released axis, and how many of each bin's
candidates that axis selects -- so this renderer needs no other input.

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
from matplotlib.patches import Patch  # noqa: E402

# The stage palette is a pure public implementation (ADR 0063): importing it
# reads no product, so the renderer still needs only its plotted-values CSV.
_API_ROOT = Path(__file__).resolve().parents[2] / 'output' / 'api'
if str(_API_ROOT) not in sys.path:
    sys.path.insert(0, str(_API_ROOT))
from bases.palette import FILL_SHADE, LINE_SHADE, color  # noqa: E402

RECORDS = ('figure', 'panel', 'hist')
CANDIDATE_COLOR = color('gray', FILL_SHADE)
SELECTED_COLOR = color('orange', LINE_SHADE)
CANDIDATE_LABEL = 'gene candidates'
SELECTED_LABEL = 'selected genes'
#: The x-axis annotation the released figure carries (2026-11-14 design request).
X_LABEL = 'rank_sum (descending is better)'
Y_LABEL = 'genes'


def load(path):
    """Read and validate the complete plotted-values table.

    Returns ``(title, panels, histograms)``: the figure title, the panel records
    in CSV order (``panel``, ``axis``, ``panel_title``) and, per panel id, the
    bin list ``(left, right, count, selected)``. The renderer refuses a table
    whose panels disagree on their bins, whose bars select more genes than they
    hold, or whose selected totals are not each panel's declared ``N``.
    """
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
    panels = []
    for row in rows:
        if row['record'] != 'panel':
            continue
        panel_id = (row['panel'] or '').strip()
        size = (row['axis'] or '').strip()
        panel_title = (row['panel_title'] or '').strip()
        if not panel_id or not size.isdigit() or int(size) <= 0 or not panel_title:
            raise ValueError('a panel record is malformed')
        if panel_id in {panel['panel'] for panel in panels}:
            raise ValueError('panel ids must be unique')
        panels.append({'panel': panel_id, 'axis': int(size),
                       'panel_title': panel_title})
    if not panels:
        raise ValueError('at least one panel record is required')
    histograms = {}
    for panel in panels:
        bins = []
        for row in rows:
            if row['record'] != 'hist':
                continue
            if (row['panel'] or '').strip() != panel['panel']:
                continue
            left = float(row['bin_left'])
            right = float(row['bin_right'])
            count = int(row['count'])
            selected = int(row['selected'])
            if right <= left:
                raise ValueError('a histogram bin is malformed')
            if count < 0 or selected < 0 or selected > count:
                raise ValueError('a histogram bin has an impossible selected count')
            bins.append((left, right, count, selected))
        if not bins:
            raise ValueError(f"panel {panel['panel']} has no histogram rows")
        histograms[panel['panel']] = bins
    edges = [(left, right) for left, right, _, _ in histograms[panels[0]['panel']]]
    for panel in panels:
        panel_edges = [(left, right)
                       for left, right, _, _ in histograms[panel['panel']]]
        if panel_edges != edges:
            raise ValueError('every panel must share the same candidate bins')
        counts = sum(count for _, _, count, _ in histograms[panel['panel']])
        if counts <= 0:
            raise ValueError('the candidate histogram is empty')
        selected = sum(chosen for _, _, _, chosen in histograms[panel['panel']])
        if selected != panel['axis']:
            raise ValueError(f"panel {panel['panel']} does not select exactly N genes")
    return figures[0]['title'], panels, histograms


def render(title, panels, histograms, output):
    """Draw the multi-panel evidence figure with fixed publication geometry."""
    output = Path(output)
    if output.suffix != '.png':
        raise ValueError('a .png output path is required')
    output.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, len(panels), figsize=(5.6 * len(panels), 4.4),
                                layout='constrained', squeeze=False)
    try:
        for axis, panel in zip(axes[0], panels):
            bins = histograms[panel['panel']]
            left = [item[0] for item in bins]
            width = [item[1] - item[0] for item in bins]
            counts = [item[2] for item in bins]
            selected = [item[3] for item in bins]
            # The gray candidate bar is the panel's base; the colored segment
            # hangs from each bar's top down to the unselected remainder.
            axis.bar(left, counts, width=width, align='edge', color=CANDIDATE_COLOR,
                     edgecolor='none')
            axis.bar(left, selected,
                     bottom=[count - chosen for count, chosen in zip(counts, selected)],
                     width=width, align='edge', color=SELECTED_COLOR, edgecolor='none')
            axis.set(title=panel['panel_title'], xlabel=X_LABEL, ylabel=Y_LABEL)
            axis.legend(handles=[Patch(color=CANDIDATE_COLOR, label=CANDIDATE_LABEL),
                                 Patch(color=SELECTED_COLOR, label=SELECTED_LABEL)],
                        frameon=False, fontsize=8, loc='upper right')
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
