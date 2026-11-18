"""Render the paired-axis Round 1 figures strictly from plotting CSVs."""
import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np

SOURCES = ('tahoe_bulk', 'tenv')
SOURCE_LABELS = dict(tahoe_bulk='Tahoe', tenv='tenv')
SCENARIOS = ('cold-drug', 'cold-drug-dose')
METHODS = ('no-dose', 'scaling', 'input-numeric', 'latent-numeric',
           'input-embedding', 'latent-embedding', 'film')
LABELS = ('No dose', 'Scaling', 'Input numeric', 'Latent numeric',
          'Input embedding', 'Latent embedding', 'FiLM')
METRICS = ('pcc', 'spearman', 'rmse', 'r2')
TITLES = ('Δ PCC ↑', 'Δ Spearman ↑', 'Δ RMSE ↓', 'Δ R² ↑')
COLORS = ('#465970', '#9c6032', '#167c94', '#ad4a5c', '#387048', '#705ba7', '#a45f00')
TRAIN_COLORS = ('#244873', '#6885a6', '#95adbd')
VALID_COLORS = ('#236e60', '#618b51', '#869344')
PNG_NAMES = ('round1_metrics.png', 'round1_learning_curves.png', 'round1_metrics_cold_drug_dose.png')


def read_rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def number(value):
    return float(value) if value not in ('', None) else float('nan')


def segment(lower, upper, start, stop, ticks, scale='linear'):
    return dict(lower=lower, upper=upper, start=start, stop=stop,
                ticks=list(map(float, ticks)), scale=scale)


def linear_segment(lower, upper, start=0., stop=1., *, show_zero=False):
    ticks = MaxNLocator(nbins=4).tick_values(lower, upper)
    ticks = [tick for tick in ticks if lower <= tick <= upper]
    if show_zero:
        ticks = sorted(set(ticks) | {0.})
    return segment(lower, upper, start, stop, ticks)


class AxisScale:
    """Piecewise placement of linear/log segments in a single scientific panel."""

    def __init__(self, segments):
        self.segments = segments

    @classmethod
    def choose(cls, values, *, logarithmic=False):
        values = [float(value) for value in values if math.isfinite(value)]
        if logarithmic and any(value > 0 for value in values):
            if any(value < 0 for value in values):
                raise ValueError('negative values cannot appear in the logarithmic MSE overview')
            positive = [value for value in values if value > 0]
            lo, hi = math.log10(min(positive)), math.log10(max(positive))
            pad = max((hi - lo) * .08, .08)
            start, stop = lo - pad, hi + pad
            ticks = [10. ** exponent for exponent in range(math.ceil(start), math.floor(stop) + 1)]
            if len(ticks) < 2:
                ticks = [float(f'{10 ** exponent:.3g}') for exponent in np.linspace(lo, hi, 3)]
                ticks = sorted(set(ticks))
            return cls([segment(0., 0., .035, .035, [0.]),
                        segment(10. ** start, 10. ** stop, .15, 1., ticks, 'log')])
        if not values or all(value == 0 for value in values):
            return cls([linear_segment(-.05, .05, show_zero=True)])
        lower, upper = min(values), max(values)
        pad = max((upper - lower) * .07, max(abs(lower), abs(upper)) * .015)
        data_lower, data_upper = lower - pad, upper + pad
        width = data_upper - data_lower
        if data_lower > 1.5 * width:
            return cls([segment(0., 0., .035, .035, [0.]),
                        linear_segment(data_lower, data_upper, .15, 1.)])
        if data_upper < -1.5 * width:
            return cls([linear_segment(data_lower, data_upper, 0., .85),
                        segment(0., 0., .965, .965, [0.])])
        return cls([linear_segment(min(data_lower, -pad), max(data_upper, pad), show_zero=True)])

    def map(self, values):
        positions = []
        for value in values:
            position = float('nan')
            if math.isfinite(value):
                for part in self.segments:
                    lower, upper = part['lower'], part['upper']
                    if lower <= value <= upper:
                        if lower == upper:
                            position = part['start']
                        else:
                            if part['scale'] == 'log':
                                value, lower, upper = map(math.log, (value, lower, upper))
                            position = part['start'] + (value - lower) / (upper - lower) * (part['stop'] - part['start'])
                        break
                if math.isnan(position):
                    raise ValueError(f'plot coordinate outside saved axis segments: {value}')
            positions.append(position)
        return positions

    def apply(self, ax, *, side='left', zero_line=True, color='#586275'):
        ticks = [tick for part in self.segments for tick in part['ticks']]
        ax.set_ylim(0, 1)
        ax.set_yticks(self.map(ticks), [f'{tick:.3g}' for tick in ticks])
        ax.tick_params(axis='y', labelsize=8, colors=color)
        ax.spines[side].set_color(color)
        if zero_line:
            ax.axhline(self.map([0])[0], color='.58', lw=.65, zorder=0)
        for lower, upper in zip(self.segments, self.segments[1:]):
            gap = (lower['stop'] + upper['start']) / 2
            x = 0 if side == 'left' else 1
            ax.spines[side].set_visible(False)
            ax.plot([x, x], [0, lower['stop']], transform=ax.transAxes, color=color, lw=.8, clip_on=False)
            ax.plot([x, x], [upper['start'], 1], transform=ax.transAxes, color=color, lw=.8, clip_on=False)
            for offset in (-.013, .013):
                ax.plot([x - .014, x + .014], [gap + offset - .009, gap + offset + .009],
                        transform=ax.transAxes, color=color, lw=.9, clip_on=False)

    def plot(self, ax, x, values, *args, label=None, **kwargs):
        """Draw only within each saved axis segment, retaining isolated points."""
        mapped = self.map(values)
        lines = []
        for part in self.segments:
            segment_y = [position if part['lower'] <= value <= part['upper'] else float('nan')
                         for value, position in zip(values, mapped, strict=True)]
            present = [math.isfinite(value) for value in segment_y]
            if not any(present):
                continue
            isolated = [index for index, found in enumerate(present) if found
                        and not (index > 0 and present[index - 1])
                        and not (index + 1 < len(present) and present[index + 1])]
            lines.extend(ax.plot(x, segment_y, *args, marker='.', markevery=isolated,
                                 label=label if not lines else None, **kwargs))
        return lines

    def bar(self, ax, x, value, color):
        for part in self.segments:
            lo = max(min(0, value), part['lower'])
            hi = min(max(0, value), part['upper'])
            if lo < hi:
                bottom, top = self.map([lo, hi])
                ax.bar(x, top - bottom, bottom=bottom, width=.68, color=color, alpha=.28)


def write_rows(path, table):
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)


def by_key(rows, *fields):
    result = defaultdict(list)
    for row in rows:
        result[tuple(row[field] for field in fields)].append(row)
    return result


def overview_range(points, spreads, metric):
    values = [number(row[metric]) for row in points]
    for row in spreads:
        mean, sd = number(row['mean']), number(row['sd'])
        values.extend((mean, mean - sd, mean + sd))
    return values


def export_axes(tables):
    folds = by_key(tables['round1_runs.csv'], 'source', 'scenario')
    spreads = by_key(tables['round1_summary.csv'], 'source', 'scenario', 'metric')
    curves = by_key(tables['round1_curves.csv'], 'source', 'scenario', 'method')
    doses = by_key(tables['round1_dose.csv'], 'source')
    axes = []

    def record(figure, source, scenario, panel, name, values, *, logarithmic=False):
        axis = AxisScale.choose(values, logarithmic=logarithmic)
        axes.append(dict(figure=figure, source=source, scenario=scenario, panel=panel,
                         ordinate=name, segments=json.dumps(axis.segments, separators=(',', ':'))))

    for source in SOURCES:
        for scenario in SCENARIOS:
            for metric in METRICS:
                spreads_for_metric = spreads[source, scenario, metric]
                record('overall', source, scenario, metric, 'metric',
                       overview_range(folds[source, scenario], spreads_for_metric, metric))
            all_mse = []
            for method in METHODS:
                rows = curves[source, scenario, method]
                losses = [number(row['train_loss']) for row in rows]
                mse = [number(row['validation_mse']) for row in rows]
                all_mse.extend(mse)
                record('learning', source, scenario, method, 'loss', losses)
                record('learning', source, scenario, method, 'mse', mse)
            record('learning', source, scenario, 'all-methods', 'mse', all_mse, logarithmic=True)
        for metric in METRICS:
            record('dose', source, 'cold-drug', metric, 'metric',
                   [number(row[metric]) for row in doses[source,]])
    return axes


def validate_coordinates(tables):
    run_keys = {(source, scenario, method, str(fold), '42')
                for source in SOURCES for scenario in SCENARIOS
                for method in METHODS for fold in range(3)}
    spread_keys = {(source, scenario, method, metric, '42')
                   for source in SOURCES for scenario in SCENARIOS
                   for method in METHODS for metric in METRICS}
    dose_keys = {(source, 'cold-drug', method, str(dose), '42')
                 for source, doses in zip(SOURCES, ((.05, .5, 5.), (.2, 1., 10.)))
                 for method in METHODS for dose in doses}
    for filename, keys, expected in [
        ('round1_runs.csv', ('source', 'scenario', 'method', 'fold', 'seed'), run_keys),
        ('round1_summary.csv', ('source', 'scenario', 'method', 'metric', 'seed'), spread_keys),
        ('round1_dose.csv', ('source', 'scenario', 'method', 'dose_um', 'seed'), dose_keys),
    ]:
        rows = tables[filename]
        if filename == 'round1_dose.csv':
            rows = [dict(row, dose_um=str(float(row['dose_um']))) for row in rows]
        index = by_key(rows, *keys)
        if set(index) != expected or any(len(group) != 1 for group in index.values()):
            raise ValueError(f'incomplete or duplicate plot coordinates: {filename}')
    learning = by_key(tables['round1_curves.csv'], 'source', 'scenario', 'method', 'fold', 'seed')
    if set(learning) != run_keys:
        raise ValueError('missing or unexpected learning curves')
    for key, rows in learning.items():
        if sorted(int(row['epoch']) for row in rows) != list(range(1, len(rows) + 1)):
            raise ValueError(f'incomplete epoch coordinates: {key}')


def prepare(reports, destination):
    """Freeze only report tables and explicit axis metadata into plotting CSVs."""
    destination = Path(destination)
    reports = Path(reports)
    tables = {}
    for original, plot in [('folds.csv', 'round1_runs.csv'),
                           ('fold_spread.csv', 'round1_summary.csv'),
                           ('learning.csv', 'round1_curves.csv'),
                           ('per_dose.csv', 'round1_dose.csv')]:
        relevant = ('cold-drug',) if original == 'per_dose.csv' else SCENARIOS
        tables[plot] = [row for scenario in relevant for row in read_rows(reports / scenario / original)]
    validate_coordinates(tables)
    tables['round1_axis_segments.csv'] = export_axes(tables)
    destination.mkdir(parents=True, exist_ok=False)
    for filename, table in tables.items():
        write_rows(destination / filename, table)


def save_figure(figure, out, name):
    figure.savefig(out / name, dpi=130, metadata={'Software': 'ARC axis-local round1 CSV renderer'})
    plt.close(figure)


def source_label(ax, source, scenario=None):
    label = SOURCE_LABELS[source]
    if scenario:
        label += f'  ·  {scenario}'
    ax.annotate(label, xy=(0, 1), xytext=(0, 29), xycoords='axes fraction',
                textcoords='offset points', fontsize=12, weight='bold', annotation_clip=False)


def draw_overall(data, axes, out):
    points = by_key(read_rows(data / 'round1_runs.csv'), 'source', 'scenario', 'method')
    spreads = by_key(read_rows(data / 'round1_summary.csv'), 'source', 'scenario', 'method', 'metric')
    figure, panels = plt.subplots(4, 4, figsize=(22, 16))
    for i, (source, scenario) in enumerate((s, t) for s in SOURCES for t in SCENARIOS):
        source_label(panels[i, 0], source, scenario)
        for ax, metric, title in zip(panels[i], METRICS, TITLES):
            scale = axes['overall', source, scenario, metric, 'metric']
            scale.apply(ax)
            names = []
            for x, (method, label, color) in enumerate(zip(METHODS, LABELS, COLORS, strict=True)):
                (spread,) = spreads[source, scenario, method, metric]
                mean, sd = number(spread['mean']), number(spread['sd'])
                if math.isfinite(mean):
                    scale.bar(ax, x, mean, color)
                if math.isfinite(sd):
                    lo, mid, hi = scale.map([mean - sd, mean, mean + sd])
                    ax.errorbar(x, mid, yerr=[[mid - lo], [hi - mid]],
                                color=color, capsize=3, lw=1.5)
                folds = sorted(points[source, scenario, method], key=lambda row: int(row['fold']))
                ax.scatter([x + delta for delta in (-.13, 0, .13)],
                           scale.map([number(row[metric]) for row in folds]),
                           color=color, s=20, edgecolors='white', linewidths=.3, zorder=4)
                names.append(f"{label}\n{spread['n_defined_folds']}/3 folds")
            ax.set_title(title)
            ax.set_xticks(range(len(METHODS)), names, rotation=32, ha='right', fontsize=8)
            ax.set_xlim(-.55, len(METHODS) - .45)
            ax.grid(axis='y', alpha=.15)
    figure.suptitle('Round 1 · equal-fold test means ± sample SD', fontsize=19)
    figure.text(.5, .022, 'Each dot is one outer fold. Bars weight folds equally; scientific pooled tables weight samples equally. '
                'SD describes fold spread and is neither a confidence interval nor biological replication.',
                ha='center', fontsize=11)
    figure.subplots_adjust(left=.06, right=.985, top=.91, bottom=.105, hspace=1.1, wspace=.27)
    save_figure(figure, out, PNG_NAMES[0])


def draw_learning(data, axes, out):
    curves = by_key(read_rows(data / 'round1_curves.csv'), 'source', 'scenario', 'method', 'fold')
    figure, panels = plt.subplots(4, 8, figsize=(42, 17), sharex=False)
    for i, (source, scenario) in enumerate((s, t) for s in SOURCES for t in SCENARIOS):
        source_label(panels[i, 0], source, scenario)
        overview = panels[i, -1]
        log_scale = axes['learning', source, scenario, 'all-methods', 'mse']
        log_scale.apply(overview, color='#2e4150')
        for ax, method, label, color in zip(panels[i, :-1], METHODS, LABELS, COLORS, strict=True):
            right = ax.twinx()
            mse_scale = axes['learning', source, scenario, method, 'mse']
            loss_scale = axes['learning', source, scenario, method, 'loss']
            loss_scale.apply(ax, color=TRAIN_COLORS[0])
            right.spines['right'].set_visible(True)
            mse_scale.apply(right, side='right', zero_line=False, color=VALID_COLORS[0])
            # The twin axis has its own independently positioned zero reference.
            right.axhline(mse_scale.map([0])[0], color=VALID_COLORS[0], lw=.55, alpha=.55)
            for fold in range(3):
                rows = sorted(curves[source, scenario, method, str(fold)], key=lambda row: int(row['epoch']))
                epochs = [int(row['epoch']) for row in rows]
                losses = [number(row['train_loss']) for row in rows]
                mse = [number(row['validation_mse']) for row in rows]
                loss_scale.plot(ax, epochs, losses, '--', color=TRAIN_COLORS[fold], lw=1)
                mse_scale.plot(right, epochs, mse, color=VALID_COLORS[fold], lw=1)
                log_scale.plot(overview, epochs, mse, color=color, alpha=.55, lw=.95,
                              label=label if fold == 0 else None)
            ax.set_title(label)
            ax.set_ylabel('Native train loss', color=TRAIN_COLORS[0], fontsize=9)
            right.set_ylabel('Validation MSE', color=VALID_COLORS[0], fontsize=9)
            ax.set_xlabel('Epoch')
        overview.set_title('All methods')
        overview.set_xlabel('Epoch')
        overview.set_ylabel('Validation MSE · log data segment', fontsize=9)
        overview.legend(loc='upper right', fontsize=7, framealpha=.86)
    handles = [Line2D([0], [0], color=color, ls=style, label=f'Fold {fold} · {name}')
               for name, colors, style in [('train', TRAIN_COLORS, '--'), ('validation', VALID_COLORS, '-')]
               for fold, color in enumerate(colors)]
    figure.legend(handles=handles, loc='lower center', ncol=6, bbox_to_anchor=(.5, .018))
    figure.suptitle('Round 1 · native training loss and validation MSE by epoch · seed 42', fontsize=21)
    figure.subplots_adjust(left=.038, right=.985, top=.92, bottom=.087, hspace=.78, wspace=.72)
    save_figure(figure, out, PNG_NAMES[1])


def draw_doses(data, axes, out):
    doses = by_key(read_rows(data / 'round1_dose.csv'), 'source', 'method')
    figure, panels = plt.subplots(2, 4, figsize=(22, 10))
    for source, row in zip(SOURCES, panels, strict=True):
        source_label(row[0], source)
        for ax, metric, title in zip(row, METRICS, TITLES):
            scale = axes['dose', source, 'cold-drug', metric, 'metric']
            scale.apply(ax)
            for method, label, color in zip(METHODS, LABELS, COLORS, strict=True):
                values = sorted(doses[source, method], key=lambda record: float(record['dose_um']))
                coordinates = [float(record['dose_um']) for record in values]
                ax.plot(coordinates, scale.map([number(record[metric]) for record in values]),
                        '-o', color=color, label=label, lw=1.5, ms=4)
            ax.set_xscale('log')
            ax.set_xticks(coordinates, [f'{dose:g}' for dose in coordinates])
            ax.set_xlabel('Measured dose (µM)')
            ax.set_title(title)
            ax.grid(axis='y', alpha=.15)
    handles, labels = panels[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='lower center', ncol=7, bbox_to_anchor=(.5, .018))
    figure.suptitle('Round 1 · cold-drug test performance by measured dose', fontsize=19)
    figure.text(.5, .077, 'Test samples have equal weight within each source and dose.',
                ha='center', fontsize=11)
    figure.subplots_adjust(left=.06, right=.985, top=.865, bottom=.15, hspace=.55, wspace=.25)
    save_figure(figure, out, PNG_NAMES[2])


def render(data, out):
    data, out = Path(data), Path(out)
    axes_index = by_key(read_rows(data / 'round1_axis_segments.csv'),
                        'figure', 'source', 'scenario', 'panel', 'ordinate')
    axes = {key: AxisScale(json.loads(row['segments'])) for key, (row,) in axes_index.items()}
    out.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False}):
        draw_overall(data, axes, out)
        draw_learning(data, axes, out)
        draw_doses(data, axes, out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    command = parser.add_mutually_exclusive_group(required=True)
    command.add_argument('--reports', type=Path, help='prepare plotting CSVs from formal report tables')
    command.add_argument('--output', type=Path, help='render three PNGs from existing plotting CSVs')
    args = parser.parse_args()
    if args.reports is not None:
        prepare(args.reports, args.data)
    else:
        render(args.data, args.output)
