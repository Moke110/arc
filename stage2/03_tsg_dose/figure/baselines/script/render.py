"""Render the Stage 03 successor report with the sealed Round 1 plot primitives."""
import argparse
import json
import math
from pathlib import Path
from types import ModuleType

ECFP4_MLP = 'ECFP4_MLP'
MEAN_CELL = 'Mean_cell'
METRIC_TABLES = ('round1_runs.csv', 'round1_summary.csv', 'round1_dose.csv')
CURVE_TABLE = 'round1_curves.csv'
AXIS_TABLE = 'round1_axis_segments.csv'


def load_reference(path):
    module = ModuleType('round1_renderer')
    exec(compile(path.read_bytes(), str(path), 'exec'), module.__dict__)
    return module


def validate_tables(base, tables):
    """Require exact fold, metric, dose and trained epoch coordinates."""
    original = base.METHODS
    # The predecessor validates the original methods and their exact coordinate coverage.
    base.validate_coordinates({name: [row for row in rows if row['method'] in original]
                               for name, rows in tables.items()})
    runs = {(source, scenario, method, str(fold), '42')
            for source in base.SOURCES for scenario in base.SCENARIOS
            for method in (MEAN_CELL, ECFP4_MLP) for fold in range(3)}
    spreads = {(source, scenario, method, metric, '42')
               for source in base.SOURCES for scenario in base.SCENARIOS
               for method in (MEAN_CELL, ECFP4_MLP) for metric in base.METRICS}
    doses = {(source, 'cold-drug', method, str(dose), '42')
             for source, values in zip(base.SOURCES, ((.05, .5, 5.), (.2, 1., 10.)), strict=True)
             for method in (MEAN_CELL, ECFP4_MLP) for dose in values}
    curves = {key for key in runs if key[2] == ECFP4_MLP}
    fields = ('source', 'scenario', 'method')
    for name, field, expected in zip((*METRIC_TABLES, CURVE_TABLE),
                                    ('fold', 'metric', 'dose_um', 'fold'),
                                    (runs, spreads, doses, curves), strict=True):
        rows = tables[name]
        if name == METRIC_TABLES[2]:
            rows = [dict(row, dose_um=str(float(row['dose_um']))) for row in rows]
        groups = base.by_key([row for row in rows if row['method'] not in original],
                             *fields, field, 'seed')
        if set(groups) != expected:
            raise ValueError(f'baseline coordinate inventory mismatch: {name}')
        for coordinates in groups.values():
            if name == CURVE_TABLE:
                epochs = sorted(int(row['epoch']) for row in coordinates)
                if epochs != list(range(1, len(epochs) + 1)):
                    raise ValueError('incomplete ECFP4_MLP epochs')
            elif len(coordinates) != 1:
                raise ValueError(f'duplicate baseline coordinate: {name}')


def separated_metric_axis(base, intervals):
    """Break only an empty range wide enough to obscure separated groups."""
    coverage = sorted([*intervals, (0., 0.)])
    merged = []
    for lower, upper in coverage:
        if merged and lower <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(upper, merged[-1][1]))
        else:
            merged.append((lower, upper))
    if len(merged) < 2:
        return None
    lower, upper = merged[0][0], merged[-1][1]
    left, right = max(((first[1], second[0]) for first, second in zip(merged, merged[1:])),
                      key=lambda gap: gap[1] - gap[0])
    left_width, right_width = left - lower, upper - right
    if right - left <= 1.25 * max(left_width, right_width):
        return None
    gap = right - left
    low_pad = min(gap * .09, max(left_width * .09, gap * .015))
    high_pad = min(gap * .09, max(right_width * .09, gap * .015))
    left_count = sum(stop <= left for _, stop in intervals)
    right_count = sum(start >= right for start, _ in intervals)
    if not left_count or not right_count:
        return None
    split = .28 if left_count < right_count else .62
    return base.AxisScale([
        base.linear_segment(lower - low_pad, left + low_pad, 0., split,
                            show_zero=lower <= 0 <= left),
        base.linear_segment(right - high_pad, upper + high_pad, split + .10, 1.,
                            show_zero=right <= 0 <= upper),
    ])


def metric_axis_records(base, tables):
    axes = base.export_axes(tables)
    folds = base.by_key(tables[METRIC_TABLES[0]], 'source', 'scenario')
    spreads = base.by_key(tables[METRIC_TABLES[1]], 'source', 'scenario', 'metric')
    doses = base.by_key(tables[METRIC_TABLES[2]], 'source', 'method')
    for record in axes:
        source, scenario, metric = (record[key] for key in ('source', 'scenario', 'panel'))
        intervals = []
        if record['figure'] == 'overall':
            for row in folds[source, scenario]:
                point = base.number(row[metric])
                if math.isfinite(point):
                    intervals.append((point, point))
            for row in spreads[source, scenario, metric]:
                mean, sd = base.number(row['mean']), base.number(row['sd'])
                if math.isfinite(mean):
                    radius = sd if math.isfinite(sd) else 0.
                    intervals.append((mean - radius, mean + radius))
        elif record['figure'] == 'dose':
            for (sample_source, _), rows in doses.items():
                if sample_source == source:
                    points = [base.number(row[metric]) for row in rows]
                    finite = [value for value in points if math.isfinite(value)]
                    if finite:
                        intervals.append((min(finite), max(finite)))
        axis = separated_metric_axis(base, intervals)
        if axis is not None:
            record['segments'] = json.dumps(axis.segments, separators=(',', ':'))
    return axes


def draw_learning(base, data, axes, output):
    curves = base.by_key(base.read_rows(data / CURVE_TABLE), 'source', 'scenario', 'method', 'fold')
    grid = tuple((source, scenario) for source in base.SOURCES for scenario in base.SCENARIOS)
    figure, panels = base.plt.subplots(len(grid), len(base.METHODS) + 1, figsize=(46, 17))
    for (source, scenario), row in zip(grid, panels, strict=True):
        base.source_label(row[0], source, scenario)
        overview = row[-1]
        log_scale = axes['learning', source, scenario, 'all-methods', 'mse']
        log_scale.apply(overview, color='#2e4150')
        for ax, method, label, color in zip(row[:-1], base.METHODS, base.LABELS, base.COLORS, strict=True):
            right = ax.twinx()
            mse = axes['learning', source, scenario, method, 'mse']
            loss = axes['learning', source, scenario, method, 'loss']
            loss.apply(ax, color=base.TRAIN_COLORS[0])
            right.spines['right'].set_visible(True)
            mse.apply(right, side='right', zero_line=False, color=base.VALID_COLORS[0])
            right.axhline(mse.map([0])[0], color=base.VALID_COLORS[0], lw=.55, alpha=.55)
            for fold in range(3):
                rows = sorted(curves[source, scenario, method, str(fold)],
                              key=lambda record: int(record['epoch']))
                epochs = [int(record['epoch']) for record in rows]
                train = [base.number(record['train_loss']) for record in rows]
                validation = [base.number(record['validation_mse']) for record in rows]
                loss.plot(ax, epochs, train, '--', color=base.TRAIN_COLORS[fold], lw=1)
                mse.plot(right, epochs, validation, color=base.VALID_COLORS[fold], lw=1)
                log_scale.plot(overview, epochs, validation, color=color, alpha=.55, lw=.95,
                               label=label if fold == 0 else None)
            ax.set(title=label, xlabel='Epoch')
            ax.set_ylabel('Training loss', color=base.TRAIN_COLORS[0], fontsize=9)
            right.set_ylabel('Validation MSE', color=base.VALID_COLORS[0], fontsize=9)
        overview.set(title='All trained methods', xlabel='Epoch',
                     ylabel='Validation MSE · log data segment')
        overview.legend(loc='upper right', fontsize=7, framealpha=.86)
    handles = [base.Line2D([0], [0], color=color, ls=style, label=f'Fold {fold} · {name}')
               for name, colors, style in [('train', base.TRAIN_COLORS, '--'),
                                          ('validation', base.VALID_COLORS, '-')]
               for fold, color in enumerate(colors)]
    figure.legend(handles=handles, loc='lower center', ncol=6, bbox_to_anchor=(.5, .018))
    figure.suptitle('Round 1 · training loss and validation MSE by epoch · seed 42', fontsize=21)
    figure.subplots_adjust(left=.038, right=.985, top=.92, bottom=.087, hspace=.78, wspace=.72)
    base.save_figure(figure, output, base.PNG_NAMES[1])


def render(reference, data, output, *, prepare_axes=False):
    base = load_reference(reference)
    tables = {name: base.read_rows(data / name) for name in (*METRIC_TABLES, CURVE_TABLE)}
    validate_tables(base, tables)
    # This instance belongs to this invocation. The sealed renderer on disk is unchanged.
    base.METHODS += (ECFP4_MLP,)
    base.LABELS += (ECFP4_MLP,)
    base.COLORS += ('#087c68',)
    if prepare_axes:
        base.write_rows(data / AXIS_TABLE, metric_axis_records(base, tables))
    index = base.by_key(base.read_rows(data / AXIS_TABLE), 'figure', 'source', 'scenario', 'panel', 'ordinate')
    axes = {key: base.AxisScale(json.loads(record['segments'])) for key, (record,) in index.items()}
    output.mkdir(parents=True, exist_ok=True)
    with base.plt.rc_context({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False}):
        draw_learning(base, data, axes, output)
        base.METHODS += (MEAN_CELL,)
        base.LABELS += (MEAN_CELL,)
        base.COLORS += ('#292929',)
        base.draw_overall(data, axes, output)
        base.draw_doses(data, axes, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path,
                        default=Path(__file__).resolve().parents[2] / 'script/render_axis_local.py')
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    render(args.reference, args.data, args.output)
