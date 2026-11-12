"""ADR 0063 plot-value tables for the tenv comparison figures (#164).

The renderers in #165 must reproduce ``metrics.png`` and ``learning_curve.png``
from these CSVs alone, so this module emits the long-form plot-value schema of
the shared plotting API (#162): one CSV per figure with the columns in
:data:`PLOT_VALUE_COLUMNS`, one ``record="figure"`` row, one ``record="panel"``
row per panel, one ``record="series"`` row per drawn series and one
``record="point"`` row per drawn value (curve mean/min/max, one curve per run,
bar mean and one bar point per run, or learning-curve line points).

The shared API is used when it is present (after #162 merges); until then this
module carries a byte-compatible fallback so the driver still emits the exact
schema. The fallback is presentation only: it never computes a scientific
metric.
"""
import numpy as np
import pandas as pd

METRICS_SCHEMA = 'stage2-tenv-metrics-values-v1'
LEARNING_CURVE_SCHEMA = 'stage2-tenv-learning-curve-values-v1'
METRICS_TITLE = 'Tenv model comparison'
LEARNING_CURVE_TITLE = 'Learning curve'
PLOT_VALUE_COLUMNS = ('record', 'schema', 'kind', 'panel', 'title', 'x_label', 'y_label',
                      'series', 'hue', 'shade', 'role', 'run', 'x', 'y', 'defined', 'reason')
METRICS_LAYOUT = (('PCC', 'Spearman', 'RMSE', 'R2'),
                  ('Precision@k', 'Spearman@k', 'max_k@P80S40', 'Recall@max_k'))
CURVE_METRICS = ('Precision@k', 'Spearman@k')
BAR_METRICS = ('PCC', 'Spearman', 'RMSE', 'R2', 'max_k@P80S40', 'Recall@max_k')
PANEL_ORDER = tuple(metric for row in METRICS_LAYOUT for metric in row)
METRIC_INPUT = {'PCC': ('direct', 'PCC'), 'Spearman': ('direct', 'Spearman'),
                'RMSE': ('direct', 'RMSE'), 'R2': ('direct', 'R2'),
                'max_k@P80S40': ('indirect', 'max_k'),
                'Recall@max_k': ('indirect', 'recall')}
CURVE_INPUT = {'Precision@k': 'precision', 'Spearman@k': 'spearman'}
# Palette constants from the shared palette (#162); used only for hue/shade
# identity, never to draw here. Entry order: trained entrants first, gray last.
ENTRY_HUES = ('blue', 'orange', 'green', 'purple', 'magenta', 'gray')
LINE_SHADE = 2
RUN_SHADE_ORDER = (2, 4, 6, 8, 1, 3, 5, 7, 0)


def shared_plot_api():
    """The shared ADR 0063 plotting API, when #162 is merged into this checkout."""
    try:
        from bases import plots
    except ImportError:
        return None
    if not hasattr(plots, 'METRICS_SCHEMA') or not hasattr(plots, 'write_plot_values'):
        return None
    return plots


def palette_entry_hue(index):
    try:
        from bases import palette
        return palette.entry_hue(index)
    except ImportError:
        return ENTRY_HUES[index % len(ENTRY_HUES)]


def _entry_reason(entry, metric):
    return next((run['reasons'].get(metric) for run in entry['runs']
                 if run['reasons'].get(metric)), None)


def _metrics_values_fallback(entries, title=METRICS_TITLE):
    axis = None
    for entry in entries:
        for run in entry['runs']:
            values = tuple(float(item) for item in run['indirect']['k'])
            if axis is None:
                axis = values
            elif values != axis:
                raise ValueError('every run must share one indirect k axis')
    panels = []
    for metric in PANEL_ORDER:
        if metric in CURVE_METRICS:
            field = CURVE_INPUT[metric]
            panel_entries = []
            for entry in entries:
                runs = tuple(tuple(float(item) for item in run['indirect'][field])
                             for run in entry['runs'])
                matrix = np.vstack([np.asarray(run, dtype=np.float64) for run in runs])
                defined = np.isfinite(matrix)
                count = defined.sum(axis=0)
                mean = np.divide(np.where(defined, matrix, 0.0).sum(axis=0), count,
                                 out=np.full(len(axis), np.nan), where=count > 0)
                minimum = np.where(count > 0, np.min(np.where(defined, matrix, np.inf), axis=0), np.nan)
                maximum = np.where(count > 0, np.max(np.where(defined, matrix, -np.inf), axis=0), np.nan)
                mean_defined = tuple(bool(item) for item in count > 0)
                reason = _entry_reason(entry, metric)
                if reason is None and not all(mean_defined):
                    reason = 'undefined'
                panel_entries.append({
                    'series': entry['label'], 'hue': entry['hue'],
                    'run_labels': tuple(run['run'] for run in entry['runs']),
                    'runs': runs, 'mean': tuple(float(item) for item in mean),
                    'minimum': tuple(float(item) for item in minimum),
                    'maximum': tuple(float(item) for item in maximum),
                    'defined_count': tuple(int(item) for item in count),
                    'mean_defined': mean_defined, 'reason': reason})
            panels.append({'panel': metric, 'kind': 'curve', 'title': metric,
                           'x_label': 'k', 'y_label': metric, 'k': axis,
                           'entries': panel_entries})
            continue
        section, key = METRIC_INPUT[metric]
        panel_entries = []
        for entry in entries:
            values = tuple(None if run[section].get(key) is None
                           else float(run[section][key]) for run in entry['runs'])
            defined = tuple(value is not None and np.isfinite(value) for value in values)
            samples = [value for value in values if value is not None and np.isfinite(value)]
            reason = None
            if samples:
                mean = float(np.mean(samples))
            else:
                mean = float('nan')
                reason = _entry_reason(entry, metric) or 'undefined'
            panel_entries.append({
                'series': entry['label'], 'hue': entry['hue'],
                'run_labels': tuple(run['run'] for run in entry['runs']),
                'values': values, 'defined': defined, 'mean': mean,
                'mean_defined': bool(samples), 'reason': reason,
                'undefined_runs': tuple(index for index, item in enumerate(defined) if not item)})
        panels.append({'panel': metric, 'kind': 'bar', 'title': metric,
                       'x_label': 'Entry', 'y_label': metric, 'k': None,
                       'entries': panel_entries})
    return {'schema': METRICS_SCHEMA, 'kind': 'metrics', 'title': title,
            'layout': {'rows': 2, 'columns': 4, 'panels': PANEL_ORDER}, 'panels': panels}


def metrics_values(entries, *, title=METRICS_TITLE):
    """Complete eight-panel plotted values; shared API first, fallback second."""
    shared = shared_plot_api()
    return (shared.metrics_values(entries, title=title) if shared
            else _metrics_values_fallback(entries, title))


def _learning_values_fallback(models, title=LEARNING_CURVE_TITLE):
    panels = []
    for index, model in enumerate(models):
        runs = []
        for position, run in enumerate(model['runs']):
            epoch = tuple(float(item) for item in (run.get('epoch') or ()))
            train = tuple(float(item) for item in (run.get('train') or ()))
            validation = tuple(float(item) for item in (run.get('validation') or ()))
            defined = bool(epoch)
            runs.append({'series': run.get('run') or f'run_{position}',
                         'shade': RUN_SHADE_ORDER[position % len(RUN_SHADE_ORDER)],
                         'train_hue': 'blue', 'validation_hue': 'green',
                         'defined': defined,
                         'reason': None if defined else (run.get('reason') or 'undefined'),
                         'epoch': epoch, 'train': train, 'validation': validation})
        panels.append({'panel': model['name'], 'title': model['name'],
                       'x_label': 'Epoch', 'y_label': 'Loss', 'runs': runs})
    return {'schema': LEARNING_CURVE_SCHEMA, 'kind': 'learning_curve', 'title': title,
            'layout': {'rows': 1, 'columns': len(panels),
                       'panels': tuple(panel['panel'] for panel in panels)},
            'panels': panels}


def learning_curve_values(models, *, title=LEARNING_CURVE_TITLE):
    """Complete learning-curve plotted values; shared API first, fallback second."""
    shared = shared_plot_api()
    return (shared.learning_curve_values(models, title=title) if shared
            else _learning_values_fallback(models, title))


def _blank(fields):
    row = dict.fromkeys(PLOT_VALUE_COLUMNS, '')
    row.update({'x': float('nan'), 'y': float('nan'), 'defined': 1, 'shade': '',
                'reason': ''})
    row.update(fields)
    return row


def plot_rows(values):
    """Flatten the plotted values into the recorded long-form CSV rows."""
    schema, kind = values['schema'], values['kind']
    rows = [_blank({'record': 'figure', 'schema': schema, 'kind': kind, 'title': values['title']})]
    if kind == 'metrics':
        for panel in values['panels']:
            rows.append(_blank({'record': 'panel', 'schema': schema, 'kind': panel['kind'],
                                'panel': panel['panel'], 'title': panel['title'],
                                'x_label': panel['x_label'], 'y_label': panel['y_label']}))
            for entry in panel['entries']:
                rows.append(_blank({'record': 'series', 'schema': schema, 'kind': kind,
                                    'panel': panel['panel'], 'series': entry['series'],
                                    'hue': entry['hue'], 'shade': LINE_SHADE}))
            if panel['kind'] == 'bar':
                for position, entry in enumerate(panel['entries']):
                    rows.append(_blank({'record': 'point', 'schema': schema, 'kind': kind,
                                        'panel': panel['panel'], 'series': entry['series'],
                                        'role': 'bar_mean', 'x': float(position),
                                        'y': entry['mean'], 'defined': int(entry['mean_defined']),
                                        'reason': entry['reason'] or ''}))
                    for run, value, defined in zip(entry['run_labels'], entry['values'],
                                                   entry['defined']):
                        rows.append(_blank({'record': 'point', 'schema': schema, 'kind': kind,
                                            'panel': panel['panel'], 'series': entry['series'],
                                            'role': 'bar_run', 'run': run,
                                            'x': float(position), 'y': value,
                                            'defined': int(defined)}))
                continue
            for entry in panel['entries']:
                for index, k in enumerate(panel['k']):
                    undefined_reason = entry['reason'] if not entry['mean_defined'][index] else ''
                    for role, column in (('curve_mean', 'mean'), ('curve_min', 'minimum'),
                                         ('curve_max', 'maximum')):
                        rows.append(_blank({'record': 'point', 'schema': schema, 'kind': kind,
                                            'panel': panel['panel'], 'series': entry['series'],
                                            'role': role, 'x': float(k),
                                            'y': entry[column][index],
                                            'defined': int(entry['mean_defined'][index]),
                                            'reason': undefined_reason}))
                for run, curve in zip(entry['run_labels'], entry['runs']):
                    for k, value in zip(panel['k'], curve):
                        rows.append(_blank({'record': 'point', 'schema': schema, 'kind': kind,
                                            'panel': panel['panel'], 'series': entry['series'],
                                            'role': 'curve_run', 'run': run, 'x': float(k),
                                            'y': value,
                                            'defined': int(bool(np.isfinite(value)))}))
    elif kind == 'learning_curve':
        for panel in values['panels']:
            rows.append(_blank({'record': 'panel', 'schema': schema, 'kind': kind,
                                'panel': panel['panel'], 'title': panel['title'],
                                'x_label': panel['x_label'], 'y_label': panel['y_label']}))
            for run in panel['runs']:
                for split, hue, trajectory in (('train', run['train_hue'], run['train']),
                                               ('validation', run['validation_hue'],
                                                run['validation'])):
                    series = f"{run['series']}|{split}"
                    rows.append(_blank({'record': 'series', 'schema': schema, 'kind': kind,
                                        'panel': panel['panel'], 'series': series, 'hue': hue,
                                        'shade': run['shade'], 'run': run['series'],
                                        'defined': int(run['defined']),
                                        'reason': run['reason'] or ''}))
                    for x, y in zip(run['epoch'], trajectory):
                        rows.append(_blank({'record': 'point', 'schema': schema, 'kind': kind,
                                            'panel': panel['panel'], 'series': series,
                                            'role': 'line', 'run': run['series'],
                                            'x': float(x), 'y': y,
                                            'defined': int(run['defined'])}))
    else:
        raise ValueError(f'unknown plotted-values kind: {kind!r}')
    return rows


def write_plot_values(values, path):
    """Write the complete plotted values as the recorded long-form CSV."""
    shared = shared_plot_api()
    if shared:
        return shared.write_plot_values(values, path)
    path = str(path)
    frame = pd.DataFrame(plot_rows(values), columns=list(PLOT_VALUE_COLUMNS))
    frame.to_csv(path, index=False)
    return path
