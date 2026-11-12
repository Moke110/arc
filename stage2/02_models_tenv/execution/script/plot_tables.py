"""ADR 0064 plot-value tables for the tenv comparison figures (#168).

The renderers in #169 must reproduce ``metrics.png`` and
``metrics_indirect.png`` from these CSVs alone, so this module emits the
long-form plot-value schema of the shared plotting API (#162): one CSV per
figure with the columns in :data:`PLOT_VALUE_COLUMNS`, one ``record="figure"``
row, one ``record="panel"`` row per panel, one ``record="series"`` row per drawn
series and one ``record="point"`` row per drawn value.

The metric set is the ADR 0064 family. ``metrics_values`` carries the four direct
panels, one row by four columns. ``indirect_metrics_values`` carries the
two-row-by-four-column reliability figure: the four ``@d`` curves over the
scale-free axis ``d'`` (with each run's ``d*`` marker and the ``F1@d`` reference
line on the reliability panel), the ``Reliability@d*`` and ``d*`` summaries, the
macro three-class confusion at ``d*`` and the null/changed Reliability split.
No retired ``@k`` metric name appears here.

This module is presentation only: it never computes a scientific metric. The
shared drawing API adapts to these panel kinds in #169.
"""
import numpy as np
import pandas as pd

METRICS_SCHEMA = 'stage2-tenv-metrics-values-v1'
METRICS_INDIRECT_SCHEMA = 'stage2-tenv-metrics-indirect-values-v1'
LEARNING_CURVE_SCHEMA = 'stage2-tenv-learning-curve-values-v1'
METRICS_TITLE = 'Tenv model comparison'
METRICS_INDIRECT_TITLE = 'Tenv reliability threshold scan'
LEARNING_CURVE_TITLE = 'Learning curve'
PLOT_VALUE_COLUMNS = ('record', 'schema', 'kind', 'panel', 'title', 'x_label', 'y_label',
                      'series', 'hue', 'shade', 'role', 'run', 'x', 'y', 'defined', 'reason')
DIRECT_LAYOUT = (('PCC', 'Spearman', 'RMSE', 'R2'),)
DIRECT_PANEL_ORDER = DIRECT_LAYOUT[0]
INDIRECT_LAYOUT = (('Precision@d', 'Recall@d', 'Spearman@d', 'Reliability@d'),
                   ('Reliability@d*', 'd*', 'Confusion@d*', 'NullSplit@d*'))
INDIRECT_PANEL_ORDER = tuple(metric for row in INDIRECT_LAYOUT for metric in row)
CURVE_METRICS = ('Precision@d', 'Recall@d', 'Spearman@d', 'Reliability@d')
CURVE_INPUT = {'Precision@d': 'precision', 'Recall@d': 'recall',
               'Spearman@d': 'spearman', 'Reliability@d': 'reliability'}
BAR_METRICS = {'Reliability@d*': 'reliability_at_d_star',
               'd*': 'd_star_prime'}
RELIABILITY_PANEL = 'Reliability@d'
CONFUSION_PANEL = 'Confusion@d*'
SPLIT_PANEL = 'NullSplit@d*'
CONFUSION_CLASSES = ('up', 'not_significant', 'down')
# Palette constants from the shared palette (#162); used only for hue/shade
# identity, never to draw here. Entry order: trained entrants first, gray last.
ENTRY_HUES = ('blue', 'orange', 'green', 'purple', 'magenta', 'gray')
LINE_SHADE = 2
RUN_SHADE_ORDER = (2, 4, 6, 8, 1, 3, 5, 7, 0)


def shared_plot_api():
    """The shared plotting API, when this checkout carries it."""
    try:
        from bases import plots
    except ImportError:
        return None
    if not hasattr(plots, 'PLOT_VALUE_COLUMNS') or not hasattr(plots, 'write_plot_values'):
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


def _d_prime_axis(entries):
    """Every run of every entry shares the scale-free ``d' = k/100`` axis."""
    axis = None
    for entry in entries:
        for run in entry['runs']:
            values = tuple(float(item) for item in run['d_prime'])
            if axis is None:
                axis = values
            elif values != axis:
                raise ValueError('every run must share one d_prime axis')
    if axis is None:
        raise ValueError('at least one plotted entry is required')
    return axis


def _bar_panel(entries, metric, *, x_label='Entry'):
    panel_entries = []
    for entry in entries:
        values = tuple(None if run['direct'].get(metric) is None
                       else float(run['direct'][metric]) for run in entry['runs'])
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
    return {'panel': metric, 'kind': 'bar', 'title': metric,
            'x_label': x_label, 'y_label': metric, 'entries': panel_entries}


def _curve_panel(entries, axis, metric, field):
    panel_entries = []
    for entry in entries:
        runs = tuple(tuple(float(item) for item in run[field]) for run in entry['runs'])
        matrix = np.vstack([np.asarray(run, dtype=np.float64) for run in runs])
        mean, minimum, maximum = matrix.mean(axis=0), matrix.min(axis=0), matrix.max(axis=0)
        entry_values = {
            'series': entry['label'], 'hue': entry['hue'],
            'run_labels': tuple(run['run'] for run in entry['runs']),
            'runs': runs, 'mean': tuple(float(item) for item in mean),
            'minimum': tuple(float(item) for item in minimum),
            'maximum': tuple(float(item) for item in maximum)}
        if field == 'reliability':
            entry_values['d_star'] = tuple(
                {'run': run['run'], 'x': float(run['d_star_prime']),
                 'y': float(run['reliability_at_d_star'])} for run in entry['runs'])
            reference = np.vstack([np.asarray(run['f1'], dtype=np.float64)
                                   for run in entry['runs']]).mean(axis=0)
            entry_values['reference'] = tuple(float(item) for item in reference)
        panel_entries.append(entry_values)
    return {'panel': metric, 'kind': 'curve', 'title': metric,
            'x_label': "d'", 'y_label': metric, 'd_prime': axis,
            'entries': panel_entries}


def _summary_bar_panel(entries, metric):
    field = BAR_METRICS[metric]
    panel_entries = []
    for entry in entries:
        values = tuple(float(run[field]) for run in entry['runs'])
        panel_entries.append({
            'series': entry['label'], 'hue': entry['hue'],
            'run_labels': tuple(run['run'] for run in entry['runs']),
            'values': values, 'defined': tuple(True for _ in values),
            'mean': float(np.mean(values)), 'mean_defined': True,
            'reason': None, 'undefined_runs': ()})
    return {'panel': metric, 'kind': 'bar', 'title': metric,
            'x_label': 'Entry', 'y_label': metric, 'entries': panel_entries}


def _confusion_panel(entries):
    """Sum each entry's per-run 3x3 counts; the macro diagnostic at ``d*``."""
    panel_entries = []
    for entry in entries:
        total = np.zeros((len(CONFUSION_CLASSES), len(CONFUSION_CLASSES)), dtype=np.int64)
        for run in entry['runs']:
            total += np.asarray(run['confusion']['matrix'], dtype=np.int64)
        panel_entries.append({'series': entry['label'], 'hue': entry['hue'],
                              'classes': CONFUSION_CLASSES,
                              'matrix': tuple(tuple(int(value) for value in row)
                                              for row in total)})
    return {'panel': CONFUSION_PANEL, 'kind': 'heatmap', 'title': CONFUSION_PANEL,
            'x_label': 'measured class', 'y_label': 'predicted class',
            'classes': CONFUSION_CLASSES, 'entries': panel_entries}


def _split_panel(entries):
    """The null/changed macro Reliability at ``d*``, one series per class."""
    panel_entries = []
    for entry in entries:
        for side, field in (('null', 'null_reliability_at_d_star'),
                            ('changed', 'changed_reliability_at_d_star')):
            values = tuple(float(run[field]) for run in entry['runs'])
            panel_entries.append({
                'series': f"{entry['label']} | {side}", 'hue': entry['hue'],
                'run_labels': tuple(run['run'] for run in entry['runs']),
                'values': values, 'defined': tuple(True for _ in values),
                'mean': float(np.mean(values)), 'mean_defined': True,
                'reason': None, 'undefined_runs': ()})
    return {'panel': SPLIT_PANEL, 'kind': 'bar', 'title': SPLIT_PANEL,
            'x_label': 'Group class', 'y_label': 'Reliability@d*', 'entries': panel_entries}


def _metrics_values(entries, title=METRICS_TITLE):
    """The four direct metric panels, one row by four columns."""
    panels = [_bar_panel(entries, metric) for metric in DIRECT_PANEL_ORDER]
    return {'schema': METRICS_SCHEMA, 'kind': 'metrics', 'title': title,
            'layout': {'rows': 1, 'columns': 4, 'panels': DIRECT_PANEL_ORDER},
            'panels': panels}


def metrics_values(entries, *, title=METRICS_TITLE):
    """Complete direct-metric plotted values (ADR 0064, one row by four columns)."""
    return _metrics_values(entries, title)


def _indirect_values(entries, title=METRICS_INDIRECT_TITLE):
    axis = _d_prime_axis(entries)
    panels = []
    for metric in INDIRECT_PANEL_ORDER:
        if metric in CURVE_METRICS:
            panels.append(_curve_panel(entries, axis, metric, CURVE_INPUT[metric]))
        elif metric in BAR_METRICS:
            panels.append(_summary_bar_panel(entries, metric))
        elif metric == CONFUSION_PANEL:
            panels.append(_confusion_panel(entries))
        elif metric == SPLIT_PANEL:
            panels.append(_split_panel(entries))
        else:  # pragma: no cover - the layout is closed
            raise ValueError(f'unknown indirect panel: {metric!r}')
    return {'schema': METRICS_INDIRECT_SCHEMA, 'kind': 'metrics_indirect',
            'title': title,
            'layout': {'rows': 2, 'columns': 4, 'panels': INDIRECT_PANEL_ORDER},
            'panels': panels}


def indirect_metrics_values(entries, *, title=METRICS_INDIRECT_TITLE):
    """Complete ``@d`` reliability plotted values (two rows by four columns)."""
    return _indirect_values(entries, title)


def _learning_values(models, title=LEARNING_CURVE_TITLE):
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
            else _learning_values(models, title))


def _blank(fields):
    row = dict.fromkeys(PLOT_VALUE_COLUMNS, '')
    row.update({'x': float('nan'), 'y': float('nan'), 'defined': 1, 'shade': '',
                'reason': ''})
    row.update(fields)
    return row


def _plot_rows_metrics(values):
    schema, kind = values['schema'], values['kind']
    rows = [_blank({'record': 'figure', 'schema': schema, 'kind': kind, 'title': values['title']})]
    for panel in values['panels']:
        rows.append(_blank({'record': 'panel', 'schema': schema, 'kind': panel['kind'],
                            'panel': panel['panel'], 'title': panel['title'],
                            'x_label': panel['x_label'], 'y_label': panel['y_label']}))
        for entry in panel['entries']:
            rows.append(_blank({'record': 'series', 'schema': schema, 'kind': kind,
                                'panel': panel['panel'], 'series': entry['series'],
                                'hue': entry['hue'], 'shade': LINE_SHADE}))
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
    return rows


def _plot_rows_indirect(values):
    schema, kind = values['schema'], values['kind']
    rows = [_blank({'record': 'figure', 'schema': schema, 'kind': kind, 'title': values['title']})]
    for panel in values['panels']:
        rows.append(_blank({'record': 'panel', 'schema': schema, 'kind': panel['kind'],
                            'panel': panel['panel'], 'title': panel['title'],
                            'x_label': panel['x_label'], 'y_label': panel['y_label']}))
        if panel['kind'] == 'heatmap':
            for entry in panel['entries']:
                rows.append(_blank({'record': 'series', 'schema': schema, 'kind': kind,
                                    'panel': panel['panel'], 'series': entry['series'],
                                    'hue': entry['hue'], 'shade': LINE_SHADE}))
                for row_index, predicted in enumerate(entry['classes']):
                    for column_index, measured in enumerate(entry['classes']):
                        rows.append(_blank({
                            'record': 'point', 'schema': schema, 'kind': kind,
                            'panel': panel['panel'], 'series': entry['series'],
                            'role': 'heat_cell', 'x': float(column_index),
                            'y': int(entry['matrix'][row_index][column_index]),
                            'reason': f'{predicted}|{measured}'}))
            continue
        for entry in panel['entries']:
            rows.append(_blank({'record': 'series', 'schema': schema, 'kind': kind,
                                'panel': panel['panel'], 'series': entry['series'],
                                'hue': entry['hue'], 'shade': LINE_SHADE}))
        if panel['kind'] == 'curve':
            for entry in panel['entries']:
                for index, x in enumerate(panel['d_prime']):
                    for role, column in (('curve_mean', 'mean'), ('curve_min', 'minimum'),
                                         ('curve_max', 'maximum')):
                        rows.append(_blank({'record': 'point', 'schema': schema, 'kind': kind,
                                            'panel': panel['panel'], 'series': entry['series'],
                                            'role': role, 'x': float(x),
                                            'y': entry[column][index], 'defined': 1}))
                for run, curve in zip(entry['run_labels'], entry['runs']):
                    for x, y in zip(panel['d_prime'], curve):
                        rows.append(_blank({'record': 'point', 'schema': schema, 'kind': kind,
                                            'panel': panel['panel'], 'series': entry['series'],
                                            'role': 'curve_run', 'run': run, 'x': float(x),
                                            'y': y, 'defined': 1}))
                if 'reference' in entry:
                    reference = f"{entry['series']} | F1@d"
                    rows.append(_blank({'record': 'series', 'schema': schema, 'kind': kind,
                                        'panel': panel['panel'], 'series': reference,
                                        'hue': entry['hue'], 'shade': LINE_SHADE}))
                    for x, y in zip(panel['d_prime'], entry['reference']):
                        rows.append(_blank({'record': 'point', 'schema': schema, 'kind': kind,
                                            'panel': panel['panel'], 'series': reference,
                                            'role': 'reference', 'x': float(x), 'y': y,
                                            'defined': 1}))
                for marker in entry.get('d_star', ()):
                    rows.append(_blank({'record': 'point', 'schema': schema, 'kind': kind,
                                        'panel': panel['panel'], 'series': entry['series'],
                                        'role': 'marker', 'run': marker['run'],
                                        'x': marker['x'], 'y': marker['y'], 'defined': 1}))
            continue
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
    return rows


def plot_rows(values):
    """Flatten the plotted values into the recorded long-form CSV rows."""
    kind = values['kind']
    if kind == 'metrics':
        return _plot_rows_metrics(values)
    if kind == 'metrics_indirect':
        return _plot_rows_indirect(values)
    if kind == 'learning_curve':
        return _learning_rows(values)
    raise ValueError(f'unknown plotted-values kind: {kind!r}')


def _learning_rows(values):
    schema, kind = values['schema'], values['kind']
    rows = [_blank({'record': 'figure', 'schema': schema, 'kind': kind, 'title': values['title']})]
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
    return rows


def write_plot_values(values, path):
    """Write the complete plotted values as the recorded long-form CSV.

    The serializer is local because the ``@d`` figure panels are new in this
    ticket; the shared drawing API adapts to them in #169. The output schema is
    the shared :data:`PLOT_VALUE_COLUMNS` schema byte for byte.
    """
    path = str(path)
    frame = pd.DataFrame(plot_rows(values), columns=list(PLOT_VALUE_COLUMNS))
    frame.to_csv(path, index=False)
    return path
