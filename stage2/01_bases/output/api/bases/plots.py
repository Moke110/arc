"""Stage 2 plotting APIs: the tenv retrieval panel and the ADR 0064 comparison figures.

Every Stage 2 comparison drawing call takes only recorded values and returns the
complete plotted values it drew, so a step renderer can rebuild a PNG from its
``figure/data/*.csv`` alone (ADR 0062, ADR 0064).  Colours come only from the
stage palette in :mod:`bases.palette`; undefined values are NaN with a reason and
are drawn explicitly, never as zero.

The comparison figures are two recorded products:

* ``metrics.png`` — one row by four columns of the direct-metric bars
  (``PCC``, ``Spearman``, ``RMSE``, ``R2``), schema :data:`METRICS_SCHEMA`.
* ``metrics_indirect.png`` — two rows by four columns (ADR 0064): the four ``@d``
  curves over the scale-free axis ``d'`` with the four-run mean line and
  min–max band, each run's ``d*`` marker and the ``F1@d`` reference line on the
  Reliability panel; then ``Reliability@d*`` and ``d*`` bars of the four-run
  mean with one dot per run, the macro three-class confusion heat map at ``d*``
  and the null/changed Reliability split, schema :data:`INDIRECT_SCHEMA`.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from . import palette
from .metrics import METRICS, PR_FIELDS


def render_retrieval(runs, output, *, source='tenv', diagnostic=False):
    """Write a three-direction PR panel and replayable CSV/JSON.

    Each run supplies label, equal-group macro summary and fixed-group pooled
    curves from score_groups. Pooled curve geometry is never labelled as its AP.
    No Tahoe significance or PR publication is supported.
    """
    if source != 'tenv':
        raise ValueError('significance retrieval panels are tenv-only')
    panel_title = ('Diagnostic · measured delta vs tenv truth (not model performance)' if diagnostic
             else 'Tenv · fixed-group pooled precision–recall')
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    output = Path(output)
    if output.suffix != '.png' or not runs:
        raise ValueError('a PNG output and at least one named run are required')
    labels = [str(run['label']) for run in runs]
    if len(labels) != len(set(labels)) or any(not x for x in labels):
        raise ValueError('run labels must be unique and nonempty')
    output.parent.mkdir(parents=True,exist_ok=True)
    rows, summaries = [], []
    fig, axes = plt.subplots(1,3,figsize=(13.5,4.8),sharex=True,sharey=True,layout='constrained')
    titles = ('Combined', 'Upregulated', 'Downregulated')
    colors = plt.get_cmap('tab10').colors
    try:
        for index, run in enumerate(runs):
            label = str(run['label'])
            summaries.append({'label':label, **{name:run['summary'].get(name) for name in METRICS}})
            for name, ax in zip(METRICS,axes):
                curve = run['curves'][name].copy()
                if not set(PR_FIELDS) <= set(curve):
                    raise ValueError('pooled curve schema differs')
                value = run['summary'].get(name)
                if value is not None and (not np.isfinite(value) or not 0 <= value <= 1):
                    raise ValueError('invalid macro AP')
                if len(curve):
                    if (not np.isfinite(curve[['recall','precision']]).all().all()
                            or not curve[['recall','precision']].ge(0).all().all()
                            or not curve[['recall','precision']].le(1).all().all()
                            or not curve.recall.is_monotonic_increasing
                            or curve.eligible_groups.nunique() != 1
                            or not curve.contributing_groups.eq(curve.eligible_groups).all()):
                        raise ValueError('expected fixed-population monotonic PR source')
                    ax.step(curve.recall,curve.precision,where='post',color=colors[index%len(colors)],
                            linewidth=1.8,label=f'{label} · macro AP {value:.3f}' if value is not None else label)
                    rows.append(curve.assign(label=label,metric=name))
                elif value is not None:
                    raise ValueError('eligible AP has no PR curve')
        for ax,title in zip(axes,titles):
            ax.set(title=title,xlabel='Recall',xlim=(0,1.01),ylim=(0,1.04))
            ax.grid(alpha=.18)
            if ax.get_legend_handles_labels()[0]:
                ax.legend(loc='upper center',bbox_to_anchor=(.5,-.23),fontsize=8,frameon=False)
            else:
                ax.text(.5,.5,'No eligible measured groups',ha='center',va='center',transform=ax.transAxes)
        axes[0].set_ylabel('Precision')
        fig.suptitle(panel_title,fontsize=13)
        fig.savefig(output,dpi=160)
    finally:
        plt.close(fig)
    frame = pd.concat(rows,ignore_index=True) if rows else pd.DataFrame(columns=[*PR_FIELDS,'label','metric'])
    frame.to_csv(output.with_suffix('.csv'),index=False)
    report = {'schema':'stage2-indexed-pr-panel-v1','source':source,'diagnostic':bool(diagnostic),'title':panel_title,'summaries':summaries,
              'curve_semantics':'fixed eligible groups pooled at shared magnitude thresholds; displayed AP is exact equal-group macro',
              'files':{output.name:hashlib.sha256(output.read_bytes()).hexdigest(),
                       output.with_suffix('.csv').name:hashlib.sha256(output.with_suffix('.csv').read_bytes()).hexdigest()}}
    output.with_suffix('.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    return report


def replay_retrieval(report, output):
    """Replay the panel solely from its sealed CSV and macro summaries."""
    report = Path(report)
    value = json.loads(report.read_text())
    if value.get('schema') != 'stage2-indexed-pr-panel-v1' or value.get('source') != 'tenv':
        raise ValueError('unknown tenv PR panel source')
    source = report.with_suffix('.csv')
    if hashlib.sha256(source.read_bytes()).hexdigest() != value['files'].get(source.name):
        raise ValueError('PR numeric source seal differs')
    frame = pd.read_csv(source)
    runs = []
    for summary in value['summaries']:
        label = summary['label']
        curves = {name:frame.loc[frame.label.eq(label)&frame.metric.eq(name),list(PR_FIELDS)] for name in METRICS}
        runs.append({'label':label,'summary':summary,'curves':curves})
    return render_retrieval(runs,output,source=value['source'],diagnostic=value.get('diagnostic',False))


# ---------------------------------------------------------------------------
# ADR 0064 comparison figures: the four direct-metric panels and the reliability
# threshold scan.  Everything below is presentation of already-computed metric
# values; no scientific quantity is recomputed here.
# ---------------------------------------------------------------------------
METRICS_LAYOUT = (("PCC", "Spearman", "RMSE", "R2"),)
INDIRECT_LAYOUT = (("Precision@d", "Recall@d", "Spearman@d", "Reliability@d"),
                   ("Reliability@d*", "d*", "Confusion@d*", "NullSplit@d*"))
PANEL_ORDER = tuple(metric for row in METRICS_LAYOUT for metric in row)
INDIRECT_PANEL_ORDER = tuple(metric for row in INDIRECT_LAYOUT for metric in row)
CURVE_METRICS = ("Precision@d", "Recall@d", "Spearman@d", "Reliability@d")
CURVE_INPUT = {"Precision@d": "precision", "Recall@d": "recall",
               "Spearman@d": "spearman", "Reliability@d": "reliability"}
SUMMARY_BARS = {"Reliability@d*": "reliability_at_d_star", "d*": "d_star_prime"}
RELIABILITY_PANEL = "Reliability@d"
CONFUSION_PANEL = "Confusion@d*"
SPLIT_PANEL = "NullSplit@d*"
CONFUSION_CLASSES = ("up", "not_significant", "down")
METRICS_TITLE = "Tenv model comparison"
INDIRECT_TITLE = "Tenv reliability threshold scan"
LEARNING_CURVE_TITLE = "Learning curve"
METRICS_SCHEMA = "stage2-tenv-metrics-values-v2"
INDIRECT_SCHEMA = "stage2-tenv-metrics-indirect-values-v1"
LEARNING_CURVE_SCHEMA = "stage2-tenv-learning-curve-values-v1"
# Every recorded plot-value CSV has exactly these columns, whatever the figure.
# A metrics figure records its grid in the figure row's ``x`` (rows) and ``y``
# (columns) columns, so the recorded layout, not a hardcoded shape, drives the
# drawing and replay calls; the panel rows record the panel list in row order.
PLOT_VALUE_COLUMNS = ("record", "schema", "kind", "panel", "title", "x_label", "y_label",
                      "series", "hue", "shade", "role", "run", "x", "y", "defined", "reason")


def _scalar(value, label):
    """Read one run-level number; ``None`` is undefined (NaN), infinity is rejected."""
    if value is None:
        return float("nan")
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{label} must be a number or undefined") from None
    if math.isinf(number):
        raise ValueError(f"{label} must be finite or undefined")
    return number


def _vector(value, length, label):
    """Read one run-level series; NaN entries are undefined, infinity is rejected."""
    if value is None:
        raise ValueError(f"{label} is required")
    try:
        values = tuple(float(item) for item in value)
    except (TypeError, ValueError):
        raise ValueError(f"{label} must be numeric") from None
    if len(values) != length:
        raise ValueError(f"{label} must match the shared d' axis")
    if any(math.isinf(item) for item in values):
        raise ValueError(f"{label} must be finite or undefined")
    return values


def _finite_series(value, label):
    """Read a learning-curve trajectory: finite values only, undefined means empty."""
    try:
        values = tuple(float(item) for item in value)
    except (TypeError, ValueError):
        raise ValueError(f"{label} must be numeric") from None
    if any(not math.isfinite(item) for item in values):
        raise ValueError(f"{label} must be finite")
    return values


def _panel_titles(panels, panel_titles):
    titles = {metric: metric for metric in panels}
    if panel_titles:
        unknown = sorted(set(panel_titles) - set(panels))
        if unknown:
            raise ValueError(f"unknown metric panel title(s): {unknown}")
        titles.update(panel_titles)
    return titles


def _validate_title(title):
    if not isinstance(title, str) or not title.strip():
        raise ValueError("a nonempty figure title is required")
    return title


def _validate_entries(entries):
    """Normalise the comparison entries and reject ambiguous identities."""
    if not isinstance(entries, (list, tuple)) or not entries:
        raise ValueError("at least one comparison entry is required")
    labels, normalised = set(), []
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise ValueError("each comparison entry must be a mapping")
        label = entry.get("label")
        if not isinstance(label, str) or not label.strip():
            raise ValueError("entry labels must be nonempty strings")
        if label in labels:
            raise ValueError("entry labels must be unique")
        labels.add(label)
        hue = entry.get("hue") or palette.entry_hue(index)
        if hue not in palette.FAMILIES:
            raise ValueError(f"entry hue {hue!r} is not a stage palette family")
        runs = entry.get("runs")
        if not isinstance(runs, (list, tuple)) or not runs:
            raise ValueError(f"entry {label!r} needs at least one run sample")
        run_labels, normalised_runs = set(), []
        for position, run in enumerate(runs):
            if not isinstance(run, dict):
                raise ValueError("each run sample must be a mapping")
            run_label = run.get("run") or f"run_{position}"
            if not isinstance(run_label, str) or not run_label.strip():
                raise ValueError("run labels must be nonempty strings")
            if run_label in run_labels:
                raise ValueError("run labels must be unique within an entry")
            run_labels.add(run_label)
            direct = run.get("direct") or {}
            if not isinstance(direct, dict):
                raise ValueError("run direct values must be a mapping")
            normalised_run = dict(run)
            normalised_run.update({"run": run_label, "direct": dict(direct),
                                   "reasons": dict(run.get("reasons") or {})})
            normalised_runs.append(normalised_run)
        normalised.append({"label": label, "hue": hue, "runs": normalised_runs})
    return normalised


def _layout_shape(layout):
    """Return the declared ``(rows, columns)`` of a layout and reject a ragged one."""
    if not layout or not all(layout):
        raise ValueError("a figure layout needs at least one nonempty row")
    columns = len(layout[0])
    if any(len(row) != columns for row in layout):
        raise ValueError("a figure layout must not be ragged")
    return len(layout), columns


def _validate_d_prime(axis):
    """A ``d' = d/M`` axis is strictly increasing with every point in ``(0, 1]``."""
    if not axis:
        raise ValueError("the d_prime (d') axis must not be empty")
    if any(not (item > 0 and item <= 1) for item in axis):
        raise ValueError("every d_prime (d') value must satisfy 0 < d' <= 1")
    if any(later <= earlier for earlier, later in zip(axis, axis[1:])):
        raise ValueError("the d_prime axis (d') must strictly increase")
    return axis


def _d_prime_axis(entries):
    """Return the one shared scale-free ``d' = k/100`` axis every run must carry."""
    axis = None
    for entry in entries:
        for run in entry["runs"]:
            if "d_prime" not in run:
                raise ValueError("every run needs a d_prime axis")
            values = tuple(float(item) for item in run["d_prime"])
            _validate_d_prime(values)
            if axis is None:
                axis = values
            elif values != axis:
                raise ValueError("every run must share one d_prime axis")
    return axis


def _band(matrix):
    """Cross-run mean, minimum, maximum and defined count at every x."""
    defined = np.isfinite(matrix)
    count = defined.sum(axis=0)
    length = matrix.shape[1]
    mean = np.divide(np.where(defined, matrix, 0.0).sum(axis=0), count,
                     out=np.full(length, np.nan), where=count > 0)
    minimum = np.where(count > 0, np.min(np.where(defined, matrix, np.inf), axis=0), np.nan)
    maximum = np.where(count > 0, np.max(np.where(defined, matrix, -np.inf), axis=0), np.nan)
    return mean, minimum, maximum, count


def _run_reason(entry, metric):
    return next((run["reasons"].get(metric) for run in entry["runs"]
                 if run["reasons"].get(metric)), None)


def _direct_bar_panel(entries, metric, title):
    panel_entries = []
    for entry in entries:
        values = tuple(_scalar(run["direct"].get(metric), f"{metric} for {run['run']}")
                       for run in entry["runs"])
        defined = tuple(bool(np.isfinite(value)) for value in values)
        samples = [value for value in values if np.isfinite(value)]
        reason = None
        if samples:
            mean = float(np.mean(samples))
        else:
            mean = float("nan")
            reason = _run_reason(entry, metric) or "undefined"
        panel_entries.append({
            "series": entry["label"], "hue": entry["hue"],
            "run_labels": tuple(run["run"] for run in entry["runs"]),
            "values": values, "defined": defined, "mean": mean,
            "mean_defined": bool(samples), "reason": reason,
            "undefined_runs": tuple(index for index, item in enumerate(defined) if not item)})
    return {"panel": metric, "kind": "bar", "title": title,
            "x_label": "Entry", "y_label": metric, "entries": panel_entries}


def _curve_panel(entries, axis, metric, field, title):
    length = len(axis)
    panel_entries = []
    for entry in entries:
        runs = tuple(_vector(run.get(field), length, f"{metric} for {run['run']}")
                     for run in entry["runs"])
        matrix = np.vstack([np.asarray(run, dtype=np.float64) for run in runs])
        mean, minimum, maximum, count = _band(matrix)
        reason = _run_reason(entry, metric)
        mean_defined = tuple(bool(item) for item in count > 0)
        if reason is None and not all(mean_defined):
            reason = "undefined"
        block = {
            "series": entry["label"], "hue": entry["hue"],
            "run_labels": tuple(run["run"] for run in entry["runs"]),
            "runs": runs, "mean": tuple(float(item) for item in mean),
            "minimum": tuple(float(item) for item in minimum),
            "maximum": tuple(float(item) for item in maximum),
            "defined_count": tuple(int(item) for item in count),
            "mean_defined": mean_defined, "reason": reason}
        if metric == RELIABILITY_PANEL:
            reference = np.vstack([np.asarray(
                _vector(run.get("f1"), length, f"F1@d for {run['run']}"),
                dtype=np.float64) for run in entry["runs"]])
            reference_mean, _, _, _ = _band(reference)
            block["reference"] = tuple(float(item) for item in reference_mean)
            block["d_star"] = tuple(
                {"run": run["run"],
                 "x": _scalar(run.get("d_star_prime"), f"d* for {run['run']}"),
                 "y": _scalar(run.get("reliability_at_d_star"),
                              f"Reliability@d* for {run['run']}")}
                for run in entry["runs"])
        panel_entries.append(block)
    return {"panel": metric, "kind": "curve", "title": title,
            "x_label": "d'", "y_label": metric, "d_prime": axis,
            "entries": panel_entries}


def _summary_bar_panel(entries, metric, field, title, x_label="Entry"):
    panel_entries = []
    for entry in entries:
        values = tuple(_scalar(run.get(field), f"{metric} for {run['run']}")
                       for run in entry["runs"])
        defined = tuple(bool(np.isfinite(value)) for value in values)
        samples = [value for value in values if np.isfinite(value)]
        panel_entries.append({
            "series": entry["label"], "hue": entry["hue"],
            "run_labels": tuple(run["run"] for run in entry["runs"]),
            "values": values, "defined": defined,
            "mean": float(np.mean(samples)) if samples else float("nan"),
            "mean_defined": bool(samples),
            "reason": None if samples else (_run_reason(entry, metric) or "undefined"),
            "undefined_runs": tuple(index for index, item in enumerate(defined) if not item)})
    return {"panel": metric, "kind": "bar", "title": title,
            "x_label": x_label, "y_label": metric, "entries": panel_entries}


def _confusion_panel(entries, title):
    """Sum each entry's per-run 3x3 counts; the macro diagnostic at ``d*``."""
    panel_entries = []
    for entry in entries:
        total = np.zeros((len(CONFUSION_CLASSES), len(CONFUSION_CLASSES)), dtype=np.int64)
        for run in entry["runs"]:
            confusion = run.get("confusion")
            if not isinstance(confusion, dict):
                raise ValueError(f"confusion counts are required for {run['run']}")
            classes = tuple(confusion.get("classes") or ())
            if classes != CONFUSION_CLASSES:
                raise ValueError("confusion classes must be up, not_significant, down")
            matrix = np.asarray(confusion.get("matrix"), dtype=np.float64)
            if matrix.shape != total.shape or not np.isfinite(matrix).all():
                raise ValueError("each confusion matrix must be a finite 3x3 table")
            if (matrix < 0).any() or (matrix != np.floor(matrix)).any():
                raise ValueError("confusion counts must be nonnegative whole numbers")
            total += matrix.astype(np.int64)
        panel_entries.append({"series": entry["label"], "hue": entry["hue"],
                              "classes": CONFUSION_CLASSES,
                              "matrix": tuple(tuple(int(value) for value in row)
                                              for row in total)})
    return {"panel": CONFUSION_PANEL, "kind": "heatmap", "title": title,
            "x_label": "measured class", "y_label": "predicted class",
            "classes": CONFUSION_CLASSES, "entries": panel_entries}


def _split_panel(entries, title):
    """The null/changed macro Reliability at ``d*``, one series per class."""
    panel_entries = []
    for entry in entries:
        for side, field in (("null", "null_reliability_at_d_star"),
                            ("changed", "changed_reliability_at_d_star")):
            values = tuple(_scalar(run.get(field), f"{side} Reliability for {run['run']}")
                           for run in entry["runs"])
            defined = tuple(bool(np.isfinite(value)) for value in values)
            samples = [value for value in values if np.isfinite(value)]
            panel_entries.append({
                "series": f"{entry['label']} | {side}", "hue": entry["hue"],
                "run_labels": tuple(run["run"] for run in entry["runs"]),
                "values": values, "defined": defined,
                "mean": float(np.mean(samples)) if samples else float("nan"),
                "mean_defined": bool(samples), "reason": None,
                "undefined_runs": tuple(index for index, item in enumerate(defined)
                                        if not item)})
    return {"panel": SPLIT_PANEL, "kind": "bar", "title": title,
            "x_label": "Group class", "y_label": "Reliability@d*",
            "entries": panel_entries}


def metrics_values(entries, *, title=METRICS_TITLE, panel_titles=None):
    """Compute the complete four-panel direct-metric plotted values, one row by four.

    The eight-panel @k figure is retired (ADR 0064); each direct panel draws a bar
    of the cross-run mean with one dot per run, and a value defined for no run is
    NaN with a reason.
    """
    title = _validate_title(title)
    titles = _panel_titles(PANEL_ORDER, panel_titles)
    entries = _validate_entries(entries)
    rows, columns = _layout_shape(METRICS_LAYOUT)
    panels = [_direct_bar_panel(entries, metric, titles[metric]) for metric in PANEL_ORDER]
    return {"schema": METRICS_SCHEMA, "kind": "metrics", "title": title,
            "layout": {"rows": rows, "columns": columns, "panels": PANEL_ORDER},
            "panels": panels}


def indirect_metrics_values(entries, *, title=INDIRECT_TITLE, panel_titles=None):
    """Compute the complete eight-panel ``@d`` reliability plotted values.

    Row one carries the four curves over ``d'`` with the cross-run mean line and
    min–max band, the per-run ``d*`` markers and the mean ``F1@d`` reference on
    the Reliability panel.  Row two carries the ``Reliability@d*`` and ``d*``
    bars, the macro three-class confusion heat map and the null/changed split.
    """
    title = _validate_title(title)
    titles = _panel_titles(INDIRECT_PANEL_ORDER, panel_titles)
    entries = _validate_entries(entries)
    axis = _d_prime_axis(entries)
    rows, columns = _layout_shape(INDIRECT_LAYOUT)
    panels = []
    for metric in INDIRECT_PANEL_ORDER:
        if metric in CURVE_METRICS:
            panels.append(_curve_panel(entries, axis, metric, CURVE_INPUT[metric],
                                       titles[metric]))
        elif metric in SUMMARY_BARS:
            panels.append(_summary_bar_panel(entries, metric, SUMMARY_BARS[metric],
                                             titles[metric]))
        elif metric == CONFUSION_PANEL:
            panels.append(_confusion_panel(entries, titles[metric]))
        elif metric == SPLIT_PANEL:
            panels.append(_split_panel(entries, titles[metric]))
        else:  # pragma: no cover - the layout is closed
            raise ValueError(f"unknown indirect panel: {metric!r}")
    return {"schema": INDIRECT_SCHEMA, "kind": "metrics_indirect", "title": title,
            "layout": {"rows": rows, "columns": columns, "panels": INDIRECT_PANEL_ORDER},
            "panels": panels}


def learning_curve_values(models, *, title=LEARNING_CURVE_TITLE):
    """Compute the complete learning-curve plotted values without drawing anything.

    One subplot per learnable model, epoch on the horizontal axis and loss on the
    vertical axis.  Training uses the blue family and validation the green family;
    one shade identifies one ``(fold, seed)`` run in both families.  A run with no
    recorded trajectory is undefined, not a zero-length curve.
    """
    title = _validate_title(title)
    if not isinstance(models, (list, tuple)) or not models:
        raise ValueError("at least one learnable model is required")
    names, panels = set(), []
    for model in models:
        if not isinstance(model, dict):
            raise ValueError("each learnable model must be a mapping")
        name = model.get("name")
        if not isinstance(name, str) or not name.strip():
            raise ValueError("model names must be nonempty strings")
        if name in names:
            raise ValueError("model names must be unique")
        names.add(name)
        runs = model.get("runs")
        if not isinstance(runs, (list, tuple)) or not runs:
            raise ValueError(f"model {name!r} needs at least one run sample")
        if len(runs) > palette.SHADE_COUNT:
            raise ValueError("at most one run per palette shade is supported")
        labels, panel_runs = set(), []
        for index, run in enumerate(runs):
            if not isinstance(run, dict):
                raise ValueError("each run sample must be a mapping")
            label = run.get("run") or f"run_{index}"
            if not isinstance(label, str) or not label.strip():
                raise ValueError("run labels must be nonempty strings")
            if label in labels:
                raise ValueError("run labels must be unique within a model")
            labels.add(label)
            epoch = _finite_series(run.get("epoch") or (), f"{label} epoch")
            train = _finite_series(run.get("train") or (), f"{label} training loss")
            validation = _finite_series(run.get("validation") or (), f"{label} validation loss")
            lengths = (len(epoch), len(train), len(validation))
            if 0 < sum(count == 0 for count in lengths) < len(lengths):
                raise ValueError(f"{label} must provide epoch, train and validation together")
            if len(train) != len(epoch) or len(validation) != len(epoch):
                raise ValueError(f"{label} trajectories must share one epoch axis")
            if any(later <= earlier for earlier, later in zip(epoch, epoch[1:])):
                raise ValueError(f"{label} epochs must strictly increase")
            defined = bool(epoch)
            panel_runs.append({
                "series": label, "shade": palette.run_shade(index),
                "train_hue": palette.TRAIN_FAMILY, "validation_hue": palette.VALIDATION_FAMILY,
                "defined": defined, "reason": None if defined else (run.get("reason") or "undefined"),
                "epoch": epoch, "train": train, "validation": validation})
        panels.append({"panel": name, "title": name, "x_label": "Epoch", "y_label": "Loss",
                       "runs": panel_runs})
    return {"schema": LEARNING_CURVE_SCHEMA, "kind": "learning_curve", "title": title,
            "layout": {"rows": 1, "columns": len(panels),
                       "panels": tuple(panel["panel"] for panel in panels)},
            "panels": panels}


def _axes_grid(axes, rows, columns):
    grid = np.asarray(axes, dtype=object)
    if grid.ndim == 0:
        grid = grid.reshape(1, 1)
    elif grid.ndim == 1:
        if rows == 1:
            grid = grid.reshape(1, -1)
        elif columns == 1:
            grid = grid.reshape(-1, 1)
        else:
            raise ValueError("one-dimensional axes cannot carry this figure layout")
    if grid.shape != (rows, columns):
        raise ValueError("axes grid does not match the figure layout")
    return grid


def _pyplot():
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    return plt


def _draw_bar_panel(panel, ax, fixed):
    positions = np.arange(len(panel["entries"]), dtype=float)
    for position, entry in zip(positions, panel["entries"]):
        line = palette.color(entry["hue"], palette.LINE_SHADE)
        fill = palette.color(entry["hue"], palette.FILL_SHADE)
        if entry["mean_defined"]:
            ax.bar(position, entry["mean"], width=.7, facecolor=fill,
                   edgecolor=line, linewidth=1.6, zorder=2)
            dots = [value for value, defined in zip(entry["values"], entry["defined"]) if defined]
            if dots:
                ax.plot([position] * len(dots), dots, linestyle="none", marker="o",
                        markersize=5, color=line,
                        markeredgecolor=fixed["dot_edge"],
                        markeredgewidth=.8, zorder=4)
        else:
            ax.plot([position], [.02], transform=ax.get_xaxis_transform(),
                    linestyle="none", marker="x", markersize=7,
                    color=fixed["undefined_marker"], zorder=4)
            ax.text(position, .04, "undefined", transform=ax.get_xaxis_transform(),
                    rotation=90, ha="center", va="bottom", fontsize=7,
                    color=fixed["undefined_label"])
    ax.axhline(0.0, color=fixed["grid"], linewidth=.8, zorder=1)
    ax.set_xticks(positions, [entry["series"] for entry in panel["entries"]],
                  rotation=20, ha="right", fontsize=8)
    ax.set_xlim(-.6, len(panel["entries"]) - .4)
    ax.grid(axis="y", color=fixed["grid"], alpha=.6, linewidth=.6)


def _draw_curve_panel(panel, ax, fixed):
    axis = np.asarray(panel["d_prime"], dtype=float)
    for entry in panel["entries"]:
        ax.fill_between(axis, np.asarray(entry["minimum"], dtype=float),
                        np.asarray(entry["maximum"], dtype=float),
                        facecolor=palette.color(entry["hue"], palette.FILL_SHADE),
                        edgecolor="none", zorder=1)
        ax.plot(axis, np.asarray(entry["mean"], dtype=float),
                color=palette.color(entry["hue"], palette.LINE_SHADE),
                linewidth=1.8, zorder=3)
    for entry in panel["entries"]:
        undefined = axis[~np.asarray(entry["mean_defined"], dtype=bool)]
        if len(undefined):
            ax.plot(undefined, np.full(len(undefined), .02),
                    transform=ax.get_xaxis_transform(), linestyle="none",
                    marker="x", markersize=6, color=fixed["undefined_marker"],
                    clip_on=False, zorder=5)
            ax.text(undefined[0], .04, "undefined", transform=ax.get_xaxis_transform(),
                    ha="left", va="bottom", fontsize=7,
                    color=fixed["undefined_label"])
    for entry in panel["entries"]:
        if "reference" not in entry:
            continue
        ax.plot(axis, np.asarray(entry["reference"], dtype=float),
                color=palette.color(entry["hue"], palette.REFERENCE_SHADE),
                linewidth=.85, linestyle="--", zorder=2)
        for marker in entry.get("d_star", ()):
            color = palette.color(entry["hue"], palette.MARKER_SHADE)
            ax.plot([marker["x"], marker["x"]], [0.0, marker["y"]],
                    color=color, linewidth=.7, linestyle=":", zorder=4)
            ax.plot([marker["x"]], [marker["y"]], linestyle="none", marker="o",
                    markersize=4.5, color=color, zorder=5)
    ax.set_xlim(float(axis.min()), float(axis.max()))
    ax.grid(axis="y", color=fixed["grid"], alpha=.6, linewidth=.6)


def _count_label(value):
    """A compact deterministic label for a confusion count."""
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if value >= 1_000:
        return f"{value / 1_000:.1f}k"
    return str(int(value))


def _draw_heatmap_panel(panel, ax, fixed):
    from matplotlib.colors import ListedColormap, LogNorm
    classes = tuple(panel["classes"])
    matrices = [np.asarray(entry["matrix"], dtype=float) for entry in panel["entries"]]
    if any(matrix.shape != (len(classes), len(classes)) for matrix in matrices):
        raise ValueError("a confusion heat map must be a square class matrix")
    vmax = max([1.0] + [float(matrix.max()) for matrix in matrices])
    cmap = ListedColormap(palette.heatmap_colors())
    norm = LogNorm(vmin=1.0, vmax=vmax)
    gap, width = 0.6, len(classes)
    offsets = [index * (width + gap) for index in range(len(matrices))]
    for entry, matrix, offset in zip(panel["entries"], matrices, offsets):
        ax.imshow(matrix, extent=(offset, offset + width, width, 0),
                  origin="upper", aspect="auto", interpolation="nearest",
                  cmap=cmap, norm=norm, zorder=2)
        ax.text(offset + width / 2, -0.35, entry["series"], ha="center", va="bottom",
                fontsize=8, color=fixed["text"], clip_on=False)
        for row in range(width):
            for column in range(width):
                value = float(matrix[row, column])
                light = value >= math.sqrt(vmax)
                ax.text(offset + column + .5, row + .5, _count_label(value),
                        ha="center", va="center", fontsize=6,
                        color=palette.color("gray", 8 if light else 0), zorder=3)
    ticks = [offset + .5 + column for offset in offsets for column in range(width)]
    labels = [classes[column] for _ in offsets for column in range(width)]
    ax.set_xticks(ticks, labels, rotation=30, ha="right", fontsize=7)
    ax.set_yticks([row + .5 for row in range(width)],
                  [name.replace("_", " ") for name in classes], fontsize=7)
    ax.set_xlim(0, offsets[-1] + width)
    ax.set_ylim(width, 0)
    ax.tick_params(length=0)


def draw_metrics(values, axes):
    """Draw the recorded direct or reliability panels; return those values.

    The panel kinds are the bar, the ``d'`` curve (optionally with a threshold
    marker and an ``F1@d`` reference line) and the three-class confusion heat map.
    The layout recorded in ``values`` drives the axes grid.
    """
    schema, kind = values.get("schema"), values.get("kind")
    if schema not in (METRICS_SCHEMA, INDIRECT_SCHEMA) or kind not in ("metrics",
                                                                      "metrics_indirect"):
        raise ValueError("unknown metrics plotted values")
    layout = values["layout"]
    grid = _axes_grid(axes, layout["rows"], layout["columns"])
    if tuple(layout["panels"]) != tuple(panel["panel"] for panel in values["panels"]):
        raise ValueError("recorded panels do not match the recorded layout")
    fixed = palette.role_colors(entry_hue=values["panels"][0]["entries"][0]["hue"], run_index=0)
    for panel, ax in zip(values["panels"], grid.flat):
        if panel["kind"] == "bar":
            _draw_bar_panel(panel, ax, fixed)
        elif panel["kind"] == "curve":
            _draw_curve_panel(panel, ax, fixed)
        elif panel["kind"] == "heatmap":
            _draw_heatmap_panel(panel, ax, fixed)
        else:
            raise ValueError(f"unknown panel kind: {panel['kind']!r}")
        ax.set(title=panel["title"], xlabel=panel["x_label"], ylabel=panel["y_label"])
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=8)
    return values


def draw_learning_curve(values, axes):
    """Draw one subplot per learnable model from its values; return those values."""
    if values.get("schema") != LEARNING_CURVE_SCHEMA or values.get("kind") != "learning_curve":
        raise ValueError("unknown learning-curve plotted values")
    grid = _axes_grid(axes, values["layout"]["rows"], values["layout"]["columns"])
    fixed = palette.role_colors(entry_hue=palette.TRAIN_FAMILY, run_index=0)
    for panel, ax in zip(values["panels"], grid.flat):
        for run in panel["runs"]:
            if not run["defined"]:
                continue
            ax.plot(run["epoch"], run["train"],
                    color=palette.color(run["train_hue"], run["shade"]),
                    linewidth=1.6, zorder=3)
            ax.plot(run["epoch"], run["validation"],
                    color=palette.color(run["validation_hue"], run["shade"]),
                    linewidth=1.6, zorder=3)
        missing = [run["series"] for run in panel["runs"] if not run["defined"]]
        if missing:
            ax.text(.98, .98, "undefined run(s): " + ", ".join(missing),
                    transform=ax.transAxes, ha="right", va="top", fontsize=8,
                    color=fixed["undefined_label"])
        ax.set(title=panel["title"], xlabel=panel["x_label"], ylabel=panel["y_label"])
        ax.grid(color=fixed["grid"], alpha=.6, linewidth=.6)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=9)
    return values


def _png(output):
    output = Path(output)
    if output.suffix != ".png":
        raise ValueError("a PNG output path is required")
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def _legend_handles(values):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    first = values["panels"][0]
    if first["kind"] == "bar":
        return [Patch(facecolor=palette.color(entry["hue"], palette.FILL_SHADE),
                      edgecolor=palette.color(entry["hue"], palette.LINE_SHADE),
                      linewidth=1.6, label=entry["series"]) for entry in first["entries"]]
    handles = [Line2D([0], [0], color=palette.color(entry["hue"], palette.LINE_SHADE),
                      linewidth=2, label=entry["series"]) for entry in first["entries"]]
    for panel in values["panels"]:
        for entry in panel["entries"]:
            if "reference" in entry:
                handles.append(Line2D(
                    [0], [0], color=palette.color(entry["hue"], palette.REFERENCE_SHADE),
                    linewidth=1, linestyle="--", label=f"{entry['series']} | F1@d"))
    return handles


def _save_metrics(values, output):
    output = _png(output)
    plt = _pyplot()
    layout = values["layout"]
    fig, axes = plt.subplots(layout["rows"], layout["columns"],
                             figsize=(4.2 * layout["columns"], 4.2 * layout["rows"]),
                             layout="constrained")
    try:
        draw_metrics(values, axes)
        handles = _legend_handles(values)
        fig.legend(handles=handles, loc="outside lower center", ncol=len(handles),
                   fontsize=9, frameon=False)
        fig.suptitle(values["title"], fontsize=13)
        fig.savefig(output, dpi=160)
    finally:
        plt.close(fig)
    return values


def _save_learning_curve(values, output):
    output = _png(output)
    plt = _pyplot()
    from matplotlib.lines import Line2D
    columns = values["layout"]["columns"]
    fig, axes = plt.subplots(1, columns, figsize=(7.5 * columns, 4.8),
                             layout="constrained", squeeze=False)
    try:
        draw_learning_curve(values, axes)
        handles = [Line2D([0], [0], color=palette.color(palette.TRAIN_FAMILY, palette.LINE_SHADE),
                          linewidth=2, label="Training"),
                   Line2D([0], [0], color=palette.color(palette.VALIDATION_FAMILY, palette.LINE_SHADE),
                          linewidth=2, label="Validation")]
        fig.legend(handles=handles, loc="outside lower center", ncol=2, fontsize=9,
                   frameon=False)
        fig.suptitle(values["title"], fontsize=13)
        fig.savefig(output, dpi=160)
    finally:
        plt.close(fig)
    return values


def render_metrics(entries, output, *, title=METRICS_TITLE, panel_titles=None):
    """Draw and save the direct-metric figure, returning every plotted value."""
    return _save_metrics(metrics_values(entries, title=title, panel_titles=panel_titles), output)


def render_indirect_metrics(entries, output, *, title=INDIRECT_TITLE, panel_titles=None):
    """Draw and save the reliability figure, returning every plotted value."""
    return _save_metrics(indirect_metrics_values(entries, title=title,
                                                 panel_titles=panel_titles), output)


def render_learning_curve(models, output, *, title=LEARNING_CURVE_TITLE):
    """Draw and save the learning-curve figure, returning every plotted value."""
    return _save_learning_curve(learning_curve_values(models, title=title), output)


def replay_metrics(path, output):
    """Rebuild a metrics PNG from its recorded plot-value CSV alone."""
    return _save_metrics(read_plot_values(path), output)


def replay_learning_curve(path, output):
    """Rebuild the learning-curve PNG from its recorded plot-value CSV alone."""
    return _save_learning_curve(read_plot_values(path), output)


def _blank(**fields):
    row = dict.fromkeys(PLOT_VALUE_COLUMNS, "")
    row.update({"x": float("nan"), "y": float("nan"), "defined": 1, "shade": "", "reason": ""})
    row.update(fields)
    return row


def _figure_row(values):
    layout = values["layout"]
    fields = {"record": "figure", "schema": values["schema"], "kind": values["kind"],
              "title": values["title"]}
    if values["kind"] in ("metrics", "metrics_indirect"):
        fields.update({"x": int(layout["rows"]), "y": int(layout["columns"])})
    return _blank(**fields)


def _panel_row(schema, kind, panel):
    return _blank(record="panel", schema=schema, kind=panel["kind"],
                  panel=panel["panel"], title=panel["title"],
                  x_label=panel["x_label"], y_label=panel["y_label"])


def _series_row(schema, kind, panel, series, hue, shade=palette.LINE_SHADE):
    return _blank(record="series", schema=schema, kind=kind, panel=panel,
                  series=series, hue=hue, shade=shade)


def _bar_rows(schema, kind, panel):
    rows = [_series_row(schema, kind, panel["panel"], entry["series"], entry["hue"])
            for entry in panel["entries"]]
    for position, entry in enumerate(panel["entries"]):
        rows.append(_blank(record="point", schema=schema, kind=kind,
                           panel=panel["panel"], series=entry["series"],
                           role="bar_mean", x=float(position), y=entry["mean"],
                           defined=int(entry["mean_defined"]),
                           reason=entry["reason"] or ""))
        for run, value, defined in zip(entry["run_labels"], entry["values"],
                                       entry["defined"]):
            rows.append(_blank(record="point", schema=schema, kind=kind,
                               panel=panel["panel"], series=entry["series"],
                               role="bar_run", run=run, x=float(position), y=value,
                               defined=int(defined)))
    return rows


def _curve_rows(schema, kind, panel):
    rows = [_series_row(schema, kind, panel["panel"], entry["series"], entry["hue"])
            for entry in panel["entries"]]
    for entry in panel["entries"]:
        for index, x in enumerate(panel["d_prime"]):
            undefined_reason = entry["reason"] if not entry["mean_defined"][index] else ""
            for role, column in (("curve_mean", "mean"), ("curve_min", "minimum"),
                                 ("curve_max", "maximum")):
                rows.append(_blank(record="point", schema=schema, kind=kind,
                                   panel=panel["panel"], series=entry["series"],
                                   role=role, x=float(x), y=entry[column][index],
                                   defined=int(entry["mean_defined"][index]),
                                   reason=undefined_reason))
        for run, curve in zip(entry["run_labels"], entry["runs"]):
            for x, value in zip(panel["d_prime"], curve):
                rows.append(_blank(record="point", schema=schema, kind=kind,
                                   panel=panel["panel"], series=entry["series"],
                                   role="curve_run", run=run, x=float(x), y=value,
                                   defined=int(bool(np.isfinite(value)))))
        if "reference" in entry:
            reference = f"{entry['series']} | F1@d"
            rows.append(_series_row(schema, kind, panel["panel"], reference, entry["hue"]))
            for x, value in zip(panel["d_prime"], entry["reference"]):
                rows.append(_blank(record="point", schema=schema, kind=kind,
                                   panel=panel["panel"], series=reference,
                                   role="reference", x=float(x), y=value, defined=1))
        for marker in entry.get("d_star", ()):
            rows.append(_blank(record="point", schema=schema, kind=kind,
                               panel=panel["panel"], series=entry["series"],
                               role="marker", run=marker["run"], x=marker["x"],
                               y=marker["y"], defined=1))
    return rows


def _heatmap_rows(schema, kind, panel):
    rows = []
    for entry in panel["entries"]:
        rows.append(_series_row(schema, kind, panel["panel"], entry["series"], entry["hue"]))
        for row_index, predicted in enumerate(entry["classes"]):
            for column_index, measured in enumerate(entry["classes"]):
                rows.append(_blank(record="point", schema=schema, kind=kind,
                                   panel=panel["panel"], series=entry["series"],
                                   role="heat_cell", x=float(column_index),
                                   y=float(entry["matrix"][row_index][column_index]),
                                   defined=1, reason=f"{predicted}|{measured}"))
    return rows


def _learning_rows(values):
    schema, kind = values["schema"], values["kind"]
    rows = [_figure_row(values)]
    for panel in values["panels"]:
        rows.append(_blank(record="panel", schema=schema, kind=kind,
                           panel=panel["panel"], title=panel["title"],
                           x_label=panel["x_label"], y_label=panel["y_label"]))
        for run in panel["runs"]:
            for split, hue, trajectory in (
                    ("train", run["train_hue"], run["train"]),
                    ("validation", run["validation_hue"], run["validation"])):
                series = f"{run['series']}|{split}"
                rows.append(_blank(record="series", schema=schema, kind=kind,
                                   panel=panel["panel"], series=series, hue=hue,
                                   shade=run["shade"], run=run["series"],
                                   defined=int(run["defined"]),
                                   reason=run["reason"] or ""))
                for x, y in zip(run["epoch"], trajectory):
                    rows.append(_blank(record="point", schema=schema, kind=kind,
                                       panel=panel["panel"], series=series,
                                       role="line", run=run["series"], x=float(x),
                                       y=y, defined=int(run["defined"])))
    return rows


def plot_rows(values):
    """Flatten the plotted values into the recorded long-form CSV rows."""
    kind = values["kind"]
    if kind == "learning_curve":
        return _learning_rows(values)
    if kind not in ("metrics", "metrics_indirect"):
        raise ValueError(f"unknown plotted-values kind: {kind!r}")
    schema = values["schema"]
    rows = [_figure_row(values)]
    for panel in values["panels"]:
        rows.append(_panel_row(schema, kind, panel))
        if panel["kind"] == "bar":
            rows.extend(_bar_rows(schema, kind, panel))
        elif panel["kind"] == "curve":
            rows.extend(_curve_rows(schema, kind, panel))
        elif panel["kind"] == "heatmap":
            rows.extend(_heatmap_rows(schema, kind, panel))
        else:
            raise ValueError(f"unknown panel kind: {panel['kind']!r}")
    return rows


def write_plot_values(values, path):
    """Write the complete plotted values as the recorded long-form CSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(plot_rows(values), columns=list(PLOT_VALUE_COLUMNS)).to_csv(path, index=False)
    return path


def _number(text):
    return float(text) if text not in ("", None) else float("nan")


def _recorded_layout(figure):
    rows, columns = _number(figure["x"]), _number(figure["y"])
    if not (math.isfinite(rows) and math.isfinite(columns)):
        raise ValueError("the metrics figure layout is not recorded")
    if not (rows.is_integer() and columns.is_integer()) or rows < 1 or columns < 1:
        raise ValueError("a recorded metrics layout needs positive whole rows and columns")
    return int(rows), int(columns)


def _panel_points(rows, name):
    return [row for row in rows if row["record"] == "point" and row["panel"] == name]


def _bar_panel_from_rows(rows, panel_row):
    name = panel_row["panel"]
    series_rows = [row for row in rows if row["record"] == "series" and row["panel"] == name]
    points = _panel_points(rows, name)
    if not series_rows:
        raise ValueError(f"panel {name!r} has no series")
    entries = []
    for series in series_rows:
        means = [row for row in points if row["role"] == "bar_mean"
                 and row["series"] == series["series"]]
        runs = [row for row in points if row["role"] == "bar_run"
                and row["series"] == series["series"]]
        if len(means) != 1:
            raise ValueError(f"panel {name!r} needs one bar mean per series")
        entries.append({
            "series": series["series"], "hue": series["hue"],
            "run_labels": tuple(row["run"] for row in runs),
            "values": tuple(_number(row["y"]) for row in runs),
            "defined": tuple(row["defined"] == "1" for row in runs),
            "mean": _number(means[0]["y"]),
            "mean_defined": means[0]["defined"] == "1",
            "reason": means[0]["reason"] or None,
            "undefined_runs": tuple(index for index, row in enumerate(runs)
                                    if row["defined"] != "1")})
    return {"panel": name, "kind": "bar", "title": panel_row["title"],
            "x_label": panel_row["x_label"], "y_label": panel_row["y_label"],
            "entries": entries}


def _curve_panel_from_rows(rows, panel_row):
    name = panel_row["panel"]
    series_rows = [row for row in rows if row["record"] == "series" and row["panel"] == name]
    points = _panel_points(rows, name)
    entry_series = [series for series in series_rows
                    if any(row["role"] == "curve_mean" and row["series"] == series["series"]
                           for row in points)]
    if not entry_series:
        raise ValueError(f"panel {name!r} has no curve series")
    entries, axis = [], None
    for series in entry_series:
        def ordered(role):
            return sorted([row for row in points if row["role"] == role
                           and row["series"] == series["series"]],
                          key=lambda row: _number(row["x"]))
        mean_rows, min_rows, max_rows = ordered("curve_mean"), ordered("curve_min"), ordered("curve_max")
        x_axis = tuple(_number(row["x"]) for row in mean_rows)
        _validate_d_prime(x_axis)
        if axis is None:
            axis = x_axis
        elif x_axis != axis:
            raise ValueError(f"panel {name!r} needs one shared d' axis")
        if len(min_rows) != len(x_axis) or len(max_rows) != len(x_axis):
            raise ValueError(f"panel {name!r} needs one min and max per d'")
        run_order, runs = [], {}
        for row in points:
            if row["role"] != "curve_run" or row["series"] != series["series"]:
                continue
            if row["run"] not in runs:
                runs[row["run"]] = []
                run_order.append(row["run"])
            runs[row["run"]].append(row)
        runs = {label: tuple(_number(row["y"]) for row in sorted(part, key=lambda row: _number(row["x"])))
                for label, part in runs.items()}
        if any(len(runs[label]) != len(x_axis) for label in run_order):
            raise ValueError(f"panel {name!r} needs one curve point per run and d'")
        reference_name = f"{series['series']} | F1@d"
        reference = [row for row in points if row["role"] == "reference"
                     and row["series"] == reference_name]
        entry = {
            "series": series["series"], "hue": series["hue"],
            "run_labels": tuple(run_order),
            "runs": tuple(runs[label] for label in run_order),
            "mean": tuple(_number(row["y"]) for row in mean_rows),
            "minimum": tuple(_number(row["y"]) for row in min_rows),
            "maximum": tuple(_number(row["y"]) for row in max_rows),
            "defined_count": tuple(sum(1 for label in run_order
                                       if np.isfinite(runs[label][index]))
                                   for index in range(len(x_axis))),
            "mean_defined": tuple(row["defined"] == "1" for row in mean_rows),
            "reason": next((row["reason"] for row in mean_rows if row["reason"]), None)}
        if reference:
            ordered_reference = sorted(reference, key=lambda row: _number(row["x"]))
            if tuple(_number(row["x"]) for row in ordered_reference) != x_axis:
                raise ValueError(f"panel {name!r} needs one reference point per d'")
            entry["reference"] = tuple(_number(row["y"]) for row in ordered_reference)
        markers = [row for row in points if row["role"] == "marker"
                   and row["series"] == series["series"]]
        if markers:
            entry["d_star"] = tuple({"run": row["run"], "x": _number(row["x"]),
                                     "y": _number(row["y"])} for row in markers)
        entries.append(entry)
    return {"panel": name, "kind": "curve", "title": panel_row["title"],
            "x_label": panel_row["x_label"], "y_label": panel_row["y_label"],
            "d_prime": axis, "entries": entries}


def _heatmap_panel_from_rows(rows, panel_row):
    name = panel_row["panel"]
    series_rows = [row for row in rows if row["record"] == "series" and row["panel"] == name]
    points = _panel_points(rows, name)
    if not series_rows:
        raise ValueError(f"panel {name!r} has no series")
    classes = None
    entries = []
    for series in series_rows:
        cells = [row for row in points if row["role"] == "heat_cell"
                 and row["series"] == series["series"]]
        if not cells:
            raise ValueError(f"panel {name!r} has no heat cells for {series['series']!r}")
        predicted, measured, counts = [], [], {}
        for row in cells:
            labels = row["reason"].split("|")
            if len(labels) != 2 or not all(labels):
                raise ValueError(f"panel {name!r} needs a predicted|measured cell label")
            if labels[0] not in predicted:
                predicted.append(labels[0])
            if labels[1] not in measured:
                measured.append(labels[1])
            counts[(labels[0], labels[1])] = int(_number(row["y"]))
        if predicted != measured:
            raise ValueError(f"panel {name!r} classes must agree on both axes")
        if classes is None:
            classes = tuple(predicted)
        elif classes != tuple(predicted):
            raise ValueError(f"panel {name!r} entries must share one class axis")
        matrix = tuple(tuple(counts[(row, column)] for column in classes)
                       for row in classes)
        entries.append({"series": series["series"], "hue": series["hue"],
                        "classes": classes, "matrix": matrix})
    return {"panel": name, "kind": "heatmap", "title": panel_row["title"],
            "x_label": panel_row["x_label"], "y_label": panel_row["y_label"],
            "classes": classes, "entries": entries}


def _metrics_from_rows(rows, figure):
    row_count, column_count = _recorded_layout(figure)
    panels = []
    for panel_row in [row for row in rows if row["record"] == "panel"]:
        kind = panel_row["kind"]
        if kind == "bar":
            panels.append(_bar_panel_from_rows(rows, panel_row))
        elif kind == "curve":
            panels.append(_curve_panel_from_rows(rows, panel_row))
        elif kind == "heatmap":
            panels.append(_heatmap_panel_from_rows(rows, panel_row))
        else:
            raise ValueError(f"unknown recorded panel kind: {kind!r}")
    if len(panels) != row_count * column_count:
        raise ValueError("recorded panel count does not fill the recorded layout")
    return {"schema": figure["schema"], "kind": figure["kind"], "title": figure["title"],
            "layout": {"rows": row_count, "columns": column_count,
                       "panels": tuple(panel["panel"] for panel in panels)},
            "panels": panels}


def _learning_from_rows(rows, figure):
    panels = []
    for panel_row in [row for row in rows if row["record"] == "panel"]:
        name = panel_row["panel"]
        series_rows = [row for row in rows if row["record"] == "series" and row["panel"] == name]
        points = [row for row in rows if row["record"] == "point" and row["panel"] == name]
        order, runs = [], {}
        for row in series_rows:
            if row["series"].endswith("|train"):
                split, label = "train", row["series"][:-len("|train")]
            elif row["series"].endswith("|validation"):
                split, label = "validation", row["series"][:-len("|validation")]
            else:
                raise ValueError(f"unknown learning-curve series: {row['series']!r}")
            if label not in runs:
                runs[label] = {"series": label, "shade": int(row["shade"]),
                               "defined": row["defined"] == "1", "reason": row["reason"] or None}
                order.append(label)
            run = runs[label]
            run[f"{split}_hue"] = row["hue"]
            trajectory = sorted([point for point in points if point["role"] == "line"
                                 and point["series"] == row["series"]],
                                key=lambda point: _number(point["x"]))
            if split == "train":
                run["epoch"] = tuple(_number(point["x"]) for point in trajectory)
            run[split] = tuple(_number(point["y"]) for point in trajectory)
        panel_runs = []
        for label in order:
            run = runs[label]
            if not run["defined"]:
                run["epoch"] = ()
                run["train"] = ()
                run["validation"] = ()
            panel_runs.append({"series": run["series"], "shade": run["shade"],
                               "train_hue": run["train_hue"],
                               "validation_hue": run["validation_hue"],
                               "defined": run["defined"], "reason": run["reason"],
                               "epoch": run["epoch"], "train": run["train"],
                               "validation": run["validation"]})
        panels.append({"panel": name, "title": panel_row["title"],
                       "x_label": panel_row["x_label"], "y_label": panel_row["y_label"],
                       "runs": panel_runs})
    return {"schema": figure["schema"], "kind": "learning_curve", "title": figure["title"],
            "layout": {"rows": 1, "columns": len(panels),
                       "panels": tuple(panel["panel"] for panel in panels)},
            "panels": panels}


def read_plot_values(path):
    """Rebuild the plotted values from their recorded long-form CSV alone."""
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    if tuple(frame.columns) != PLOT_VALUE_COLUMNS:
        raise ValueError(f"unexpected plot-value columns: {tuple(frame.columns)}")
    rows = frame.to_dict("records")
    figures = [row for row in rows if row["record"] == "figure"]
    if len(figures) != 1:
        raise ValueError("exactly one figure record is required")
    figure = figures[0]
    if figure["kind"] in ("metrics", "metrics_indirect"):
        expected = METRICS_SCHEMA if figure["kind"] == "metrics" else INDIRECT_SCHEMA
        if figure["schema"] != expected:
            raise ValueError(f"unknown {figure['kind']} plot-value schema: "
                             f"{figure['schema']!r}")
        return _metrics_from_rows(rows, figure)
    if figure["kind"] == "learning_curve" and figure["schema"] == LEARNING_CURVE_SCHEMA:
        return _learning_from_rows(rows, figure)
    raise ValueError(f"unknown plot-value source: {figure['kind']!r}")
