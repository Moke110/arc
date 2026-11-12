"""Stage 2 plotting APIs: the tenv retrieval panel and the ADR 0063 comparison figures.

Every Stage 2 comparison drawing call takes only recorded values and returns the
complete plotted values it drew, so a step renderer can rebuild a PNG from its
``figure/data/*.csv`` alone (ADR 0062, ADR 0063).  Colours come only from the
stage palette in :mod:`bases.palette`; undefined values are NaN with a reason and
are drawn explicitly, never as zero.
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
# ADR 0063 comparison figures: the eight-panel metrics figure and the learning
# curve.  Everything below is presentation of already-computed metric values; no
# scientific quantity is recomputed here.
# ---------------------------------------------------------------------------
METRICS_LAYOUT = (("PCC", "Spearman", "RMSE", "R2"),
                  ("Precision@k", "Spearman@k", "max_k@P80S40", "Recall@max_k"))
CURVE_METRICS = ("Precision@k", "Spearman@k")
BAR_METRICS = ("PCC", "Spearman", "RMSE", "R2", "max_k@P80S40", "Recall@max_k")
PANEL_ORDER = tuple(metric for row in METRICS_LAYOUT for metric in row)
METRICS_TITLE = "Tenv model comparison"
LEARNING_CURVE_TITLE = "Learning curve"
METRICS_SCHEMA = "stage2-tenv-metrics-values-v1"
LEARNING_CURVE_SCHEMA = "stage2-tenv-learning-curve-values-v1"
# Every recorded plot-value CSV has exactly these columns, whatever the figure.
PLOT_VALUE_COLUMNS = ("record", "schema", "kind", "panel", "title", "x_label", "y_label",
                      "series", "hue", "shade", "role", "run", "x", "y", "defined", "reason")
_METRIC_INPUT = {"PCC": ("direct", "PCC"), "Spearman": ("direct", "Spearman"),
                 "RMSE": ("direct", "RMSE"), "R2": ("direct", "R2"),
                 "max_k@P80S40": ("indirect", "max_k"),
                 "Recall@max_k": ("indirect", "recall")}
_CURVE_INPUT = {"Precision@k": "precision", "Spearman@k": "spearman"}


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
        raise ValueError(f"{label} must match the shared k axis")
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


def _panel_titles(panel_titles):
    titles = {metric: metric for metric in PANEL_ORDER}
    if panel_titles:
        unknown = sorted(set(panel_titles) - set(PANEL_ORDER))
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
            direct, indirect = run.get("direct") or {}, run.get("indirect") or {}
            if not isinstance(direct, dict) or not isinstance(indirect, dict):
                raise ValueError("run direct and indirect values must be mappings")
            normalised_runs.append({"run": run_label, "direct": direct,
                                    "indirect": indirect,
                                    "reasons": dict(run.get("reasons") or {})})
        normalised.append({"label": label, "hue": hue, "runs": normalised_runs})
    return normalised


def _k_axis(entries):
    """Return the one shared indirect k axis every run must carry."""
    axis = None
    for entry in entries:
        for run in entry["runs"]:
            if "k" not in run["indirect"]:
                raise ValueError("every run needs an indirect k axis")
            axis_values = run["indirect"]["k"]
            if not isinstance(axis_values, (list, tuple)) or not axis_values:
                raise ValueError("every run needs a nonempty indirect k axis")
            values = _vector(axis_values, len(axis_values), "the indirect k axis")
            if any(not (item > 0 and item.is_integer()) for item in values):
                raise ValueError("the indirect k axis must be positive whole numbers")
            if any(later <= earlier for earlier, later in zip(values, values[1:])):
                raise ValueError("the indirect k axis must strictly increase")
            if axis is None:
                axis = values
            elif values != axis:
                raise ValueError("every run must share one indirect k axis")
    return axis


def _run_reason(entry, metric):
    return next((run["reasons"].get(metric) for run in entry["runs"]
                 if run["reasons"].get(metric)), None)


def metrics_values(entries, *, title=METRICS_TITLE, panel_titles=None):
    """Compute the complete eight-panel plotted values without drawing anything.

    Indirect panels draw ``k`` with the cross-run mean line and min-max band; the
    direct and threshold panels draw bars of the cross-run mean with one dot per
    run.  A value defined for no run is NaN with a reason.
    """
    title = _validate_title(title)
    titles = _panel_titles(panel_titles)
    entries = _validate_entries(entries)
    k_axis = _k_axis(entries)
    panels = []
    for metric in PANEL_ORDER:
        if metric in CURVE_METRICS:
            field, panel_entries = _CURVE_INPUT[metric], []
            for entry in entries:
                runs = tuple(_vector(run["indirect"].get(field), len(k_axis),
                                     f"{metric} for {run['run']}") for run in entry["runs"])
                matrix = np.vstack([np.asarray(run, dtype=np.float64) for run in runs])
                defined = np.isfinite(matrix)
                count = defined.sum(axis=0)
                mean = np.divide(np.where(defined, matrix, 0.0).sum(axis=0), count,
                                 out=np.full(len(k_axis), np.nan), where=count > 0)
                minimum = np.where(count > 0, np.min(np.where(defined, matrix, np.inf), axis=0), np.nan)
                maximum = np.where(count > 0, np.max(np.where(defined, matrix, -np.inf), axis=0), np.nan)
                mean_defined = tuple(bool(item) for item in count > 0)
                reason = _run_reason(entry, metric)
                if reason is None and not all(mean_defined):
                    reason = "undefined"
                panel_entries.append({
                    "series": entry["label"], "hue": entry["hue"],
                    "run_labels": tuple(run["run"] for run in entry["runs"]),
                    "runs": runs, "mean": tuple(float(item) for item in mean),
                    "minimum": tuple(float(item) for item in minimum),
                    "maximum": tuple(float(item) for item in maximum),
                    "defined_count": tuple(int(item) for item in count),
                    "mean_defined": mean_defined, "reason": reason})
            panels.append({"panel": metric, "kind": "curve", "title": titles[metric],
                           "x_label": "k", "y_label": metric, "k": k_axis,
                           "entries": panel_entries})
            continue
        section, key = _METRIC_INPUT[metric]
        panel_entries = []
        for entry in entries:
            values = tuple(_scalar(run[section].get(key), f"{metric} for {run['run']}")
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
        panels.append({"panel": metric, "kind": "bar", "title": titles[metric],
                       "x_label": "Entry", "y_label": metric, "k": None,
                       "entries": panel_entries})
    return {"schema": METRICS_SCHEMA, "kind": "metrics", "title": title,
            "layout": {"rows": 2, "columns": 4, "panels": PANEL_ORDER},
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


def draw_metrics(values, axes):
    """Draw the eight ADR 0063 panels from their plotted values; return those values."""
    if values.get("schema") != METRICS_SCHEMA or values.get("kind") != "metrics":
        raise ValueError("unknown metrics plotted values")
    grid = _axes_grid(axes, values["layout"]["rows"], values["layout"]["columns"])
    fixed = palette.role_colors(entry_hue=values["panels"][0]["entries"][0]["hue"], run_index=0)
    for panel, ax in zip(values["panels"], grid.flat):
        if panel["kind"] == "bar":
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
        else:
            k = np.asarray(panel["k"], dtype=float)
            for entry in panel["entries"]:
                ax.fill_between(k, np.asarray(entry["minimum"], dtype=float),
                                np.asarray(entry["maximum"], dtype=float),
                                facecolor=palette.color(entry["hue"], palette.FILL_SHADE),
                                edgecolor="none", zorder=1)
                ax.plot(k, np.asarray(entry["mean"], dtype=float),
                        color=palette.color(entry["hue"], palette.LINE_SHADE),
                        linewidth=1.8, zorder=3)
            for entry in panel["entries"]:
                undefined = k[~np.asarray(entry["mean_defined"], dtype=bool)]
                if len(undefined):
                    ax.plot(undefined, np.full(len(undefined), .02),
                            transform=ax.get_xaxis_transform(), linestyle="none",
                            marker="x", markersize=6, color=fixed["undefined_marker"],
                            clip_on=False, zorder=5)
                    ax.text(undefined[0], .04, "undefined", transform=ax.get_xaxis_transform(),
                            ha="left", va="bottom", fontsize=7,
                            color=fixed["undefined_label"])
            ax.set_xlim(float(k.min()), float(k.max()))
            ax.grid(axis="y", color=fixed["grid"], alpha=.6, linewidth=.6)
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


def _save_metrics(values, output):
    output = _png(output)
    plt = _pyplot()
    from matplotlib.patches import Patch
    fig, axes = plt.subplots(values["layout"]["rows"], values["layout"]["columns"],
                             figsize=(16, 7.5), layout="constrained")
    try:
        draw_metrics(values, axes)
        entries = values["panels"][0]["entries"]
        handles = [Patch(facecolor=palette.color(entry["hue"], palette.FILL_SHADE),
                         edgecolor=palette.color(entry["hue"], palette.LINE_SHADE),
                         linewidth=1.6, label=entry["series"]) for entry in entries]
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
    """Draw and save the metrics figure, returning every plotted value."""
    return _save_metrics(metrics_values(entries, title=title, panel_titles=panel_titles), output)


def render_learning_curve(models, output, *, title=LEARNING_CURVE_TITLE):
    """Draw and save the learning-curve figure, returning every plotted value."""
    return _save_learning_curve(learning_curve_values(models, title=title), output)


def replay_metrics(path, output):
    """Rebuild the metrics PNG from its recorded plot-value CSV alone."""
    return _save_metrics(read_plot_values(path), output)


def replay_learning_curve(path, output):
    """Rebuild the learning-curve PNG from its recorded plot-value CSV alone."""
    return _save_learning_curve(read_plot_values(path), output)


def _blank(**fields):
    row = dict.fromkeys(PLOT_VALUE_COLUMNS, "")
    row.update({"x": float("nan"), "y": float("nan"), "defined": 1, "shade": "", "reason": ""})
    row.update(fields)
    return row


def plot_rows(values):
    """Flatten the plotted values into the recorded long-form CSV rows."""
    schema, kind = values["schema"], values["kind"]
    rows = [_blank(record="figure", schema=schema, kind=kind, title=values["title"])]
    if kind == "metrics":
        for panel in values["panels"]:
            rows.append(_blank(record="panel", schema=schema, kind=panel["kind"],
                               panel=panel["panel"], title=panel["title"],
                               x_label=panel["x_label"], y_label=panel["y_label"]))
            for entry in panel["entries"]:
                rows.append(_blank(record="series", schema=schema, kind=kind,
                                   panel=panel["panel"], series=entry["series"],
                                   hue=entry["hue"], shade=palette.LINE_SHADE))
            if panel["kind"] == "bar":
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
                                           role="bar_run", run=run, x=float(position),
                                           y=value, defined=int(defined)))
                continue
            for entry in panel["entries"]:
                for index, k in enumerate(panel["k"]):
                    undefined_reason = entry["reason"] if not entry["mean_defined"][index] else ""
                    for role, column in (("curve_mean", "mean"), ("curve_min", "minimum"),
                                         ("curve_max", "maximum")):
                        rows.append(_blank(record="point", schema=schema, kind=kind,
                                           panel=panel["panel"], series=entry["series"],
                                           role=role, x=float(k),
                                           y=entry[column][index],
                                           defined=int(entry["mean_defined"][index]),
                                           reason=undefined_reason))
                for run, curve in zip(entry["run_labels"], entry["runs"]):
                    for k, value in zip(panel["k"], curve):
                        rows.append(_blank(record="point", schema=schema, kind=kind,
                                           panel=panel["panel"], series=entry["series"],
                                           role="curve_run", run=run, x=float(k), y=value,
                                           defined=int(bool(np.isfinite(value)))))
    elif kind == "learning_curve":
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
    else:
        raise ValueError(f"unknown plotted-values kind: {kind!r}")
    return rows


def write_plot_values(values, path):
    """Write the complete plotted values as the recorded long-form CSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(plot_rows(values), columns=list(PLOT_VALUE_COLUMNS)).to_csv(path, index=False)
    return path


def _number(text):
    return float(text) if text not in ("", None) else float("nan")


def _metrics_from_rows(rows, figure):
    panels = []
    for panel_row in [row for row in rows if row["record"] == "panel"]:
        name = panel_row["panel"]
        series_rows = [row for row in rows if row["record"] == "series" and row["panel"] == name]
        points = [row for row in rows if row["record"] == "point" and row["panel"] == name]
        if not series_rows:
            raise ValueError(f"panel {name!r} has no series")
        entries = []
        if panel_row["kind"] == "bar":
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
            panel = {"panel": name, "kind": "bar", "title": panel_row["title"],
                     "x_label": panel_row["x_label"], "y_label": panel_row["y_label"],
                     "k": None, "entries": entries}
        else:
            axis = None
            for series in series_rows:
                def ordered(role):
                    return sorted([row for row in points if row["role"] == role
                                   and row["series"] == series["series"]],
                                  key=lambda row: _number(row["x"]))
                mean_rows, min_rows, max_rows = ordered("curve_mean"), ordered("curve_min"), ordered("curve_max")
                k = tuple(_number(row["x"]) for row in mean_rows)
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
                if any(len(runs[label]) != len(k) for label in run_order):
                    raise ValueError(f"panel {name!r} needs one curve point per run and k")
                if axis is None:
                    axis = k
                elif k != axis:
                    raise ValueError(f"panel {name!r} needs one shared k axis")
                entries.append({
                    "series": series["series"], "hue": series["hue"],
                    "run_labels": tuple(run_order),
                    "runs": tuple(runs[label] for label in run_order),
                    "mean": tuple(_number(row["y"]) for row in mean_rows),
                    "minimum": tuple(_number(row["y"]) for row in min_rows),
                    "maximum": tuple(_number(row["y"]) for row in max_rows),
                    "defined_count": tuple(sum(1 for label in run_order
                                               if np.isfinite(runs[label][index]))
                                           for index in range(len(k))),
                    "mean_defined": tuple(row["defined"] == "1" for row in mean_rows),
                    "reason": next((row["reason"] for row in mean_rows if row["reason"]), None)})
            panel = {"panel": name, "kind": "curve", "title": panel_row["title"],
                     "x_label": panel_row["x_label"], "y_label": panel_row["y_label"],
                     "k": axis, "entries": entries}
        panels.append(panel)
    return {"schema": figure["schema"], "kind": "metrics", "title": figure["title"],
            "layout": {"rows": 2, "columns": 4,
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
    if figure["kind"] == "metrics" and figure["schema"] == METRICS_SCHEMA:
        return _metrics_from_rows(rows, figure)
    if figure["kind"] == "learning_curve" and figure["schema"] == LEARNING_CURVE_SCHEMA:
        return _learning_from_rows(rows, figure)
    raise ValueError(f"unknown plot-value source: {figure['kind']!r}")
