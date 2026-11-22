"""Draw the accepted baseline comparison in baseline-first order."""
import argparse
import json
from pathlib import Path
from runpy import run_path

FIGURE = Path(__file__).resolve().parents[1]
HISTORY = FIGURE.parent / 'execution/historical/report_axis_local'
MEAN_GRAY = '#474747'
MLP_GRAY = '#898989'


def render(data, output, *, history=HISTORY):
    previous = run_path(str(history / 'baselines/script/render.py'))
    base = previous['load_reference'](history / 'script/render_axis_local.py')
    read_rows = base.read_rows

    def read_plot_rows(path):
        return read_rows(path.with_name(path.name.removeprefix('round1_')))

    base.read_rows = read_plot_rows
    base.TITLES = tuple(title.rsplit(' ', 1)[0] for title in base.TITLES)
    tables = {name: base.read_rows(data / name)
              for name in (*previous['METRIC_TABLES'], previous['CURVE_TABLE'])}
    previous['validate_tables'](base, tables)
    index = base.by_key(base.read_rows(data / previous['AXIS_TABLE']),
                        'figure', 'source', 'scenario', 'panel', 'ordinate')
    axes = {key: base.AxisScale(json.loads(record['segments'])) for key, (record,) in index.items()}
    for (kind, _, _, metric, _), scale in axes.items():
        if kind != 'overall' or metric == 'r2':
            continue
        first = scale.segments[0]
        if first['upper'] == 0:
            # An isolated zero would leave every bar floating above the axis break.
            scale.segments = [base.linear_segment(0., scale.segments[-1]['upper'], show_zero=True)]
        else:
            first['lower'] = 0.
            first['ticks'] = [tick for tick in first['ticks'] if tick >= 0]

    mean_cell, mlp = previous['MEAN_CELL'], previous['ECFP4_MLP']
    save = base.save_figure

    def save_styled(figure, destination, name):
        title = figure._suptitle
        heading = title.get_text().removeprefix('Round 1 · ')
        title.set_text(heading[0].upper() + heading[1:])
        if name == base.PNG_NAMES[0]:
            for panel in figure.axes:
                labels = [tick.get_text().split('\n', 1)[0] for tick in panel.get_xticklabels()]
                panel.set_xticks(panel.get_xticks(), labels, rotation=32, ha='right', fontsize=8)
        if name == base.PNG_NAMES[1]:
            palette = dict(zip(base.TRAIN_COLORS, ('#353535', '#626262', '#8b8b8b'), strict=True))
            palette.update(zip(base.VALID_COLORS, ('#494949', '#757575', '#9b9b9b'), strict=True))
            for panel in figure.axes:
                if panel.get_title() == mlp:
                    for axis in panel.get_shared_x_axes().get_siblings(panel):
                        for artist in (*axis.lines, axis.yaxis.label):
                            color = artist.get_color()
                            if color in palette:
                                artist.set_color(palette[color])
                        axis.tick_params(axis='y', colors=axis.yaxis.label.get_color())
        if name == base.PNG_NAMES[2]:
            handles, labels = figure.axes[0].get_legend_handles_labels()
            figure.legends[0].remove()
            figure.legend(handles, labels, loc='lower center', ncol=len(labels),
                          bbox_to_anchor=(.5, .018))
        save(figure, destination, name.removeprefix('round1_'))

    base.save_figure = save_styled
    base.METHODS = (mlp, *base.METHODS)
    base.LABELS = (mlp, *base.LABELS)
    base.COLORS = (MLP_GRAY, *base.COLORS)
    output.mkdir(parents=True, exist_ok=True)
    with base.plt.rc_context({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False}):
        previous['draw_learning'](base, data, axes, output)
        base.METHODS = (mean_cell, *base.METHODS)
        base.LABELS = (mean_cell, *base.LABELS)
        base.COLORS = (MEAN_GRAY, *base.COLORS)
        base.draw_overall(data, axes, output)
        base.draw_doses(data, axes, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=FIGURE / 'data')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    render(args.data, args.output)
