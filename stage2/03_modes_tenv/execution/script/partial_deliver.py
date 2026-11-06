"""Publish checked partial figures at the fixed interim/final visible paths."""
import argparse
from pathlib import Path
import shutil
import textwrap
from types import FunctionType
from artifacts import read_json, rows, sha256, write_json
from partial_benchmark import PartialPlot, WARNING, POLICY
import plots


def deliver_partial(snapshot, publish):
    snapshot, publish = Path(snapshot).resolve(), Path(publish).resolve()
    evidence = read_json(snapshot / 'evaluation.json')
    manifest = read_json(snapshot / 'predictions.json')
    if evidence.get('scope') != 'partial_finished_models' or manifest.get('sealed') is not False:
        raise ValueError('checked nonsealed partial evaluation required')
    if publish != Path(manifest['source_root']).parent:
        raise ValueError('publish only to the fixed step output directory')
    for name, digest in evidence['artifacts'].items():
        if sha256(snapshot / name) != digest:
            raise ValueError('partial publication checksum mismatch: ' + name)
    summaries = rows(snapshot / 'summary.csv')
    counts = '; '.join(f"{r['model']} {r['completed_runs']}/9" for r in summaries)
    dated = manifest['captured_at_utc'] + ' | completed runs: ' + counts

    class DatedPlot(PartialPlot):
        def subplots(self, *args, **kwargs):
            figure, axes = super().subplots(*args, **kwargs)
            figure.text(.5, .035, textwrap.fill(dated, 115), ha='center', va='bottom', fontsize=7,
                        bbox={'facecolor': 'white', 'alpha': .9, 'edgecolor': 'none'})
            figure.text(.5, .085, 'Nine-run entries: primary | Incomplete entries: descriptive only',
                        ha='center', va='bottom', fontsize=8, color='darkred')
            save = figure.savefig
            decorated = False

            def save_partial(path, *save_args, **save_kwargs):
                nonlocal decorated
                if not decorated:
                    decorated = True
                    stem = Path(path).stem
                    figure.set_layout_engine('constrained', rect=(0, .115, 1, .885))
                    if stem == 'metrics':
                        for axis in figure.axes:
                            axis.axvspan(4.5, 6.5, color='orange', alpha=.12, zorder=-1)
                            axis.axvline(4.5, color='darkorange', linestyle=':')
                            metric = axis.get_title()
                            labels = [f"{r['model']}\n{r['completed_runs']}/9 done; "
                                      f"{r[f'{metric}_valid_runs']}/9 defined" for r in summaries]
                            axis.set_xticks(range(len(summaries)), labels, rotation=35, ha='right', fontsize=8)
                        figure.text(.5, .105, 'Reference seed copies are deterministic; orange region is incomplete',
                                    ha='center', fontsize=8)
                    elif stem == 'pr_panels':
                        for axis, row in zip(figure.axes, summaries):
                            if int(row['completed_runs']) < 9:
                                axis.set_facecolor('oldlace')
                                axis.set_title(axis.get_title() + '\nINCOMPLETE: descriptive only')
                    elif stem == 'best_run':
                        figure.set_size_inches(11, 8)
                        axis = figure.axes[0]
                        for line, row in zip(axis.lines, summaries):
                            count = int(row['completed_runs'])
                            line.set_label(line.get_label() + f" | {count}/9 completed"
                                           + ('; incomplete, descriptive only' if count < 9 else '')
                                           + ('; deterministic copies' if row['deterministic_reference'] == 'True' else ''))
                            if count < 9:
                                line.set_linestyle('--')
                        axis.set_title('Post-hoc best TEST run among available completed runs only\n'
                                       'Incomplete entries dashed; entries may select different folds', fontsize=11)
                        axis.legend(fontsize=8)
                return save(path, *save_args, **save_kwargs)

            figure.savefig = save_partial
            return figure, axes

    FunctionType(plots.render.__code__, dict(plots.render.__globals__, plt=DatedPlot(plots.plt)))(snapshot)
    names = [f'{stem}.{extension}' for stem in ('metrics', 'pr_panels', 'best_run')
             for extension in ('png', 'svg')]
    publish.mkdir(parents=True, exist_ok=True)
    for name in names:
        shutil.copyfile(snapshot / name, publish / name)
    publication = {'scope': 'partial_finished_models', 'warning': WARNING, 'use_policy': POLICY,
                   'captured_at_utc': manifest['captured_at_utc'], 'completed_counts': counts,
                   'snapshot': str(snapshot), 'evaluation_sha256': sha256(snapshot / 'evaluation.json'),
                   'scored_artifacts': evidence['artifacts'],
                   'rendering_sources': {p.name: sha256(p) for p in (Path(__file__), Path(plots.__file__),
                                                                  Path(__file__).with_name('partial_benchmark.py'))},
                   'overwrite_policy': 'Final full-terminal delivery overwrites these same six visible paths; '
                                       'this dated partial snapshot remains separate. No canonical seal is created.',
                   'artifacts': {name: sha256(publish / name) for name in names}}
    write_json(publish / 'partial_publication.json', publication)
    write_json(snapshot / 'partial_publication.json', publication)
    return publication


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--publish', type=Path, required=True)
    args = parser.parse_args()
    deliver_partial(args.snapshot, args.publish)
