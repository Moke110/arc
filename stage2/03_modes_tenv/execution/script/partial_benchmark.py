"""Read-only finished-run scoring, explicitly separate from final 63-run evaluation.

Reuse the frozen evaluator function in an isolated globals namespace, substituting
only the manifest gate, source run root and partial rendering callback. This keeps
all group scoring, run aggregation and selection code exactly the formal code.
No module globals are patched and no final seal is created or accepted here.
"""
import argparse
from datetime import datetime, timezone
from pathlib import Path
from types import FunctionType
import numpy as np
import benchmark
from artifacts import (read_json, write_json, sha256, run_directory, validate_run,
                       verify_input_files, SEEDS)

WARNING = 'PARTIAL finished-run snapshot — training ongoing; not final 63-run results'
POLICY = ('Descriptive test scores only. Never use these results for training, checkpoint '
          'selection, hyperparameters, retries, scheduling or entrant changes. '
          'Nine-run completed entries are primary; incomplete entries are descriptive only. '
          'Best-run curves are post-hoc test-selected and may use different folds.')


def snapshot(input_root, source, expected_path, commands):
    verify_input_files(input_root)
    expected = read_json(expected_path)
    keys = [(r['model'], r['fold'], r['seed']) for r in expected]
    models = list(dict.fromkeys(k[0] for k in keys))
    required = [(m, f, s) for m in models for f in range(3) for s in SEEDS]
    if keys != required or len(keys) != 63:
        raise ValueError('expected frozen 63 unique ordered identities')
    entries = []
    for item in expected:
        model, fold, seed = (item[k] for k in ('model', 'fold', 'seed'))
        directory = run_directory(source, model, fold, seed)
        path = directory / 'run.json'
        run = read_json(path) if path.exists() else None
        entry = {'model': model, 'fold': fold, 'seed': seed,
                 'status': run['status'] if run else 'not_dispatched',
                 'command': commands.get(model)}
        if run and run['status'] == 'completed':
            if run.get('identity') != item['identity']:
                raise ValueError('frozen expected identity mismatch')
            validated = validate_run(input_root, source, model, fold, seed, commands.get(model))
            entry.update(identity=validated['identity'], artifacts=validated['artifacts'],
                         run_sha256=sha256(path))
        entries.append(entry)
    if not any(e['status'] == 'completed' for e in entries):
        raise ValueError('no completed predictions to score')
    return {'scope': 'partial_finished_models', 'sealed': False, 'final_barrier_passed': False,
            'captured_at_utc': datetime.now(timezone.utc).isoformat(), 'warning': WARNING,
            'use_policy': POLICY, 'source_root': str(Path(source).resolve()),
            'input_sha256': sha256(Path(input_root) / 'manifest.json'),
            'expected_identities_sha256': sha256(expected_path), 'models': models, 'runs': entries}


def verify_snapshot(input_root, source, manifest):
    if manifest.get('sealed') is not False or manifest.get('scope') != 'partial_finished_models':
        raise ValueError('partial snapshot required; final seal is forbidden')
    if sha256(Path(input_root) / 'manifest.json') != manifest['input_sha256']:
        raise ValueError('partial input digest changed')
    for entry in manifest['runs']:
        if entry['status'] != 'completed':
            continue
        directory = run_directory(source, entry['model'], entry['fold'], entry['seed'])
        if sha256(directory / 'run.json') != entry['run_sha256']:
            raise ValueError('completed snapshot run changed')
        run = validate_run(input_root, source, entry['model'], entry['fold'], entry['seed'], entry['command'])
        if run['identity'] != entry['identity'] or run['artifacts'] != entry['artifacts']:
            raise ValueError('completed snapshot artifacts changed')
    return manifest


class PartialPlot:
    """Local pyplot facade; add durable warning text before frozen renderer saves."""
    def __init__(self, pyplot):
        self.pyplot = pyplot

    def __getattr__(self, name):
        return getattr(self.pyplot, name)

    def subplots(self, *args, **kwargs):
        figure, axes = self.pyplot.subplots(*args, **kwargs)
        figure.text(.5, .005, WARNING, ha='center', va='bottom', color='darkred', fontsize=10,
                    bbox={'facecolor': 'white', 'alpha': .9, 'edgecolor': 'darkred'})
        return figure, axes


def render_partial(output):
    import plots
    namespace = dict(plots.render.__globals__, plt=PartialPlot(plots.plt))
    FunctionType(plots.render.__code__, namespace)(output)


def evaluate_partial(input_root, source, output, expected_path, commands, workers=2):
    input_root, source, output = map(lambda p: Path(p).resolve(), (input_root, source, output))
    if output == source or output.is_relative_to(source) or source.is_relative_to(output):
        raise ValueError('partial output must be separate from canonical source tree')
    if output.exists() and any(output.iterdir()):
        raise ValueError('partial output must be a fresh directory')
    if workers not in (1, 2):
        raise ValueError('workers must be 1 or 2')
    manifest = snapshot(input_root, source, expected_path, commands)
    output.mkdir(parents=True, exist_ok=True)
    # The frozen evaluator expects this basename, but sealed:false guarantees it
    # cannot satisfy benchmark.verify_barrier, even if copied accidentally.
    write_json(output / 'predictions.json', manifest)
    namespace = dict(benchmark.evaluate.__globals__)
    namespace.update(verify_barrier=lambda *_: verify_snapshot(input_root, source, manifest),
                     run_directory=lambda _output, m, f, s: run_directory(source, m, f, s),
                     replay=lambda *_: None)
    evaluator = FunctionType(benchmark.evaluate.__code__, namespace)
    summary = evaluator(input_root, output, workers)
    verify_snapshot(input_root, source, manifest)
    evidence = read_json(output / 'evaluation.json')
    evidence.update(scope=manifest['scope'], sealed=False, final_barrier_passed=False,
                    warning=WARNING, use_policy=POLICY,
                    sources={name: sha256(Path(__file__).with_name(name)) for name in
                             ('partial_benchmark.py', 'benchmark.py', 'artifacts.py', 'metrics.py', 'plots.py')})
    write_json(output / 'evaluation.json', evidence)
    render_partial(output)
    write_json(output / 'partial_publication.json', {
        'scope': manifest['scope'], 'warning': WARNING, 'use_policy': POLICY,
        'evaluation_sha256': sha256(output / 'evaluation.json'),
        'artifacts': {p.name: sha256(p) for p in sorted(output.iterdir()) if p.is_file()}})
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('input', 'source', 'output', 'expected', 'commands'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--workers', type=int, choices=(1, 2), default=2)
    args = parser.parse_args()
    evaluate_partial(args.input, args.source, args.output, args.expected,
                     read_json(args.commands), args.workers)


if __name__ == '__main__':
    main()
