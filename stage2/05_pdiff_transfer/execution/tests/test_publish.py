"""Publication replay with a missing Transfer prediction and undefined direction."""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

STEP = Path(__file__).resolve().parents[1]
CLI = STEP / 'script/publish.py'
METRICS = ('AUPRC_positive@100', 'AUPRC_negative@100', 'Recall_positive@P80',
           'Recall_negative@P80', 'PCC', 'Spearman', 'RMSE', 'R2')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture(root):
    root.mkdir()
    def write(name, data):
        (root / name).write_text(json.dumps(data))

    runs, curves = [], []
    for scope, arms in {
        'warm-start': ('CM', 'Tahoe', 'Transfer', 'Mix', 'no_change', 'delta_mean'),
        'cold-shared': ('CM', 'no_change', 'delta_mean'),
        'cold-cm-drug': ('Tahoe', 'Transfer', 'Mix'),
        'cold-drug': ('Tahoe', 'Transfer', 'Mix'),
    }.items():
        for arm in arms:
            for seed in (None,) if arm in ('no_change', 'delta_mean') else (17, 42, 2026):
                failure = (scope, arm, seed) == ('warm-start', 'Transfer', 42)
                undefined = (scope, arm, seed) == ('warm-start', 'CM', 2026)
                slot = {'scope': scope, 'dataset': 'tenv', 'fold': 0, 'model': arm, 'seed': seed}
                base = 0.75 if arm == 'Transfer' else 0.25
                metrics = {metric: None if failure or (undefined and metric == METRICS[1])
                           else base for metric in METRICS}
                runs.append({**slot, 'status': 'failed' if failure else 'completed',
                             'reason': 'prediction timed out' if failure else '',
                             'group_count': 0 if failure else 3, **metrics,
                             **{f'{metric}_valid_groups': 0 if value is None else 3
                                for metric, value in metrics.items()}})
                curves.append({**slot, 'precision': [
                    [value, value, value] if value is not None else None
                    for value in (metrics[METRICS[0]], metrics[METRICS[1]])]})
    with (root / 'per_run.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=runs[0])
        writer.writeheader()
        writer.writerows(runs)
    write('pr_points.json', {'recall': [0.0, 0.5, 1.0], 'runs': curves})
    write('learning_points.json', {'checkpoints': {}, 'slots': {}})
    write('undefined_reasons.json', [])
    write('scores.json', {'schema': 'step05-scores-v1', 'dose_limitation': 'The 10 µM endpoint is clamped.',
                         'coverage': {}, 'overlap': {}, 'targets': {}, 'inputs': {},
                         'files': {path.name: digest(path) for path in root.iterdir()}})
    return digest(root / 'scores.json')


class PublicationTest(unittest.TestCase):
    def cli(self, *args, success=True):
        result = subprocess.run([sys.executable, str(CLI), *map(str, args)], capture_output=True, text=True)
        self.assertEqual(result.returncode == 0, success, result.stdout + result.stderr)
        return result

    def test_incomplete_seed_pairs_and_offline_replay(self):
        with tempfile.TemporaryDirectory() as scratch:
            root = Path(scratch)
            scores, release = root / 'scores', root / 'release'
            checksum = fixture(scores)
            self.cli('render', scores, release, '--scores-sha256', checksum)
            pairs = json.loads((release / 'paired_ap100.json').read_text())
            warm = pairs['warm-start']
            self.assertEqual(warm['complete_pairs'], 1)
            self.assertIsNone(warm['mean_difference'])
            self.assertEqual(warm['pairs'][0]['difference'], 0.5)
            self.assertIn('prediction timed out', warm['pairs'][1]['reason'])
            self.assertIn('AUPRC_negative@100', warm['pairs'][2]['reason'])
            self.assertEqual(pairs['cold-drug']['mean_difference'], 0.5)
            self.assertEqual(pairs['cold-cm-drug']['complete_pairs'], 3)
            with (release / 'summary.csv').open(newline='') as stream:
                cold_refs = [row for row in csv.DictReader(stream)
                             if row['scope'].startswith('cold') and row['model'] in ('no_change', 'delta_mean')]
            self.assertEqual({(row['scope'], row['model']) for row in cold_refs},
                             {('cold-shared', 'no_change'), ('cold-shared', 'delta_mean')})
            self.assertEqual(len(cold_refs), 2)
            publication = digest(release / 'publication.json')
            self.cli('verify-replay', scores, release, '--publication-sha256', publication)
            step = root / 'step'
            record = step / 'execution/formal_record'
            record.mkdir(parents=True)
            paths = {}
            for source in release.iterdir():
                name = 'results.md' if source.name == 'README.md' else source.name
                destination = (step if source.suffix in ('.png', '.csv') else record) / name
                shutil.copyfile(source, destination)
                paths[source.name] = destination.relative_to(step).as_posix()
            run = {'publication': {'manifest': {'sha256': publication}, 'tracked_artifact_paths': paths}}
            (record / 'run.json').write_text(json.dumps(run))
            self.cli('verify-release', step)
            self.cli('verify-release', step, '--scores', scores)
            (record / 'results.md').write_bytes(b'tampered')
            self.cli('verify-release', step, success=False)
            (release / 'pr_curves.png').write_bytes(b'tampered')
            self.cli('verify-replay', scores, release, '--publication-sha256', publication, success=False)


if __name__ == '__main__':
    unittest.main()
