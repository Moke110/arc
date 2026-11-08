"""Publication CLI tests against the persisted compact formal evidence."""
import json
from pathlib import Path
import shutil
import subprocess
import sys

STEP = Path(__file__).resolve().parents[2]
RECORD = STEP / 'execution/formal_record'
CLI = STEP / 'execution/script/publish.py'


def publish(record, target):
    return subprocess.run(
        [sys.executable, str(CLI), str(record), str(target / 'output'),
         str(target / 'publication.json')],
        capture_output=True, text=True,
    )


def test_replays_original_histories_and_scores(tmp_path):
    for name in ('first', 'second'):
        run = publish(RECORD, tmp_path / name)
        assert run.returncode == 0, run.stderr
    first, second = (tmp_path / name for name in ('first', 'second'))
    assert (first / 'publication.json').read_bytes() == (second / 'publication.json').read_bytes()
    for name in ('learning_curve.png', 'metrics.png'):
        left = (first / 'output' / name).read_bytes()
        assert left.startswith(b'\x89PNG\r\n\x1a\n')
        assert left == (second / 'output' / name).read_bytes()

    report = json.loads((first / 'publication.json').read_text())
    curves = report['learning_curve']
    assert list(curves) == ['tahoe', 'tenv']
    assert [len(panel) for panel in curves.values()] == [9, 9]
    example = curves['tenv'][0]
    assert (example['fold'], example['seed']) == (0, 17)
    assert example['best']['point'] == [6500, 0.8688014149665833]
    for panel in curves.values():
        for series in panel:
            assert [series[k]['color'] for k in ('training', 'validation', 'best')] == [
                '#2166ac', '#24804b', '#bd252d',
            ]
            assert series['best']['point'] in series['validation']['points']

    scores = report['metrics']
    assert len(scores) == 8
    assert all(len(panel) == 6 for panel in scores.values())
    assert [row['color'] for row in scores['RMSE']] == [
        '#2166ac', '#a5aab1', '#636c77', '#24804b', '#a5aab1', '#636c77',
    ]
    assert scores['RMSE'][3]['height'] == 0.13031433815481058
    for panel in ('PCC', 'Spearman'):
        assert [(row['dataset'], row['status'], row['height'], row['valid_groups'])
                for row in scores[panel] if row['model'] == 'no_change'] == [
            ('tahoe', 'undefined', None, 0), ('tenv', 'undefined', None, 0),
        ]
    assert scores['Recall_positive@P80'][1]['height'] == 0.0
    assert scores['R2'][4]['height'] == -0.019548444264354092
    assert len(report['checkpoint_index']) == 9


def test_rejects_missing_original_history_before_rendering(tmp_path):
    record = shutil.copytree(RECORD, tmp_path / 'record')
    (record / 'tenv/slots/fold-0/seed-17/run.json').unlink()
    target = tmp_path / 'target'
    run = publish(record, target)
    assert run.returncode != 0
    assert 'run.json' in run.stderr
    assert not (target / 'output').exists()


def test_rejects_duplicate_checkpoint_identity_before_rendering(tmp_path):
    record = shutil.copytree(RECORD, tmp_path / 'record')
    path = record / 'tenv/checkpoints.json'
    index = json.loads(path.read_text())
    index['checkpoints'][1] = index['checkpoints'][0]
    path.write_text(json.dumps(index))
    target = tmp_path / 'target'
    run = publish(record, target)
    assert run.returncode != 0
    assert 'nine unique checkpoint identities' in run.stderr
    assert not (target / 'output').exists()
