"""Publication CLI tests against the persisted compact formal evidence."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

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
    for name in ('learning_curve.png', 'metrics.png', 'pr_curves.png'):
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
    pr = report['pr_curves']
    assert len(pr['recall']) == 101
    for direction in ('positive', 'negative'):
        assert [(row['dataset'], row['model'], row['color'], row['linestyle'], row['valid_runs'])
                for row in pr[direction]] == [
            ('tahoe', 'perturbdiff', '#2166ac', '-', 9),
            ('tenv', 'perturbdiff', '#2166ac', '--', 9),
            ('tahoe', 'delta_mean', '#606771', '-', 3),
            ('tenv', 'delta_mean', '#606771', '--', 3),
        ]
        assert all(0 <= lo <= mid <= hi <= 1 for row in pr[direction]
                   for lo, mid, hi in zip(row['min'], row['median'], row['max'], strict=True))
    assert report['figures_sha256']['pr_curves.png'] == hashlib.sha256(
        (first / 'output/pr_curves.png').read_bytes()).hexdigest()
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


def test_rejects_pr_source_from_a_different_evaluation(tmp_path):
    record = shutil.copytree(RECORD, tmp_path / 'record')
    path = record / 'tenv/pr_curves.json'
    source = json.loads(path.read_text())
    source['per_run_sha256'] = '0' * 64
    path.write_text(json.dumps(source))
    target = tmp_path / 'target'
    run = publish(record, target)
    assert run.returncode != 0
    assert 'tenv PR source does not match the formal evaluation' in run.stderr
    assert not (target / 'output').exists()


def test_local_checkpoint_rejects_in_place_change_to_the_formal_hard_link(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(CLI.parent))
    from publish import verify_local_release

    release = tmp_path / 'release'
    slot = release / 'fold-0/seed-17'
    slot.mkdir(parents=True)
    index = release / 'checkpoints.json'
    index.write_text('{}')
    entry = {'directory': 'fold-0/seed-17', 'identity': {'fold': 0, 'seed': 17}}
    for filename, key in [('run.json', 'run_sha256'), ('replay.json', 'replay_sha256'),
                          ('request.json', 'request_sha256'), ('recipe.json', 'raw_recipe_sha256'),
                          ('best.pt', 'checkpoint_sha256')]:
        content = filename.encode()
        (slot / filename).write_bytes(content)
        entry[key] = hashlib.sha256(content).hexdigest()
    formal = tmp_path / 'formal'
    fit = formal / 'slots/fold-0-seed-17/fit'
    fit.mkdir(parents=True)
    (fit / 'best.pt').hardlink_to(slot / 'best.pt')
    verify_local_release(release, formal, index, [entry])

    (fit / 'best.pt').write_bytes(b'altered')
    assert (slot / 'best.pt').samefile(fit / 'best.pt')
    with pytest.raises(ValueError, match='best.pt.*hash mismatch'):
        verify_local_release(release, formal, index, [entry])


def test_rejects_fractional_history_step_in_rehashed_formal_record(tmp_path):
    dataset = shutil.copytree(RECORD / 'tahoe', tmp_path / 'record/tahoe')
    slot = dataset / 'slots/fold-0/seed-17'

    def update(path, change):
        value = json.loads(path.read_text())
        change(value)
        data = json.dumps(value).encode()
        path.write_bytes(data)
        return hashlib.sha256(data).hexdigest()

    def fractional_step(run):
        run['history'][0]['step'] -= 0.5

    run_sha = update(slot / 'run.json', fractional_step)
    prediction_sha = update(slot / 'prediction.json', lambda prediction: prediction.update(run_sha256=run_sha))
    outcomes_sha = update(dataset / 'outcomes.json', lambda inventory: next(
        row for row in inventory['outcomes'] if (row['fold'], row['seed']) == (0, 17)
    ).update(prediction_manifest_sha256=prediction_sha))
    update(dataset / 'evaluation.json', lambda report: report.update(outcomes_sha256=outcomes_sha))

    result = publish(dataset.parent, tmp_path / 'target')
    assert result.returncode != 0
    assert 'invalid tahoe/slots/fold-0/seed-17 validation history' in result.stderr
