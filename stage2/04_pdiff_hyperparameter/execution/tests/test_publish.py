"""History/checkpoint guards retained after retrieval publication consolidation."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import tempfile
import unittest

STEP = Path(__file__).resolve().parents[2]
SCRIPT = STEP / 'execution/script'
RECORD = STEP / 'execution/formal_record'
sys.path.insert(0, str(SCRIPT))
from checkpoint_release import verify_checkpoint_index, verify_local_release
from training_policy import verify_best_checkpoint


def learning(record, target):
    return subprocess.run([sys.executable, str(SCRIPT / 'learning_curves.py'), '--record', str(record),
                           '--output', str(target / 'learning.png'), '--evidence', str(target / 'learning.json')],
                          capture_output=True, text=True)

class PublicationHistoryTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_replays_original_histories_and_selected_coordinates(self):
        tmp_path = self.root
        for name in ('first', 'second'):
            result = learning(RECORD, tmp_path / name)
            assert result.returncode == 0, result.stderr
        first, second = (tmp_path / name for name in ('first', 'second'))
        assert (first / 'learning.png').read_bytes() == (second / 'learning.png').read_bytes()
        assert (first / 'learning.json').read_bytes() == (second / 'learning.json').read_bytes()
        assert (first / 'learning.png').read_bytes() == (STEP / 'output/learning_curve.png').read_bytes()
        report = json.loads((first / 'learning.json').read_text())
        for dataset in ('tenv', 'tahoe'):
            assert len(report['plotted'][dataset]) == 9
            for run in report['plotted'][dataset]:
                assert run['best']['point'] in run['validation']['points']
                assert run['best']['point'][1] == min(y for _, y in run['validation']['points'])
                assert len(run['training']['points']) == len(run['validation']['points'])

    def test_rejects_missing_original_history_before_rendering(self):
        tmp_path = self.root
        record = shutil.copytree(RECORD, tmp_path / 'record')
        (record / 'tenv/slots/fold-0/seed-17/run.json').unlink()
        target = tmp_path / 'target'
        result = learning(record, target)
        assert result.returncode != 0 and 'run.json' in result.stderr
        assert not (target / 'learning.png').exists()

    def test_rejects_duplicate_checkpoint_identity(self):
        entries = json.loads((RECORD / 'tenv/checkpoints.json').read_text())['checkpoints']
        entries[1] = entries[0]
        with self.assertRaisesRegex(ValueError, 'nine unique checkpoint identities'):
            verify_checkpoint_index(entries, [], {})

    def test_rejects_fractional_history_step(self):
        run = json.loads((RECORD / 'tahoe/slots/fold-0/seed-17/run.json').read_text())
        run['history'][0]['step'] -= .5
        with self.assertRaisesRegex(ValueError, 'validation history'):
            verify_best_checkpoint(run, 'tahoe fixture')

    def test_local_checkpoint_rejects_changed_formal_hard_link(self):
        tmp_path = self.root
        release = tmp_path / 'release'
        slot = release / 'fold-0/seed-17'; slot.mkdir(parents=True)
        index = release / 'checkpoints.json'; index.write_text('{}')
        entry = {'directory': 'fold-0/seed-17', 'identity': {'fold': 0, 'seed': 17}}
        for filename, key in [('run.json', 'run_sha256'), ('replay.json', 'replay_sha256'),
                              ('request.json', 'request_sha256'), ('recipe.json', 'raw_recipe_sha256'),
                              ('best.pt', 'checkpoint_sha256')]:
            content = filename.encode(); (slot / filename).write_bytes(content)
            entry[key] = hashlib.sha256(content).hexdigest()
        formal = tmp_path / 'formal'; fit = formal / 'slots/fold-0-seed-17/fit'; fit.mkdir(parents=True)
        (fit / 'best.pt').hardlink_to(slot / 'best.pt')
        verify_local_release(release, formal, index, [entry])
        (fit / 'best.pt').write_bytes(b'altered')
        with self.assertRaisesRegex(ValueError, 'best.pt.*hash mismatch'):
            verify_local_release(release, formal, index, [entry])
