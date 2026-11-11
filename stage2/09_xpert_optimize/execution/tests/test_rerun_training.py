"""Observable two-phase CLI contract for ADR 0060."""
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parent))
import test_landmark
run = test_landmark.run


class RerunTrainingTest(unittest.TestCase):
    setUp = test_landmark.LandmarkTest.setUp
    prepare = test_landmark.LandmarkTest.prepare
    def test_rerun_retains_one_checkpoint_and_records_exposure(self):
        policy = self.tmp_path / 'rerun-policy.json'
        policy.write_text(json.dumps({
            'kind': 'qualification', 'training_contract': 'step09-round1-rerun-v1',
            'device': 'cpu', 'kernel': 'math', 'dtype': 'float32',
            'microbatch': 4, 'effective_batch': 16,
            'phases': {name: {'min_passes': 1, 'max_passes': passes,
                              'validation_interval': cadence, 'learning_rate': .0001,
                              'max_updates': 2}
                       for name, passes, cadence in [('tahoe', 40, 50), ('tenv', 5, 10)]}}))
        target = self.tmp_path / 'rerun'
        run('script/xpert.py', 'fit', '--predictor', self.release / 'predictor',
            '--policy', policy, '--fold', 0, '--seed', 17, '--output', target)
        report = json.loads((target / 'run.json').read_text())
        self.assertEqual([p.name for p in target.glob('*.pt')], ['best.pt'])
        self.assertEqual(report['status'], 'complete')
        self.assertEqual(report['purpose'], 'smoke')
        for phase in report['phases'].values():
            self.assertEqual(phase['updates'], 2)
            self.assertEqual(phase['training_well_presentations'], 32)
            self.assertEqual([r['update'] for r in phase['history']], [2])
            self.assertIn(phase['stop_reason'], ('qualification_limit', 'maximum_passes'))
            self.assertEqual(phase['history'][0]['validation_wells'], phase['validation_wells'])
        self.assertEqual(report['phases']['tenv']['initial_checkpoint_sha256'],
                         report['phases']['tahoe']['checkpoint_sha256'])
        run('script/xpert.py', 'evaluate-loss', '--predictor', self.release / 'predictor',
            '--checkpoint', target / 'best.pt', '--evaluator', self.release / 'evaluator',
            '--output', self.tmp_path / 'loss.json')
        loss = json.loads((self.tmp_path / 'loss.json').read_text())
        self.assertEqual(loss['gene_count'], 3)
        self.assertGreater(loss['test_wells'], 0)
        complete = target / 'run.json'
        original = complete.read_bytes()
        complete.unlink()
        rejected = run('script/xpert.py', 'predict', '--predictor', self.release / 'predictor',
                       '--checkpoint', target / 'best.pt', '--role', 'test',
                       '--output', self.tmp_path / 'incomplete-pred', check=False)
        self.assertNotEqual(rejected.returncode, 0)
        self.assertIn('incomplete run', rejected.stderr)
        complete.write_bytes(original)
        run('script/xpert.py', 'predict', '--predictor', self.release / 'predictor',
            '--checkpoint', target / 'best.pt', '--role', 'test', '--output', self.tmp_path / 'pred')
        import shutil
        publication = self.tmp_path / 'publication'
        (publication / 'checkpoints').mkdir(parents=True)
        (publication / 'logs').mkdir()
        final = publication / 'checkpoints' / 'XPert_base-fold-0-seed-17.pt'
        shutil.copyfile(target / 'best.pt', final)
        log = publication / 'logs' / f'{final.stem}.json'
        log.write_text(json.dumps({'status': 'complete', 'training_record': report,
                                   'checkpoint_sha256': report['checkpoint_sha256']}))
        run('script/xpert.py', 'predict', '--predictor', self.release / 'predictor',
            '--checkpoint', final, '--role', 'test', '--output', self.tmp_path / 'published-pred')
        self.assertEqual((self.tmp_path / 'pred/delta.npy').read_bytes(),
                         (self.tmp_path / 'published-pred/delta.npy').read_bytes())
        payload = json.loads(log.read_text())
        payload['checkpoint_sha256'] = 'corrupt'
        log.write_text(json.dumps(payload))
        rejected = run('script/xpert.py', 'predict', '--predictor', self.release / 'predictor',
                       '--checkpoint', final, '--role', 'test',
                       '--output', self.tmp_path / 'bad-published-pred', check=False)
        self.assertNotEqual(rejected.returncode, 0)




    def test_real_optimizer_stops_after_four_increases_and_restores_early_best(self):
        import numpy as np
        import torch
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
        from xpert import train_rerun_phase

        class ScalarModel(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.value = torch.nn.Parameter(torch.tensor(1.))
            def forward(self, control):
                values = self.value.expand_as(control)
                return values, values, values

        class Cohort:
            def __init__(self, size, target):
                self.size, self.target = size, target
            def __len__(self):
                return self.size
            def batch(self, rows, device):
                zero = torch.zeros((len(rows), 2), device=device)
                return (zero,), torch.full_like(zero, self.target)

        model = ScalarModel()
        policy = {'device': 'cpu', 'dtype': 'float32', 'kernel': 'math',
                  'effective_batch': 16, 'microbatch': 4,
                  'phases': {'tahoe': {'max_passes': 40, 'validation_interval': 50,
                                        'learning_rate': .0001}}}
        target = self.tmp_path / 'scalar'
        target.mkdir()
        report = train_rerun_phase(model, 'tahoe', Cohort(2000, 0.), Cohort(7, 2.),
                                  target, {'seed': 17}, policy)
        self.assertEqual(report['stop_reason'], 'three_strict_increases')
        self.assertEqual(report['updates'], 200)
        self.assertEqual(report['minimum_updates'], 125)
        self.assertEqual(report['best_update'], 50)
        self.assertEqual(report['best_training_well_presentations'], 800)
        self.assertEqual([r['validation_wells'] for r in report['history']], [7] * 4)
        second = self.tmp_path / 'gated'
        second.mkdir()
        gated = train_rerun_phase(ScalarModel(), 'tahoe', Cohort(8000, 0.), Cohort(7, 2.),
                                  second, {'seed': 17}, policy)
        self.assertEqual(gated['updates'], 500)
        self.assertFalse(gated['history'][3]['stopping_eligible'])
        self.assertEqual(gated['best_update'], 50)
        saved = torch.load(target / 'best.pt', weights_only=True)
        self.assertEqual(model.value.item(), saved['model']['value'].item())
        self.assertEqual(saved['update'], 50)

    def test_phase_stop_replay_resets_equalities_and_counts_all_points(self):
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
        from xpert import stopping_state
        self.assertEqual(stopping_state([1., 2., 3., 4.]), (3, 3))
        self.assertEqual(stopping_state([1., 2., 2., 3.]), (1, 3))
        self.assertEqual(stopping_state([1., 2., 1.]), (0, 0))
        self.assertEqual(stopping_state([1.] + [2.] * 10), (0, 10))
        with self.assertRaises(FloatingPointError):
            stopping_state([1., float('nan')])


if __name__ == '__main__':
    unittest.main()
