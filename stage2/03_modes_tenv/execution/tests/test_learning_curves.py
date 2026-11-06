"""Lightweight history mapping contracts; no model/data/evaluator imports."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location('learning_curves', Path(__file__).parents[1] / 'script/learning_curves.py')
curves = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(curves)


class HistoryMappingTests(unittest.TestCase):
    def test_native_keys_and_actual_coordinates(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            for model, (x, train, validation, _, _) in curves.MAPPINGS.items():
                directory = root / 'runs' / model / 'fold_0/seed_17'
                directory.mkdir(parents=True)
                (directory / 'checkpoint.json').write_text(json.dumps({
                    'identity': {'model': model, 'fold': 0, 'seed': 17},
                    'protocol': curves.V2_MODELS.get(model), 'validation_group_count':138,
                    'history': [{x: 98, train: -2, validation: -3, 'validation_conditions_evaluated':138},
                                {x: 203, train: -4, validation: -5, 'validation_conditions_evaluated':138}]}))
            records = curves.collect(root)
            self.assertEqual(len(records), 63)
            available = [r for r in records if r['history_status'] == 'available']
            self.assertEqual(len(available), 5)
            for record in available:
                self.assertEqual(record['train'], [(98., -2.), (203., -4.)])
                self.assertEqual(record['validation'], [(98., -3.), (203., -5.)])
            self.assertEqual(sum(r['history_status'] == 'missing' for r in records), 40)
            self.assertEqual(sum(r['history_status'] == 'not applicable' for r in records), 18)

    def test_v1_is_separate_from_nine_v2_slots_and_terminal_barrier(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            old = root / 'runs/perturbdiff/fold_0/seed_17'
            old.mkdir(parents=True)
            (old / 'run.json').write_text('{"status":"completed"}')
            (old / 'checkpoint.json').write_text(json.dumps({
                'identity': {'model':'perturbdiff','fold':0,'seed':17},
                'history':[{'step':500,'training_loss':2,'validation_loss':3}]}))
            records = curves.collect(root)
            formal = [r for r in records if r['model']=='perturbdiff']
            self.assertEqual(len(formal), 9)
            self.assertTrue(all(not r['train'] for r in formal))
            self.assertEqual(formal[0]['status'], 'legacy_v1_only')
            legacy = curves.collect_legacy(root)
            self.assertEqual(len(legacy), 1)
            self.assertEqual(legacy[0]['validation_population'], 'first 40 groups / 14 of 46 drugs')
            self.assertEqual(legacy[0]['train'], [(500.,2.)])
            archive = root / 'legacy/perturbdiff-v1/runs/fold_0'
            archive.mkdir(parents=True)
            old.rename(archive / 'seed_17')
            old.mkdir()
            (old / 'run.json').write_text('{"status":"completed"}')
            (old / 'checkpoint.json').write_text(json.dumps({
                'identity': {'model':'perturbdiff','fold':0,'seed':17},
                'protocol':'perturbdiff-v2-all-validation-ema', 'validation_group_count':138,
                'history':[{'step':1000,'training_loss':4,'validation_loss':5}]}))
            formal = [r for r in curves.collect(root) if r['model']=='perturbdiff']
            self.assertEqual(formal[0]['train'], [(1000.,4.)])
            self.assertEqual(sum(bool(r['train']) for r in formal), 1)
            self.assertEqual(curves.collect_legacy(root)[0]['train'], [(500.,2.)])

    def test_interrupted_legacy_without_history_cannot_satisfy_v2_slot(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            directory = root / 'runs/perturbdiff/fold_0/seed_2026'
            directory.mkdir(parents=True)
            (directory / 'run.json').write_text('{"status":"timeout"}')
            record = next(r for r in curves.collect(root) if r['model']=='perturbdiff' and r['seed']==2026 and r['fold']==0)
            self.assertEqual(record['status'], 'legacy_v1_only')
            self.assertEqual(curves.collect_legacy(root)[0]['history_status'], 'missing')

    def test_cmonge_durable_v2_and_archived_v1_are_not_spliced(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            live = root / 'runs/cmonge/fold_0/seed_17'
            old = root / 'legacy/cmonge-v1/runs/fold_0/seed_17'
            live.mkdir(parents=True)
            old.mkdir(parents=True)
            identity = dict(model='cmonge', fold=0, seed=17)
            (old / 'checkpoint.json').write_text(json.dumps(dict(identity=identity,
                history=[dict(iteration=50, training_loss=9, validation_loss=8)])))
            event = dict(iteration=50, training_loss=4, validation_loss=3, validation_conditions_evaluated=138)
            data = dict(identity=identity, protocol=curves.V2_MODELS['cmonge'], history=[event])
            (live / 'history.json').write_text(json.dumps(data))
            formal = [r for r in curves.collect(root) if r['model']=='cmonge']
            self.assertEqual(len(formal), 9)
            self.assertEqual(formal[0]['train'], [(50., 4.)])
            self.assertEqual(formal[0]['status'], 'not_dispatched')
            self.assertTrue(formal[0]['source'].endswith('history.json'))
            legacy = [r for r in curves.collect_legacy(root) if r['model']=='cmonge']
            self.assertEqual(legacy[0]['train'], [(50., 9.)])
            event['validation_conditions_evaluated'] = 1
            (live / 'history.json').write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, 'population mismatch'):
                curves.collect(root)

    def test_cmonge_live_legacy_completion_is_not_formal_v2(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            live = root / 'runs/cmonge/fold_0/seed_17'
            live.mkdir(parents=True)
            (live / 'run.json').write_text('{"status":"completed"}')
            (live / 'checkpoint.json').write_text(json.dumps(dict(
                identity=dict(model='cmonge',fold=0,seed=17),
                history=[dict(iteration=50,training_loss=9,validation_loss=8)])))
            formal = [r for r in curves.collect(root) if r['model']=='cmonge']
            self.assertEqual(len(formal), 9)
            self.assertTrue(all(not r['train'] for r in formal))
            self.assertEqual(formal[0]['status'], 'legacy_v1_only')
            self.assertEqual(len([r for r in curves.collect_legacy(root) if r['model']=='cmonge']), 1)

    def test_no_fabricated_coordinates_or_nonfinite_values(self):
        for history in ([{'step': 3, 'loss': 1}, {'step': 2, 'loss': 1}],
                        [{'step': 3, 'loss': float('nan')}]):
            with self.assertRaises(ValueError):
                curves.points(history, 'step', 'loss')
        self.assertEqual(curves.points([], 'epoch', 'loss'), [])

    def test_checkpoint_identity_is_required(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            directory = root / 'runs/prnet/fold_0/seed_17'
            directory.mkdir(parents=True)
            (directory / 'checkpoint.json').write_text('{"identity": {"model": "cpa"}}')
            with self.assertRaisesRegex(ValueError, 'identity mismatch'):
                curves.collect(root)


if __name__ == '__main__':
    unittest.main()
