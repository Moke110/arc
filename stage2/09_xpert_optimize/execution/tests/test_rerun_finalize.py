"""Transaction and provenance boundary of Round 1 replacement."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import unittest
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from rerun_finalize import finalize


class FinalizeTest(unittest.TestCase):
    def fixture(self, root):
        step = root / 'stage2/09_xpert_optimize'
        old = step / 'output/round1'
        (old / 'labels/test').mkdir(parents=True)
        (old / 'labels/test/genes.csv').write_text('hgnc_id\nH1\n')
        (old / 'summary.csv').write_text('historical evidence\n')
        record = step / 'execution/round1_record'
        record.mkdir(parents=True)
        (record / 'plan.json').write_text('{"historical":true}\n')
        inventory = root / 'stage2/evaluation/inventory.json'
        inventory.parent.mkdir(parents=True)
        payload = {'runs': [{'prediction': {'path': str((old/'summary.csv').relative_to(root)),
                                            'sha256': 'kept-digest'}}]}
        inventory.write_text(json.dumps(payload))
        retrieval = step / 'output/retrieval/manifest.json'
        retrieval.parent.mkdir()
        retrieval.write_text(json.dumps({'inventory_runs': payload['runs'], 'files': {'scores.csv':'same'}}))
        (retrieval.parent / 'scores.csv').write_bytes(b'immutable scores')
        staging = step / 'output/staging'
        staging.mkdir()
        for name in ('checkpoints', 'logs', 'figure_data'):
            (staging/name).mkdir()
        for name in ('r1_learning_curves.png', 'r1_metrics.png', 'r1_PR_curves.png'):
            (staging/name).write_bytes(b'replayed')
        (staging / 'checkpoints/selected.pt').write_bytes(b'new weights')
        return step, old, inventory, retrieval, staging

    def test_archive_rebinds_only_paths_and_retains_all_old_evidence(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            step, old, inventory, retrieval, staging = self.fixture(root)
            original_inventory, original_manifest = inventory.read_bytes(), retrieval.read_bytes()
            archive = step / 'execution/history/round1-exploratory'
            with patch('rerun_finalize.replay', return_value=True) as replay:
                result = finalize(root, staging, old, archive, inventory, retrieval,
                                  step / 'execution/round1_record')
            self.assertEqual(replay.call_count, 2)
            self.assertEqual((old/'checkpoints/selected.pt').read_bytes(), b'new weights')
            self.assertEqual((archive/'round1/summary.csv').read_text(), 'historical evidence\n')
            self.assertEqual((archive/'original/inventory.json').read_bytes(), original_inventory)
            self.assertEqual((archive/'original/retrieval-manifest.json').read_bytes(), original_manifest)
            self.assertEqual((archive/'round1_record/plan.json').read_text(), '{"historical":true}\n')
            rebound=json.loads(inventory.read_text())['runs'][0]['prediction']
            self.assertEqual(rebound['sha256'], 'kept-digest')
            self.assertIn('execution/history/round1-exploratory/round1/', rebound['path'])
            self.assertEqual((retrieval.parent/'scores.csv').read_bytes(), b'immutable scores')
            self.assertEqual(result['status'], 'complete')
            self.assertTrue((archive/'cleanup-proposal.json').exists())
            self.assertIn('new 8-run checkpoints', json.loads((archive/'cleanup-proposal.json').read_text())['retain'])

    def test_failed_postswap_verification_rolls_back_publication_and_bindings(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            step, old, inventory, retrieval, staging = self.fixture(root)
            original_inventory, original_manifest = inventory.read_bytes(), retrieval.read_bytes()
            archive=step/'execution/history/round1-exploratory'
            with patch('rerun_finalize.replay', side_effect=[True, ValueError('replay changed')]):
                with self.assertRaisesRegex(ValueError, 'replay changed'):
                    finalize(root, staging, old, archive, inventory, retrieval,
                             step/'execution/round1_record')
            self.assertEqual((old/'summary.csv').read_text(), 'historical evidence\n')
            self.assertEqual(inventory.read_bytes(), original_inventory)
            self.assertEqual(retrieval.read_bytes(), original_manifest)
            self.assertTrue((staging/'checkpoints/selected.pt').exists())
