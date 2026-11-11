import tempfile
import unittest
from pathlib import Path
import sys
REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "stage2/06_pdiff_transfer_optimize/execution/script"))

from paired_resources import resources
from artifacts import sha256, write_json


class PairedPublication(unittest.TestCase):
    def test_scheduler_rejects_completed_prediction_for_an_old_plan(self):
        with tempfile.TemporaryDirectory() as temp:
            runtime = Path(temp)
            step = runtime / 'stage2/06_pdiff_transfer_optimize'
            root = step / 'output/paired-transfer'
            cell = root / 'reference/seed-17'
            write_json(cell / 'plan.json', {'version': 'current'})
            prediction = cell / 'Transfer/prediction.json'
            write_json(prediction, {'version': 'current'})
            slot = root / 'schedule/reference-predict'
            scheduler = slot / 'scheduler'
            scheduler.mkdir(parents=True)
            job = 'reference-17-predict'
            queued = {
                'id': job, 'cwd': str(runtime), 'kind': 'gpu', 'memory_gib': 10,
                'disk_gib': 2, 'timeout_seconds': 1800,
                'command': [str(runtime / 'stage2/_envs/perturbdiff/bin/python'),
                            str(step / 'execution/script/paired_reference.py'), 'predict', str(cell),
                            '--runtime', str(runtime), '--seed', '17', '--plan-sha256', '0' * 64],
            }
            write_json(slot / 'queue.json', {'record_dir': 'scheduler', 'jobs': [queued]})
            log = scheduler / f'{job}.attempt1.log'
            log.write_text(f'{prediction} SHA-256 {sha256(prediction)}\n')
            write_json(scheduler / 'status.json', {'jobs': {job: {
                'identity': queued, 'status': 'completed',
                'attempts': [{'log': str(log), 'outcome': 'completed', 'returncode': 0,
                              'started_at': 1, 'finished_at': 2}],
            }}})
            (scheduler / 'telemetry.jsonl').touch()
            with self.assertRaisesRegex(ValueError, f'{job}.*sealed plan'):
                resources(root, runtime)
