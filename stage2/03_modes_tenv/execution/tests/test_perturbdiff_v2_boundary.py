"""Boundary evidence via the approved artifact CLI seam, without live processes."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'script/perturbdiff_v2_boundary.py'


class BoundaryEvidence(unittest.TestCase):
    def test_preserves_files_and_debits_all_attempt_wall_time(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / 'formal'
            states = output / 'scheduling'
            states.mkdir(parents=True)
            jobs = [dict(model=model, fold=fold, seed=seed)
                    for model in ('perturbdiff', 'cpa', 'prnet', 'chemcpa', 'cmonge', 'no_change', 'delta_mean')
                    for fold in range(3) for seed in (17, 42, 2026)]
            for job in jobs[:9]:
                key = f"perturbdiff_fold{job['fold']}_seed{job['seed']}"
                (states / (key + '.json')).write_text(json.dumps(dict(status='completed', identity=job,
                    training_seconds=100, attempts=[dict(wall_seconds=120), dict(wall_seconds=30)])))
            binary = output / 'native.pt'
            binary.write_bytes(b'legacy binary')
            queue = root / 'queue.json'
            queue.write_text(json.dumps(dict(jobs=jobs, output_root=str(output))))
            before = {str(p.relative_to(output)): p.read_bytes() for p in output.rglob('*') if p.is_file()}
            result = subprocess.run([sys.executable, str(SCRIPT), 'capture', '--queue', str(queue),
                '--destination', str(root / 'capture')], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            evidence = json.loads((root / 'capture' / 'boundary.json').read_text())
            self.assertEqual(evidence['ledger']['perturbdiff_fold0_seed17']['consumed_seconds'], 150)
            self.assertEqual(evidence['ledger']['perturbdiff_fold0_seed17']['remaining_seconds'], 14250)
            self.assertFalse(evidence['activation_authorized'])
            self.assertEqual(before, {str(p.relative_to(output)): p.read_bytes() for p in output.rglob('*') if p.is_file()})
            self.assertEqual(evidence['other_jobs'], jobs[9:])
