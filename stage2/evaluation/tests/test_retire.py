"""Retirement must preserve canonical figures at the Step 06 output root."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from stage2.evaluation import retire


class RetirementTest(unittest.TestCase):
    def test_manifest_owned_step06_figures_survive_cleanup(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / '06_example/output'
            release = output / 'retrieval'
            release.mkdir(parents=True)
            for name in ('PR_curve.png', 'PR_curve.svg', 'AUPRC.png', 'AUPRC.svg', 'old_PR.png'):
                (output / name).write_bytes(b'figure')
            (release / 'manifest.json').write_text(json.dumps({'files': {}, 'figures': {
                name: 'bound-hash' for name in ('PR_curve.png', 'PR_curve.svg', 'AUPRC.png', 'AUPRC.svg')}}))
            with patch.object(retire, 'ROOT', root):
                self.assertEqual(retire.candidates(), [output / 'old_PR.png'])


if __name__ == '__main__':
    unittest.main()
