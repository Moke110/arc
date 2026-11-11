"""Atomic publishing checks physical artifacts from both independent scorers."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from stage2.evaluation import publish
from stage2.evaluation.inventory import sha


class PublishingTest(unittest.TestCase):
    def produce(self, command, check):
        root = Path(command[command.index('--destination') + 1])
        root.mkdir()
        (root / 'per_run.csv').write_text('AUPRC\n0.5\n')
        (root / 'PR_curve.png').write_bytes(b'figure')
        (root / 'manifest.json').write_text(json.dumps({
            'files': {'per_run.csv': sha(root / 'per_run.csv')},
            'figures': {'PR_curve.png': sha(root / 'PR_curve.png')}}))

    def test_step06_figures_are_published_directly_in_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            step = Path(temporary)
            with patch.dict(publish.STEPS, {6: step}), patch.object(publish.subprocess, 'run', self.produce):
                publish.publish(6, step / 'inventory.json')
                destination = step / 'output/retrieval'
                self.assertTrue((step / 'output/PR_curve.png').is_file())
                self.assertFalse((destination / 'PR_curve.png').exists())
                publish.verify_existing(destination, step / 'output')

    def test_mutated_first_build_is_rejected_before_replacement(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / 'retrieval'
            def scorer(command, check):
                self.produce(command, check)
                root = Path(command[command.index('--destination') + 1])
                if root.name == 'score':
                    (root / 'per_run.csv').write_text('AUPRC\n0.9\n')
            with patch.object(publish.subprocess, 'run', scorer):
                with self.assertRaisesRegex(ValueError, 'independent replay differs'):
                    publish.publish(9, Path(temporary) / 'inventory.json', destination=destination)
            self.assertFalse(destination.exists())


if __name__ == '__main__':
    unittest.main()
