"""Public input/fit boundaries reject incomplete identities and false replay."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'output/api'))
from tsg_dose.inputs import Inputs,REQUIRED
from tsg_dose.protocol import INPUT_SCHEMA
from tsg_dose.learner import accept_replay,train

class IntegrityTests(unittest.TestCase):
    def test_freeze_rejects_missing_and_changed_payloads_before_opening_matrices(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'manifest.json').write_text(json.dumps(dict(schema=INPUT_SCHEMA,outputs={})))
            with self.assertRaisesRegex(ValueError,'coverage'):Inputs(root)
            (root/'manifest.json').write_text(json.dumps(dict(schema=INPUT_SCHEMA,outputs={name:'0'*64 for name in REQUIRED})))
            for name in REQUIRED:(root/name).write_bytes(b'corrupted')
            with self.assertRaisesRegex(ValueError,'digest'):Inputs(root)

    def test_fit_requires_attempt_registration_and_replay_requires_new_process(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            with self.assertRaisesRegex(ValueError,'preregistered execution'):
                train(root/'missing-inputs',root/'formal',root/'missing-source')
            (root/'terminal.json').write_text(json.dumps(dict(status='validation-trained-awaiting-replay',training_pid=os.getpid())))
            with self.assertRaisesRegex(ValueError,'fresh process'):
                accept_replay(root/'missing-inputs',root,root/'missing-source')

if __name__=='__main__':unittest.main()
