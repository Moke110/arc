"""Outer-test access rejects forged validation choices before opening data."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'output/api'))
from tsg_dose.inputs import sha256
from tsg_dose.learner import evaluate,manifest,write_json
from tsg_dose.protocol import METHODS,SEEDS


class SelectionGate(unittest.TestCase):
    def origins(self,root):
        inputs=root/'inputs';inputs.mkdir();write_json(inputs/'manifest.json',{'fixture':'identity only'})
        digest=sha256(inputs/'manifest.json');origins={}
        for method in METHODS:
            for seed in SEEDS:
                path=root/method/str(seed);path.mkdir(parents=True)
                write_json(path/'config.json',dict(scenario='cold-drug',fold=0,method=method,seed=seed,smoke_rows=None,execution_binding={}))
                write_json(path/'terminal.json',dict(status='validation-complete',parameter_count=10))
                write_json(path/'validation_metrics.json',dict(pcc=.5 if method=='scaling' else .1))
                write_json(path/'validation_manifest.json',dict(schema='tsg-dose-validation-v1',input_manifest_sha256=digest,parameter_count=10,outputs={name:sha256(path/name) for name in ['config.json','validation_metrics.json']}))
                manifest(path,inputs);origins[str(path.resolve())]=sha256(path/'validation_manifest.json')
        selection=root/'receipt.json'
        write_json(selection,dict(scenario='cold-drug',fold=0,input_manifest_sha256=digest,method='scaling',validation_runs=origins))
        return inputs,selection,root/'no-dose'/str(SEEDS[0])

    def test_forged_winner_is_rejected_before_any_input_or_model_load(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);inputs,selection,candidate=self.origins(root)
            receipt=json.loads(selection.read_text());receipt['method']='film';write_json(selection,receipt)
            # This fixture deliberately has no expression matrices/checkpoint.
            # An invalid choice must fail before accessing either of them.
            with self.assertRaisesRegex(ValueError,'validation winner'):
                evaluate(inputs,candidate,root/'unavailable-source',selection)

    def test_current_checkpoint_must_be_one_of_the_bound_validation_origins(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);inputs,selection,candidate=self.origins(root)
            substitute=root/'substitute';substitute.mkdir()
            for name in ['config.json','terminal.json','validation_metrics.json','validation_manifest.json']:
                (substitute/name).write_bytes((candidate/name).read_bytes())
            manifest(substitute,inputs)
            with self.assertRaisesRegex(ValueError,'evaluated candidate'):
                evaluate(inputs,substitute,root/'unavailable-source',selection)

    def test_parameter_tie_break_is_bound_to_sealed_validation_origin(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);inputs,selection,candidate=self.origins(root)
            write_json(candidate/'terminal.json',dict(status='validation-complete',parameter_count=1))
            manifest(candidate,inputs)
            with self.assertRaisesRegex(ValueError,'parameter count changed'):
                evaluate(inputs,candidate,root/'unavailable-source',selection)


if __name__=='__main__':unittest.main()
