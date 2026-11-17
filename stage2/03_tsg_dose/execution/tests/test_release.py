"""Downstream release boundary: no original matrices or fit directory needed."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd

STEP=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(STEP/'output/api'))

SOURCE=Path('/home/chang/arc/pert_model/TranSiGen/source')


def fixture(root,method='no-dose'):
    import torch
    from tsg_dose import model
    from tsg_dose.learner import runtime_identity,environment,predict
    from tsg_dose.inputs import sha256
    torch.set_num_threads(2)
    root.mkdir(exist_ok=True);(root/'checkpoints').mkdir()
    genes=[str(i) for i in range(5000)]
    pd.DataFrame(dict(gene_index=genes)).to_csv(root/'genes.csv',index=False)
    contract=dict(axis_size=5000,input_manifest_sha256='a'*64,morgan=dict(radius=2,bits=2048,includeChirality=True,binary=True))
    (root/'input_contract.json').write_text(json.dumps(contract))
    (root/'selection_receipt.json').write_text(json.dumps(dict(schema='tsg-dose-selection-v1',scenario='final',fold=0,method=method,input_manifest_sha256='a'*64)))
    controls=np.zeros((2,5000),np.float32);morgan=np.zeros((2,2048),np.float32);morgan[0,0]=1;morgan[1,1]=1
    counts={}
    class Data:
        control=controls
        def fingerprints(self,rows):return morgan[rows]
        def doses(self,rows):return np.array([.2,10.])[rows]
    expected=None
    transform=model.fit_transform(method,[.2,1.,10.])
    for seed in [17,42,2026]:
        native=model.adapter(SOURCE,method,transform);network=native.build(seed,'cpu')
        config=dict(model.configuration(method,transform),scenario='final',fold=0,seed=seed,device='cpu',smoke_rows=None,max_epochs=1000,patience=50,validation_seed=22617,inference_seed=seed+100000,native_source_sha256=model.SOURCE_DIGEST,runtime_identity=runtime_identity(),backend=model.configure_backend('cpu'))
        torch.save(dict(state_dict=network.state_dict(),config=config,input_manifest_sha256='a'*64),root/'checkpoints'/f'{seed}.pt')
        counts[str(seed)]=sum(p.numel() for p in network.parameters())
        if seed==17:expected=predict(native,network,Data(),np.arange(2),config['inference_seed'])
        del network
    manifest=dict(schema='tsg-dose-release-v1',n_formal_fits=147,scenario_counts={'cold-drug':63,'held-dose':63,'final':21},method=method,seeds=[17,42,2026],role='internal-validation-selected restored checkpoints; no new independent test',input_manifest_sha256='a'*64,parameter_counts=counts,environments={str(s):environment('cpu') for s in [17,42,2026]},outputs={str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file()})
    (root/'manifest.json').write_text(json.dumps(manifest))
    return genes,controls,morgan,expected


class ReleaseTests(unittest.TestCase):
    def test_incomplete_formal_matrix_cannot_be_consumed(self):
        from tsg_dose.release import Release
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'manifest.json').write_text(json.dumps(dict(schema='tsg-dose-release-v1',n_formal_fits=146)))
            with self.assertRaisesRegex(ValueError,'147'):
                Release(root,source=root/'missing-source')

    def test_three_checkpoint_release_replays_without_original_training_inputs(self):
        from tsg_dose.release import Release
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'model';genes,controls,morgan,expected=fixture(root)
            release=Release(root,SOURCE)
            result=release.predict(17,controls=controls,morgan=morgan,doses=[.2,10.],gene_ids=genes,profile_ids=['second','first'])
            self.assertEqual(list(result.index),['second','first'])
            self.assertEqual(result.shape,(2,5000))
            self.assertEqual(result.to_numpy().dtype,np.float32)
            np.testing.assert_array_equal(result.to_numpy(),expected)

    def test_saved_dose_transform_replays_and_wrong_inputs_or_payload_are_rejected(self):
        from tsg_dose.release import Release
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'model';genes,controls,morgan,expected=fixture(root,'film')
            release=Release(root,SOURCE)
            args=dict(controls=controls,morgan=morgan,doses=[.2,10.],gene_ids=genes,profile_ids=['a','b'])
            np.testing.assert_array_equal(release.predict(17,**args).to_numpy(),expected)
            for changed in [dict(gene_ids=genes[::-1]),dict(profile_ids=['a','a']),dict(doses=[0,1]),dict(morgan=np.full((2,2048),.5)),dict(controls=np.full((2,5000),np.nan))]:
                with self.assertRaises(ValueError):release.predict(17,**{**args,**changed})
            with (root/'checkpoints/42.pt').open('ab') as stream:stream.write(b'changed')
            with self.assertRaisesRegex(ValueError,'digest'):release.predict(17,**args)

    def test_rehashed_checkpoint_with_wrong_transform_or_parameter_count_is_rejected(self):
        import torch
        from tsg_dose.release import Release
        from tsg_dose.inputs import sha256
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'model';fixture(root,'film')
            manifest=json.loads((root/'manifest.json').read_text());manifest['parameter_counts']['42']=1
            (root/'manifest.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,'parameter count'):Release(root,SOURCE)
            manifest['parameter_counts']['42']=manifest['parameter_counts']['17']
            path=root/'checkpoints/17.pt';checkpoint=torch.load(path,weights_only=False)
            checkpoint['config']['transform']['std']=0;torch.save(checkpoint,path)
            manifest['outputs']['checkpoints/17.pt']=sha256(path)
            (root/'manifest.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,'transform'):Release(root,SOURCE)


if __name__=='__main__':unittest.main()
