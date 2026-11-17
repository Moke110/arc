"""Approved dose adapter → native prediction/learning boundary."""
import sys
from pathlib import Path
import unittest
import numpy as np
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'output/api'))
from tsg_dose import model
from tenv5k.learning import fixed_rng
SOURCE=Path('/home/chang/arc/pert_model/TranSiGen/source')
CONCAT=['input-numeric','input-embedding','latent-numeric','latent-embedding']

class MethodTests(unittest.TestCase):
    def test_training_only_transform_rejects_invalid_and_roundtrips(self):
        transform=model.fit_transform('input-numeric',[.1,1.,10.])
        self.assertAlmostEqual(transform['mean'],0.,places=6)
        self.assertAlmostEqual(transform['std'],0.8164965809,places=6)
        np.testing.assert_allclose(model.transform_doses([.1,1.,10.],transform),[-1.2247448714,0.,1.2247448714],rtol=1e-6)
        for doses in [[0,1],[-1,1],[np.nan,1],[1,1]]:
            with self.assertRaises(ValueError):model.fit_transform('input-numeric',doses)
        with self.assertRaises(ValueError):model.transform_doses([1],dict(transform,std=0))

    def test_each_concatenation_uses_dose_in_both_predictive_heads(self):
        transform=model.fit_transform('input-numeric',[.1,1.,10.])
        control=torch.linspace(0,1,4*5000).reshape(4,5000)
        morgan=torch.zeros(4,2048);morgan[:,::7]=1
        for method in CONCAT:
            with self.subTest(method=method):
                adapter=model.adapter(SOURCE,method,transform);network=adapter.build(17,'cpu').eval()
                scalar=torch.tensor(model.transform_doses([.1]*4,transform))[:,None]
                low=torch.cat([morgan,scalar],1);high=torch.cat([morgan,-scalar],1)
                with torch.no_grad(),fixed_rng(123,'cpu'):a=network(control,low)
                with torch.no_grad(),fixed_rng(123,'cpu'):b=network(control,high)
                self.assertEqual(a[3].shape,(4,5000))
                self.assertTrue(all(torch.isfinite(value).all() for value in a))
                self.assertFalse(torch.equal(a[4],b[4]))
                self.assertFalse(torch.equal(a[5],b[5]))
                # Prediction depends only on x1/features; targets never enter forward.
                with torch.no_grad(),fixed_rng(123,'cpu'):again=network(control,low)
                self.assertTrue(torch.equal(a[3],again[3]))
                del network,adapter

    def test_adapter_preserves_requested_order_and_learns_with_finite_gradients(self):
        class Data:
            control=np.linspace(0,1,4*5000,dtype=np.float32).reshape(4,5000)
            expression=control+.1
            def fingerprints(self,rows):return np.asarray([[float(i%2)]*2048 for i in rows],dtype=np.float32)
            def doses(self,rows):return np.array([.1,1.,10.,1.],dtype=np.float32)[rows]
        data=Data();rows=np.array([2,0,3,1])
        for method in CONCAT:
            with self.subTest(method=method):
                transform=model.fit_transform(method,data.doses(rows))
                adapter=model.adapter(SOURCE,method,transform);network=adapter.build(17,'cpu')
                features=np.concatenate([data.fingerprints(rows),model.transform_doses(data.doses(rows),transform)[:,None]],1)
                network.eval()
                with torch.no_grad(),fixed_rng(321,'cpu'):
                    expected=network(torch.tensor(data.control[rows]),torch.tensor(features))[3].numpy()
                with fixed_rng(321,'cpu'):actual=adapter.predict(network,data,rows)
                np.testing.assert_array_equal(actual,expected)
                data.expression=np.full_like(data.expression,999)
                with fixed_rng(321,'cpu'):np.testing.assert_array_equal(actual,adapter.predict(network,data,rows))
                data.expression=data.control+.1
                before=next(network.parameters()).detach().clone()
                loss=adapter.train_epoch(network,data,rows,adapter.optimizer(network),np.random.default_rng(17))
                self.assertTrue(np.isfinite(loss))
                self.assertFalse(torch.equal(before,next(network.parameters())))
                self.assertTrue(all(torch.isfinite(p.grad).all() for p in network.parameters() if p.grad is not None))
                del network,adapter

    def test_disabled_conditioning_matches_native_at_equal_state_and_rng(self):
        native=model.adapter(SOURCE).build(17,'cpu').eval()
        base=native.state_dict();control=torch.ones(2,5000);morgan=torch.ones(2,2048)
        with torch.no_grad(),fixed_rng(44,'cpu'):expected=native(control,morgan)
        transform=model.fit_transform('input-numeric',[.1,1.,10.])
        for method in CONCAT:
            with self.subTest(method=method):
                network=model.adapter(SOURCE,method,transform).build(17,'cpu').eval()
                state=network.state_dict()
                for key,value in state.items():
                    if key.startswith('dose_embedding.'):state[key]=torch.zeros_like(value)
                    elif value.shape==base[key].shape:state[key]=base[key]
                    else:
                        expanded=torch.zeros_like(value)
                        expanded[:,:base[key].shape[1]]=base[key]
                        state[key]=expanded
                network.load_state_dict(state)
                with torch.no_grad(),fixed_rng(44,'cpu'):
                    actual=network(control,torch.cat([morgan,torch.zeros(2,1)],1))
                # Wider GEMMs can change float32 accumulation order; this checks
                # disabled-dose equivalence, not cross-architecture bitwise identity.
                for a,b in zip(actual,expected):torch.testing.assert_close(a,b,rtol=1e-5,atol=1e-6)
                del network

    def test_scaling_uses_unstandardized_log10_one_plus_dose(self):
        class Data:
            control=np.ones((3,5000),dtype=np.float32)
            def fingerprints(self,rows):return np.ones((len(rows),2048),dtype=np.float32)
            def doses(self,rows):return np.array([.2,1.,10.],dtype=np.float32)[rows]
        data=Data();rows=np.arange(3)
        self.assertIsNone(model.fit_transform('scaling',[1,1]))
        adapter=model.adapter(SOURCE,'scaling',None);network=adapter.build(17,'cpu').eval()
        # Independent literal factors for the three registered doses.
        scaled=torch.tensor(np.repeat(np.array([[.079181246,.301029996,1.041392685]],dtype=np.float32).T,2048,axis=1))
        with torch.no_grad(),fixed_rng(55,'cpu'):expected=network(torch.tensor(data.control),scaled)[3].numpy()
        with fixed_rng(55,'cpu'):actual=adapter.predict(network,data,rows)
        np.testing.assert_array_equal(actual,expected)
        with self.assertRaises(ValueError):model.adapter(SOURCE,'scaling',{'std':1})
        with self.assertRaises(ValueError):model.fit_transform('scaling',[0,1])

    def test_identity_film_preserves_native_and_learns_dose_modulation(self):
        transform=model.fit_transform('film',[.1,1.,10.])
        adapter=model.adapter(SOURCE,'film',transform);network=adapter.build(17,'cpu').eval()
        native=model.adapter(SOURCE).build(17,'cpu').eval()
        # Equal native state, preserve FiLM's own declared zero-final initialization.
        network.load_state_dict(native.state_dict(),strict=False)
        x=torch.ones(4,5000);morgan=torch.ones(4,2048)
        low=torch.cat([morgan,torch.full((4,1),-1.)],1)
        high=torch.cat([morgan,torch.full((4,1),1.)],1)
        with torch.no_grad(),fixed_rng(71,'cpu'):expected=native(x,morgan)
        for features in [low,high]:
            with torch.no_grad(),fixed_rng(71,'cpu'):actual=network(x,features)
            for a,b in zip(actual,expected):torch.testing.assert_close(a,b,rtol=0,atol=0)
        network.train();optimizer=adapter.optimizer(network);optimizer.zero_grad()
        result=network(x,high);loss=result[3].square().mean()+result[4].square().mean()+result[5].square().mean()
        loss.backward()
        self.assertTrue(all(torch.isfinite(p.grad).all() for p in network.parameters() if p.grad is not None))
        optimizer.step();network.eval()
        with torch.no_grad(),fixed_rng(71,'cpu'):a=network(x,low)
        with torch.no_grad(),fixed_rng(71,'cpu'):b=network(x,high)
        self.assertFalse(torch.equal(a[4],b[4]));self.assertFalse(torch.equal(a[5],b[5]))
        with self.assertRaises(ValueError):model.adapter(SOURCE,'film',dict(transform,std=0))

if __name__=='__main__':unittest.main()
