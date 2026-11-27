"""Cold-dose acceptance at the public command/artifact seam."""
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

CLI = Path(__file__).resolve().parents[1] / 'workflow.py'

class ColdDoseAcceptance(unittest.TestCase):
    def test_sparse_solver_receives_single_thread_limit(self):
        sys.path.insert(0, str(CLI.parent))
        from cold_dose import _validation
        keys = [('p', 0.2), ('p', 1.0)]
        member = ('baseline', 'cell_viability', 'p')
        with patch('scipy.optimize.milp', return_value=SimpleNamespace(x=[1, 0])) as solver:
            selected = _validation(keys, {k: {member} for k in keys}, {member}, [member[:2]], 1, 0)
        self.assertEqual(selected, {keys[0]})
        self.assertEqual(solver.call_args.kwargs['options'], {'time_limit': 60, 'threads': 1})

    def invoke(self, *args, ok=True):
        result = subprocess.run([sys.executable, str(CLI), *map(str, args)], capture_output=True, text=True)
        self.assertEqual(result.returncode == 0, ok, result.stderr)
        return result

    def fixture(self, path, sparse=False, reverse=False):
        rows = [dict(condition=c, phenotype='cell_viability', canonical_smiles=s,
                     dose_um=d, direction=p%3-1, source_label=f'{p}-{a}', source_asset='fixture',
                     source_cell=f'{p}-{d}-{a}')
                for c in ['baseline','H-R'] for p,s in enumerate(['C','CC','CCC','CCO','CCN','CO'])
                for d in ([0.04] if sparse and p==0 else [0.04,0.2,1,5]) for a in range(2)]
        with path.open('w') as f:
            w=csv.DictWriter(f,fieldnames=rows[0]); w.writeheader(); w.writerows(reversed(rows) if reverse else rows)

    def test_four_global_folds_are_deterministic_with_retained_parents_and_smoke(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            for name,reverse in [('a',False),('b',True)]:
                self.fixture(root/f'{name}.csv',reverse=reverse)
                self.invoke('prepare','--input',root/f'{name}.csv','--output',root/name,'--settings','cold_dose')
            data=json.loads((root/'a/prepared.json').read_text())
            self.assertEqual(data,json.loads((root/'b/prepared.json').read_text()))
            self.assertEqual(len(data['expected_runs']),24)
            tested=[]
            for split,dose in zip(data['splits'],[0.04,0.2,1,5]):
                self.assertEqual(split['held_dose_um'],dose)
                self.assertEqual(split['dose_generalization'],'extrapolation' if dose in [0.04,5] else 'interpolation')
                grouped={}
                selected={role:[] for role in ['train','validation','test']}
                for row,role in zip(data['rows'],split['roles']):
                    key=(row['parent_id'],row['dose_um'])
                    self.assertEqual(grouped.setdefault(key,role),role)
                    self.assertEqual(role=='test',row['dose_um']==dose)
                    selected[role].append(row)
                    if role=='test': tested.append(row['measurement_id'])
                self.assertEqual(len({(r['parent_id'],r['dose_um']) for r in selected['validation']}),4)
                trainparents={r['parent_id'] for r in selected['train']}
                self.assertTrue({r['parent_id'] for r in selected['test']+selected['validation']}<=trainparents)
            self.assertEqual(len(tested),len(set(tested)))
            self.assertEqual(len(tested),len(data['rows']))
            run=root/'run'
            self.invoke('train','--prepared',root/'a/prepared.json','--output',run,'--condition','baseline','--phenotype','cell_viability','--setting','cold_dose','--fold',3,'--seed',17,'--smoke-epochs',1)
            stats=json.loads((run/'run.json').read_text())['dose_statistics']
            train=[math.log10(r['dose_um']) for r,role in zip(data['rows'],data['splits'][3]['roles']) if role=='train' and r['condition']=='baseline']
            self.assertAlmostEqual(stats['mean'],sum(train)/len(train))
            self.assertLess(stats['mean'],0)
            self.invoke('evaluate','--run',run)
            self.invoke('verify','--run',run,'--prepared',root/'a/prepared.json')
            # Rehashing an altered held-dose declaration cannot authorize a run.
            data['splits'][3]['held_dose_um'] = 1
            canonical = json.dumps(data['splits'],sort_keys=True,separators=(',',':'),allow_nan=False).encode()
            data['split_fingerprint'] = hashlib.sha256(canonical).hexdigest()
            (root/'a/prepared.json').write_text(json.dumps(data))
            result=self.invoke('verify','--run',run,'--prepared',root/'a/prepared.json',ok=False)
            self.assertIn('canonical partitions',result.stderr)


    def test_impossible_parent_coverage_has_actionable_failure(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t); self.fixture(root/'s.csv',sparse=True)
            result=self.invoke('prepare','--input',root/'s.csv','--output',root/'out','--settings','cold_dose',ok=False)
            self.assertIn('parent',result.stderr)
            self.assertIn('coverage',result.stderr)
            self.assertFalse((root/'out').exists())
