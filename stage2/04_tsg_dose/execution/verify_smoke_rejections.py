"""Check saved-transform and serialization rejection through the public learner."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import tempfile
import numpy as np
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'output/api'))
from tsg_dose.learner import load,serialize_prediction


def verify(inputs,original,source):
    results={}
    with tempfile.TemporaryDirectory(prefix='tsg-dose-transform-rejection-') as tmp:
        candidate=Path(tmp)
        for name in ['config.json','environment.json','registration.json']:
            shutil.copyfile(original/name,candidate/name)
        config=json.loads((candidate/'config.json').read_text())
        config['transform']['mean']+=1
        saved=torch.load(original/'checkpoint.pt',map_location='cpu',weights_only=False)
        saved['config']=config;torch.save(saved,candidate/'checkpoint.pt')
        (candidate/'config.json').write_text(json.dumps(config))
        try:load(inputs,candidate,source)
        except ValueError as error:
            if 'transform' not in str(error):raise
            results['mismatched_persisted_transform']=str(error)
        else:raise AssertionError('mismatched transform accepted')
    for value in [np.inf,np.nan,1e100]:
        try:serialize_prediction(np.array([[value]]))
        except ValueError as error:results['nonfinite_'+str(value)]=str(error)
        else:raise AssertionError('nonfinite serialized prediction accepted')
    return results


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['inputs','candidate','source','output']:parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();result=verify(args.inputs,args.candidate,args.source)
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
