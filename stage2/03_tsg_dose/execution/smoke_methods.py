"""Bounded dose-method lifecycle checks; never counted as formal origins."""
import argparse
import json
from pathlib import Path
import subprocess
import sys


def run(inputs,source,candidates,log,methods):
    execution=Path(__file__).resolve().parent
    log.mkdir(parents=True,exist_ok=True)
    records=[]
    for method in methods:
        candidate=candidates/('smoke-'+method)
        base=[sys.executable,str(execution/'transigen.py'),'--inputs',str(inputs),'--source',str(source),'--candidate',str(candidate),'--device','cuda']
        for phase,arguments in [('fit',['--method',method,'--smoke-rows','8','--max-epochs','2']),('test',['--evaluate-only'])]:
            with (log/(method+'-'+phase+'.txt')).open('w') as stream:
                subprocess.run([str(execution/'capped.sh'),'/usr/bin/time','-v',*base,*arguments],stdout=stream,stderr=subprocess.STDOUT,check=True)
        terminal=json.loads((candidate/'terminal.json').read_text())
        if terminal['status']!='complete' or not terminal.get('smoke'):raise ValueError('bounded smoke failed')
        evidence=log/method;evidence.mkdir()
        for name in ['config.json','terminal.json','manifest.json','validation_manifest.json','environment.json','validation_metrics.json','learning_curve.csv','validation_replay.json','test_replay.json']:
            path=candidate/name
            if path.exists():(evidence/name).write_bytes(path.read_bytes())
        records.append(dict(method=method,candidate=str(candidate),terminal=terminal))
        print(method+' exact fresh-process validation/test replay passed',flush=True)
    (log/'acceptance.json').write_text(json.dumps(dict(smoke_rows=8,max_epochs=2,formal_fits=0,records=records),indent=2,sort_keys=True)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ['inputs','source','candidates','log']:parser.add_argument('--'+arg,type=Path,required=True)
    parser.add_argument('--methods',nargs='+',required=True)
    args=parser.parse_args();run(args.inputs,args.source,args.candidates,args.log,args.methods)
