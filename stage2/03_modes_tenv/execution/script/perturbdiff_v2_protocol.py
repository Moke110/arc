"""Predeclared train/validation-only PerturbDiff v2 stopping protocol."""
import argparse
import json
import math
from pathlib import Path

PROTOCOL='perturbdiff-v2-all-validation-ema'
class EarlyStopping:
    def __init__(self):
        self.reference=None; self.counter=0; self.best=float('inf'); self.best_step=None
        self.last_step=0
    def observe(self,step,loss,max_steps=200000):
        if step<=self.last_step or (step%1000 and step!=max_steps) or not math.isfinite(loss):
            raise ValueError('invalid v2 checkpoint validation event')
        self.last_step=step
        improved=loss<self.best
        if improved: self.best=loss; self.best_step=step
        if step>=20000:
            if self.reference is None or self.reference-loss>=.01:
                self.reference=loss; self.counter=0
            else: self.counter+=1
        return improved,self.counter>=10

def weighted_validation(values,counts):
    if not values or len(values)!=len(counts) or any(type(n) is not int or not 1<=n<=8 for n in counts) or not all(math.isfinite(v) for v in values):
        raise ValueError('invalid bounded validation batches')
    return math.fsum(v*n for v,n in zip(values,counts))/sum(counts)

def remaining_budget(request):
    # Native v1 best snapshots are prediction weights, never restart state.
    if any(request.get(k) for k in ('resume','resume_checkpoint','restart_checkpoint')):
        raise ValueError('exact native resume unsupported; restart from seed only')
    debit=request.get('prior_training_seconds')
    if isinstance(debit,bool) or not isinstance(debit,(int,float)) or not math.isfinite(debit) or not 0<=debit<14400:
        raise ValueError('explicit cumulative prior_training_seconds in [0,14400) required')
    return 14400-debit

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--history')
    mode.add_argument('--request')
    mode.add_argument('--batches')
    args=parser.parse_args()
    if args.batches:
        batches=json.loads(Path(args.batches).read_text())
        print(json.dumps({'validation_loss':weighted_validation([b['loss'] for b in batches],[b['groups'] for b in batches])}))
        raise SystemExit(0)
    if args.request:
        print(json.dumps({'remaining_training_seconds':remaining_budget(json.loads(Path(args.request).read_text()))}))
        raise SystemExit(0)
    policy=EarlyStopping(); stop_step=None
    for event in json.loads(Path(args.history).read_text()):
        _,stop=policy.observe(event['step'],event['validation_loss'])
        if stop: stop_step=event['step']; break
    print(json.dumps({'protocol':PROTOCOL,'stop_step':stop_step,'best_step':policy.best_step,'counter':policy.counter}))
