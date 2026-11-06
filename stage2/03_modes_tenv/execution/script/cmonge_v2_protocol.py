"""CMonge v2 predeclared stopping policy; no model/data imports."""
import json
import math
import sys

PROTOCOL='cmonge-v2-full-validation-fixed-draws'
EVAL_EVERY=50
MIN_ITERATIONS=500
PATIENCE=10
RELATIVE_MIN_DELTA=.001
CAP_SECONDS=14400.

def remaining_seconds(request):
    if any(k in request for k in ('resume','resume_checkpoint','restart_checkpoint')):
        raise ValueError('v2 restarts from seed; exact resume unsupported')
    value=request.get('prior_training_seconds')
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not 0<=value<CAP_SECONDS:
        raise ValueError('explicit finite cumulative prior_training_seconds required below cap')
    return CAP_SECONDS-value

class Selection:
    def __init__(self):
        self.best_loss=math.inf; self.best_iteration=None
        self.reference=None; self.bad_checks=0; self.last_iteration=0
    def observe(self, iteration, loss, *, final=False):
        if iteration<=self.last_iteration or (iteration%EVAL_EVERY and not final):
            raise ValueError('validation requires increasing cadence or final iteration')
        if not math.isfinite(loss):raise ValueError('nonfinite validation objective')
        self.last_iteration=iteration
        improved=loss<self.best_loss
        if improved:self.best_loss=loss;self.best_iteration=iteration
        if iteration>=MIN_ITERATIONS:
            if self.reference is None:
                self.reference=loss;self.bad_checks=0
            elif self.reference-loss>=RELATIVE_MIN_DELTA*max(abs(self.reference),1e-12):
                self.reference=loss;self.bad_checks=0
            else:self.bad_checks+=1
        return improved,self.bad_checks>=PATIENCE

def main():
    payload=json.load(sys.stdin);remaining=remaining_seconds(payload);selection=Selection();stop=None
    for event in payload['events']:
        _,done=selection.observe(event['iteration'],event['loss'],final=event.get('final',False))
        if done:stop=event['iteration'];break
    print(json.dumps({'stop_iteration':stop,'best_iteration':selection.best_iteration,'remaining_seconds':remaining}))
if __name__=='__main__':main()
