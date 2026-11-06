import json,os,signal,sys,time
from pathlib import Path
sys.path.insert(0,'/home/chang/arc/stage2/03_modes_tenv/execution/script')
from schedule import process_snapshot,pressure
root=Path('/home/chang/arc/stage2/03_modes_tenv/output/formal')
rec=root/'scheduling'; pid=972642; dispatcher=938762
state=json.loads((rec/'cmonge_fold0_seed42.json').read_text())
last=(rec/'cmonge_fold0_seed42.json').stat().st_mtime
training=state.get('training_seconds',0); inference=state.get('inference_seconds',0)
history=[]; target=Path('/proc')/str(pid); ident=(target/'stat').read_text().split()[21]
def alive():
 try:
  fields=(target/'stat').read_text().split()
  return fields[21]==ident and fields[2]!='Z'
 except OSError:return False
Path('/tmp/tenv_priority_guardian.ready').write_text(json.dumps({'pid':os.getpid(),'worker_pid':pid,'timestamp':time.time(),'training_seconds':training,'inference_seconds':inference})+'\n')
with (rec/'cmonge_fold0_seed42.priority_guardian.jsonl').open('a',buffering=1) as log:
 while alive():
  now=time.time(); delta=max(0,now-last);last=now
  try: progress=json.loads((root/'runs/cmonge/fold_0/seed_42/progress.json').read_text())
  except (OSError,ValueError):progress={}
  if progress.get('phase')=='inference':inference+=delta
  else:training+=delta
  snap=process_snapshot([pid]);snap.update(timestamp=now,guardian_training_seconds=training,guardian_inference_seconds=inference,phase=progress.get('phase'),iteration=progress.get('iteration'))
  history.append(snap);history=history[-2:]
  reason='timeout' if training>=14400 else 'inference_timeout' if inference>=1800 else 'resource_pressure' if pressure(history,[dict(id='worker',owned=True,started=0)])['stop'] else None
  snap['forced_reason']=reason;log.write(json.dumps(snap)+'\n')
  with (rec/'telemetry.jsonl').open('a') as stream:stream.write(json.dumps(snap)+'\n')
  state.update(training_seconds=training,inference_seconds=inference)
  from schedule import save
  save(rec/'cmonge_fold0_seed42.json',state)
  if reason:
   os.killpg(pid,signal.SIGTERM)
   until=time.monotonic()+15
   while alive() and time.monotonic()<until:time.sleep(.2)
   if alive():os.killpg(pid,signal.SIGKILL)
   break
  time.sleep(5)
 log.write(json.dumps({'timestamp':time.time(),'event':'guardian_exit','worker_alive':alive(),'training_seconds':training,'inference_seconds':inference})+'\n')
