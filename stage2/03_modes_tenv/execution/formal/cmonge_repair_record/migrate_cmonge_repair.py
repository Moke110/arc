"""Guarded one-time migration. Requires parent-reviewed adapter SHA and --execute."""
import argparse,copy,hashlib,json,os,shutil,subprocess,sys,time
from pathlib import Path
R=Path('/home/chang/arc');sys.path.insert(0,str(R/'stage2/03_modes_tenv/execution/script'))
from artifacts import run_identity
p=argparse.ArgumentParser();p.add_argument('--approved-adapter-sha256',required=True);p.add_argument('--execute',action='store_true');a=p.parse_args()
f=R/'stage2/03_modes_tenv/execution/formal';o=R/'stage2/03_modes_tenv/output/formal';record=f/'cmonge_repair_record';old=record/'old_attempts'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
def hashes(d):return {str(x.relative_to(d)):sha(x) for x in sorted(d.rglob('*')) if x.is_file()}
assert Path('/proc/1/comm').read_text().strip() in ('systemd','init','wsl-init')
assert sha(R/'stage2/03_modes_tenv/execution/script/cmonge_adapter.py')==a.approved_adapter_sha256
for rel in ('stage2/03_modes_tenv/execution/models/cmonge/formal.json','stage2/03_modes_tenv/execution/script/cmonge_gpu_adapter.py','stage2/03_modes_tenv/execution/script/cmonge_gpu_adapter.identity.json'):
 assert sha(R/rel)==sha(old/rel),rel
assert subprocess.check_output(['git','branch','--show-current'],cwd=R,text=True).strip()=='dev'
for q in Path('/proc').iterdir():
 if not q.name.isdigit():continue
 try:cmd=(q/'cmdline').read_bytes().replace(b'\0',b' ').decode();st=(q/'stat').read_text().split()[2]
 except (OSError,UnicodeDecodeError):continue
 if st!='Z':assert not ('schedule.py run ' in cmd or 'cmonge_gpu_adapter.py --config ' in cmd or 'schedule.py worker ' in cmd),cmd
original=json.loads((old/'stage2/03_modes_tenv/execution/formal/expected_identities.json').read_text())
assert json.loads((f/'expected_identities.json').read_text())==original
assert sha(f/'queue.json')==sha(old/'stage2/03_modes_tenv/execution/formal/queue.json')
assert sha(f/'commands.json')==sha(old/'stage2/03_modes_tenv/execution/formal/commands.json')
queue=json.loads((f/'queue.json').read_text());jobs={(j['model'],j['fold'],j['seed']):j for j in queue['jobs']}
new=copy.deepcopy(original);changed=[]
for x in new:
 key=(x['model'],x['fold'],x['seed']);j=jobs[key]
 current=run_identity(queue['input_root'],*key,j.get('command'))
 if key[0]!='cmonge':assert current==x['identity'],key
 else:
  assert {k for k in current if current[k]!=x['identity'].get(k)}=={'config_sha256'},key
  changed.append({'model':key[0],'fold':key[1],'seed':key[2],'old_identity':x['identity'],'new_identity':current})
  x['identity']=current
assert len(changed)==9 and len(new)==63
states={};archives={}
for seed,charge in((17,6810.677218487021),(42,1467.71)):
 key=f'cmonge_fold0_seed{seed}';sp=o/'scheduling'/f'{key}.json';s=json.loads(sp.read_text())
 assert s['identity']==next(x['identity'] for x in original if x['model']=='cmonge' and x['fold']==0 and x['seed']==seed)
 assert s['status']==('failed' if seed==17 else 'running')
 s=copy.deepcopy(s);before=copy.deepcopy(s)
 assert charge>=before.get('training_seconds',0)
 assert s.get('attempts',[])==before.get('attempts',[])
 if seed==42:s.setdefault('attempts',[]).append({'log':str(o/'scheduling'/f'{key}.attempt0.log'),'reason':'implementation_repair_graceful_interruption','returncode':1,'adapter_returncode':143,'wall_seconds':1467.71,'training_seconds_charged':1467.71})
 s.update(identity=next(x['identity'] for x in new if x['model']=='cmonge' and x['fold']==0 and x['seed']==seed),status='pending',training_seconds=charge,stop_reason='implementation_repair_authorized',repair_revision='cmonge_behavior_gate_v2',repair_old_identity=before['identity'],repair_evidence=str(old.relative_to(R)))
 s.pop('pid',None)
 for name in('oom_retries','environment_retries','batch_scale'):assert s.get(name)==before.get(name)
 states[seed]=(sp,s)
 directory=o/f'runs/cmonge/fold_0/seed_{seed}';archives[seed]=hashes(directory)
result={'timestamp':time.time(),'approved_adapter_sha256':a.approved_adapter_sha256,'changed_identities':changed,'unchanged_identity_count':54,'training_carry_seconds':{'17':6810.677218487021,'42':1467.71},'old_run_full_file_hashes':archives,'partial_checkpoint_reuse':False,'retries_reset':False,'original_queue_and_commands_unchanged':True}
if not a.execute:
 print(json.dumps({'validation':'passed; dry-run only','changed_identities':len(changed),'unchanged':54,'carry':result['training_carry_seconds']}));sys.exit()
archive=o/'implementation_repair_archive/cmonge_behavior_gate_v1';assert not archive.exists();archive.mkdir(parents=True)
for seed in(17,42):
 src=o/f'runs/cmonge/fold_0/seed_{seed}';dst=archive/f'fold_0/seed_{seed}';dst.parent.mkdir(parents=True,exist_ok=True);shutil.copytree(src,dst);assert hashes(dst)==archives[seed]
 for q in src.iterdir():
  if q.is_dir():shutil.rmtree(q)
  else:q.unlink()
 sp,s=states[seed];write(sp,s);write(src/'run.json',{'identity':s['identity'],'status':'pending','stop_reason':s['stop_reason'],'training_seconds':s['training_seconds'],'repair_revision':s['repair_revision'],'repair_evidence':s['repair_evidence']})
write(f/'expected_identities.json',new);write(record/'migration.json',result)
print('MIGRATED ONLY TWO ATTEMPT STATES AND NINE CMONGE EXPECTED IDENTITIES; NO DISPATCH')
