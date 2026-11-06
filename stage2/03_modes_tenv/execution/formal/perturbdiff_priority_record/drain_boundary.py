#!/usr/bin/env python
"""Finalize the guarded orphan at its natural boundary, then restart reordered queue."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

ROOT = Path('/home/chang/arc')
SCRIPT = ROOT / 'stage2/03_modes_tenv/execution/script'
sys.path.insert(0, str(SCRIPT))
from schedule import process_snapshot, pressure, save, transition, stop_owned, scope_failure_evidence

OUT = ROOT / 'stage2/03_modes_tenv/output/formal'
REC = OUT / 'scheduling'
AUDIT = ROOT / 'stage2/03_modes_tenv/execution/formal/perturbdiff_priority_record'
QUEUE = AUDIT.parent / 'queue.json'
KEY = 'cmonge_fold0_seed42'
PID = 972642
DISPATCHER = 938762
BASE_TRAINING = 1467.71


def proc(pid):
    try:
        # comm can contain spaces; fields after final ')' begin with state.
        fields = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()
        return fields
    except FileNotFoundError:
        return None


def main():
    assert sys.executable == '/home/chang/miniconda3/envs/arc/bin/python'
    assert subprocess.check_output(['git','branch','--show-current'], cwd=ROOT, text=True).strip() == 'dev'
    audit = json.loads((AUDIT / 'order_audit.json').read_text())
    assert hashlib.sha256(QUEUE.read_bytes()).hexdigest() == audit['after_sha256']
    ready = json.loads(Path('/tmp/tenv_priority_guardian.ready').read_text())
    guardian_pid = ready['pid']
    assert guardian_pid == 1000033 and proc(guardian_pid)[19] == '15921976'
    # The original dispatcher disappeared before this controller existed.
    # Keep the established guardian running until worker natural completion.
    assert proc(DISPATCHER) is None
    assert proc(PID)[19] == '15437966'
    assert proc(972669)[19] == '15438461'
    guard_start = proc(guardian_pid)[19]
    save(AUDIT / 'boundary_waiter_ready.json', dict(timestamp=time.time(),pid=os.getpid(),guardian_pid=guardian_pid))
    while True:
        worker = proc(PID)
        guard = proc(guardian_pid)
        alive_worker = bool(worker and worker[19] == '15437966' and worker[0] != 'Z')
        alive_guard = bool(guard and guard[19] == guard_start and guard[0] != 'Z')
        if not alive_worker and not alive_guard:
            break
        if alive_worker and not alive_guard:
            raise RuntimeError('guardian lost while worker live; parent must restore watchdog')
        time.sleep(5)
    lock = (SCRIPT.parent / '.scheduler.lock').open('a')
    fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert proc(DISPATCHER) is None
    # No other queued worker or surviving member of the original scope may overlap.
    for entry in Path('/proc').iterdir():
        if not entry.name.isdigit():
            continue
        fields = proc(int(entry.name))
        if fields and fields[0] != 'Z':
            assert fields[2] != str(PID), 'original process group still alive'
    state_path = REC / (KEY + '.json')
    state = json.loads(state_path.read_text())
    assert state['status'] == 'running' and state['pid'] == PID
    for path in REC.glob('*.json'):
        if path.name.endswith('.job.json') or path == state_path:
            continue
        value = json.loads(path.read_text())
        assert not (isinstance(value,dict) and value.get('status') == 'running'), str(path)
    save(AUDIT / 'takeover_state.json',state)
    last_event = json.loads((REC / (KEY + '.priority_guardian.jsonl')).read_text().splitlines()[-1])
    assert last_event['event'] == 'guardian_exit' and not last_event['worker_alive']
    state['training_seconds'] = last_event['training_seconds']
    state['inference_seconds'] = last_event['inference_seconds']
    events = [json.loads(line) for line in (REC / (KEY + '.priority_guardian.jsonl')).read_text().splitlines()]
    for event in events:
        state['peak_host_rss_gib'] = max(state.get('peak_host_rss_gib',0),event.get('owned_rss_gib',{}).get(str(PID),0))
        state['peak_workstation_gpu_used_gib'] = max(state.get('peak_workstation_gpu_used_gib',0),event.get('gpu_used_gib',0))
    forced_events = [event.get('forced_reason') for event in events]
    forced = next((reason for reason in reversed(forced_events) if reason),None)
    directory = OUT / 'runs/cmonge/fold_0/seed_42'
    job = json.loads((REC / (KEY + '.job.json')).read_text())
    attempt = len(state.get('attempts', []))
    log_path = REC / f'{KEY}.attempt{attempt}.log'
    worker_start = '15437966'
    # time-v's final exit status is the original scope's process outcome.
    text = log_path.read_text(errors='replace')
    matches = re.findall(r'Exit status:\s*(\d+)',text)
    status = int(matches[-1]) if matches else 1
    if status:
        text += scope_failure_evidence(text)
    if (directory / 'adapter.log').exists():
        text += (directory / 'adapter.log').read_text(errors='replace')
    save(AUDIT / 'worker_exit_evidence.json', dict(timestamp=time.time(),time_v_returncode=status,guardian_forced_reason=forced,guardian_exit=last_event))
    reason = forced
    if reason is None:
        if status == 0:
            code = "import sys,json;sys.path.insert(0,sys.argv[1]);from artifacts import validate_run;j=json.load(open(sys.argv[2]));validate_run(j['input_root'],j['output_root'],'cmonge',0,42,j.get('command'))"
            with (AUDIT / 'boundary_validation.log').open('w') as log:
                try:
                    result = subprocess.run([sys.executable,'-c',code,str(SCRIPT),str(REC/(KEY+'.job.json'))],stdout=log,stderr=subprocess.STDOUT,timeout=300)
                    validation = dict(returncode=result.returncode,timed_out=False)
                except subprocess.TimeoutExpired:
                    validation = dict(returncode=None,timed_out=True,timeout_seconds=300)
            save(AUDIT / 'boundary_validation.json',validation)
            reason = 'completed' if validation['returncode'] == 0 else 'invalid_artifacts'
        elif 'TENV_RETRY_ENVIRONMENT' in text:
            reason = 'environment'
        elif ('RESOURCE_EXHAUSTED' in text and 'out of memory' in text.lower()) or any(marker in text for marker in ('CUDA out of memory','Memory cgroup out of memory','oom-kill','killed by the OOM killer')):
            reason = 'oom'
        else:
            reason = 'process_exit_' + str(status)
    if reason == 'completed':
        checkpoint = json.loads((directory / 'checkpoint.json').read_text())
        state['training_seconds'] = BASE_TRAINING + checkpoint['training_seconds']
        state['inference_seconds'] = checkpoint['inference_seconds']
    rss = re.findall(r'Maximum resident set size \(kbytes\):\s*(\d+)',text)
    if rss:
        state['peak_host_rss_gib'] = max(state.get('peak_host_rss_gib',0),max(map(int,rss))/1024**2)
    state = transition(state,reason,0)
    state.setdefault('attempts',[]).append(dict(returncode=status,reason=reason,
        wall_seconds=max(0,time.monotonic()-int(worker_start)/os.sysconf('SC_CLK_TCK')),
        log=str(log_path),boundary_guardian=True))
    save(state_path,state)
    if state['status'] not in ('completed','pending'):
        save(directory / 'run.json', dict(identity=state['identity'],status=state['status'],stop_reason=reason,training_seconds=state['training_seconds']))
    save(AUDIT / 'boundary_result.json',dict(timestamp=time.time(),state=state,worker_finished=True,
         original_scheduler_stopped=True,next_dispatch_model='perturbdiff',test_scoring=False))
    assert hashlib.sha256(QUEUE.read_bytes()).hexdigest() == audit['after_sha256']
    fcntl.flock(lock.fileno(),fcntl.LOCK_UN); lock.close()
    os.execv(sys.executable,[sys.executable,str(SCRIPT/'schedule.py'),'run',str(QUEUE)])


if __name__ == '__main__':
    main()
