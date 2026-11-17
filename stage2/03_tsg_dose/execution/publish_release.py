"""Full-matrix candidate acceptance and rollback-safe dose-model publication.

This module deliberately imports no model/matrix libraries at import time.
The expensive existing matrix acceptance and downstream replay run only after
all147 formal origins finish, under the repository cap and monitor.
"""
import json
import hashlib
import os
import shutil
import subprocess
import uuid
from pathlib import Path

ACCEPTANCE='acceptance-cold-drug-held-dose-final.json'


def read(path):return json.loads(Path(path).read_text())


def digest(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):result.update(block)
    return result.hexdigest()


def file_snapshot(root,paths):
    root=Path(root)
    return {str(path):digest(root/path) for path in sorted(set(paths))}


def verify_file_snapshot(root,snapshot):
    if file_snapshot(root,snapshot)!=snapshot:raise ValueError('protected product preservation mismatch')


def capture_preservation(root):
    """Hash protected tracked bytes, acceptance manifests and native code only.

    This intentionally avoids rereading52GB of immutable data during training.
    Payload manifests pin their existing digests; the final offline health gate
    verifies large assets after the long trainer has stopped.
    """
    root=Path(root).resolve()
    protected=['stage2/01_bases','stage2/02_models_tenv-5k','stage3','pert_model/TranSiGen']
    raw=subprocess.check_output(['git','-C',str(root),'ls-files','-z','--',*protected])
    tracked=[p for p in raw.decode().split('\0') if p]
    manifests=[]
    for name in protected[:3]:
        manifests.extend(str(p.relative_to(root)) for p in (root/name/'output').rglob('manifest.json'))
    source=root/'pert_model/TranSiGen/source'
    native=[str(p.relative_to(root)) for p in source.rglob('*.py')]
    return dict(schema='tsg-dose-preservation-v1',scope='protected tracked bytes, accepted output manifests and immutable native Python source; large payload verification deferred to final offline health gate',files=file_snapshot(root,tracked+manifests+native),n_tracked=len(tracked),n_acceptance_manifests=len(manifests),n_native_source_files=len(native))


def verify_matrix(matrix):
    matrix=Path(matrix);accepted=read(matrix/ACCEPTANCE)
    if accepted.get('n_formal_fits')!=147:raise ValueError('publication requires complete147-fit acceptance')
    scenarios=['cold-drug','held-dose','final']
    methods=['no-dose','scaling','input-numeric','latent-numeric','input-embedding','latent-embedding','film']
    expected=[dict(scenario=s,fold=f,method=m,seed=t) for s in scenarios for f in range(1 if s=='final' else 3) for m in methods for t in [17,42,2026]]
    key=lambda value:json.dumps(value,sort_keys=True)
    if accepted.get('schema')!='tsg-dose-matrix-acceptance-v1' or accepted.get('scenarios')!=scenarios or sorted(map(key,accepted.get('runs',[])))!=sorted(map(key,expected)):
        raise ValueError('formal matrix unique coverage mismatch')
    if accepted.get('plan_sha256')!=digest(matrix/'plan.json') or read(matrix/'plan.json').get('runs')!=expected:raise ValueError('formal plan coverage/digest mismatch')
    for run in expected:
        path=matrix/'runs'/run['scenario']/run['method']/str(run['fold'])/str(run['seed'])
        terminal=read(path/'terminal.json')
        if terminal.get('status')!=('validation-complete' if run['scenario']=='final' else 'complete'):raise ValueError('incomplete or failed formal origin')
        manifest=read(path/'manifest.json')
        actual={p.name for p in path.iterdir() if p.is_file() and p.name!='manifest.json'}
        if manifest.get('schema')!='tsg-dose-fit-v1' or set(manifest.get('outputs',{}))!=actual:raise ValueError('formal origin manifest coverage mismatch')
        for name,value in manifest['outputs'].items():
            if Path(name).name!=name or digest(path/name)!=value:raise ValueError('formal origin payload digest mismatch')
    for terminal in (matrix/'attempts').rglob('terminal.json'):
        if read(terminal).get('status')=='failed':raise ValueError('failed scientific attempt blocks publication')
    return accepted


def install(candidate,destination,acceptance):
    """Install only an accepted candidate; rollback any failed replacement.

    The acceptance receipt is produced after actual downstream replay and
    isolated CSV-only figure rendering. Model API code is never replaced.
    """
    candidate=Path(candidate);destination=Path(destination)
    if acceptance.get('status')!='passed' or acceptance.get('n_formal_fits')!=147 or acceptance.get('downstream_exact_replay') is not True or acceptance.get('csv_only_figures') is not True or acceptance.get('candidate_manifest_sha256')!=digest(candidate/'manifest.json'):
        raise ValueError('candidate has no complete downstream/figure acceptance')
    if candidate.resolve()==destination.resolve() or destination.is_symlink():raise ValueError('unsafe release installation destination')
    destination.parent.mkdir(parents=True,exist_ok=True)
    backup=destination.with_name('.'+destination.name+'.backup-'+uuid.uuid4().hex)
    had_previous=destination.exists();installed=False
    try:
        if had_previous:os.replace(destination,backup)
        os.replace(candidate,destination);installed=True
        if digest(destination/'manifest.json')!=acceptance['candidate_manifest_sha256']:raise ValueError('installed release digest changed')
    except BaseException:
        if installed and destination.exists():os.replace(destination,candidate)
        if had_previous and backup.exists():os.replace(backup,destination)
        raise
    if backup.exists():shutil.rmtree(backup)
    return dict(destination=str(destination.resolve()),manifest_sha256=acceptance['candidate_manifest_sha256'],replaced_previous=had_previous)
