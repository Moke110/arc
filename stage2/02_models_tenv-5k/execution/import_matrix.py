"""Explicit, read-only replay and provenance-preserving import of frozen PRnet runs.

Import is a new evidence record, never a cross-plan resume or a provenance rewrite.
Only the approved original 15 cold PRnet entries are eligible.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

FROZEN_RUNTIME = '551a3c9'
FROZEN_COMMIT = '551a3c91c217c83440bae8e5b3b3ca9502a5de10'
ORIGINAL_PLAN_SHA256 = '9b747f83868e305361897dadd76e50afa9f339a87653a3f5a7abd8d98d545858'
# Exact reviewed metadata-only execution plumbing. The PRnet adapter, inputs,
# stopping, random draws and numerical path stay unchanged. No other code delta
# can be authorized by free-form audit text; new deltas require a new review.
APPROVED_CODE_TRANSITIONS = {
    'output/api/tenv5k/learner.py': (
        'a22f62f3163a670800b78f677a360c64b1e62647be29ce54d6fc6aae4ef92bdd',
        '2ba24ff84112b50bce66e0d78d6241b77ee6cab18380273b23d23ec2e5c6070b'),
    'execution/prnet.py': (
        '25489570955d3310176fc9afea0eb7d9678591762e475e6af39ed652d7d332db',
        '64a430a36416ff5951eb8602c4bcabf8bf067bbb3b9b5ae9ee5fee394a1cef4c'),
    'output/api/tenv5k/attempt_identity.py': (
        None, '58d13fa2b0426941bc29c4772c64e36eb9bb965a2bf26d95b31d872a1bf20641'),
}
STEP_RELATIVE = Path('stage2/02_models_tenv-5k')


def read(path):
    return json.loads(Path(path).read_text())


def validate_compatibility(original, target, audit, source_plan_sha256):
    """Require unchanged science except reviewed, hash-bound learner plumbing."""
    if source_plan_sha256 != ORIGINAL_PLAN_SHA256:
        raise ValueError('import original plan is not the approved frozen plan')
    if original['schema'] != 'tenv5k-matrix-v1' or target['schema'] != 'tenv5k-matrix-v2':
        raise ValueError('unsupported import plan schema')
    for name in ('input_manifest_sha256', 'means_manifest_sha256', 'device', 'runs', 'max_epochs', 'patience'):
        if original[name] != target[name]:
            raise ValueError('import shared identity changed: '+name)
    for name, expected in target['shared']['evaluation_code_sha256'].items():
        if original['code_sha256'].get(name) != expected:
            raise ValueError('import scoring semantics changed: '+name)
    method = target['methods']['prnet']
    if method['scientific_identity'] != 'prnet-native-f32-v1':
        raise ValueError('PRnet scientific identity changed')
    differences = {name for name, value in method['code_sha256'].items()
                   if original['code_sha256'].get(name) != value}
    if differences - set(APPROVED_CODE_TRANSITIONS):
        raise ValueError('PRnet implementation/input identity changed')
    for name in differences:
        if (original['code_sha256'].get(name),method['code_sha256'][name]) != APPROVED_CODE_TRANSITIONS[name]:
            raise ValueError('unreviewed metadata dependency transition: '+name)
    if audit.get('source_plan_sha256') != source_plan_sha256 or audit.get('target_method') != method:
        raise ValueError('compatibility audit is not bound to these plans')
    if audit.get('source_runtime') != FROZEN_RUNTIME or not audit.get('reviewed_by'):
        raise ValueError('explicit reviewed compatibility audit required')
    if set(audit.get('changed_code', {})) != differences or any(not audit['changed_code'][n].strip() for n in differences):
        raise ValueError('every changed method dependency requires equivalence rationale')


def verify_runtime(runtime, original, digest):
    """Git object and disk bytes both bind the complete frozen Python runtime."""
    runtime = Path(runtime).resolve()
    commit = subprocess.check_output(['git', '-C', str(runtime), 'rev-parse', 'HEAD'], text=True).strip()
    expected = subprocess.check_output(['git', '-C', str(runtime), 'rev-parse', FROZEN_RUNTIME], text=True).strip()
    if commit != expected or commit != FROZEN_COMMIT:
        raise ValueError('import requires preserved frozen runtime551a3c9')
    for name, value in original['code_sha256'].items():
        path = Path(name)
        if path.is_absolute() or '..' in path.parts:
            raise ValueError('unsafe runtime path')
        if digest(runtime/STEP_RELATIVE/path) != value:
            raise ValueError('frozen runtime bytes changed: '+name)
        blob = subprocess.check_output(['git', '-C', str(runtime), 'show', commit+':'+str(STEP_RELATIVE/path)])
        import hashlib
        if hashlib.sha256(blob).hexdigest() != value:
            raise ValueError('original plan does not match frozen Git runtime')
    return commit


def verify_import(root, run, target, matrix):
    """Revalidate portable imported evidence at resume, evaluation and publication."""
    root = Path(root)
    if run not in matrix.plan() or run['adapter']!='prnet' or run['scenario']!='cold-drug':
        raise ValueError('only registered cold PRnet imports are accepted')
    receipt = read(root/'import-receipt.json')
    if receipt.get('schema') != 'tenv5k-import-v1' or receipt.get('run') != run or receipt.get('origin') != 'imported; not retrained':
        raise ValueError('invalid imported run identity')
    original_path = root/'import-provenance/original-plan.json'
    original = read(original_path)
    if matrix.digest(original_path) != receipt['source_plan_sha256']:
        raise ValueError('original import plan changed')
    validate_compatibility(original, target, read(root/'import-provenance/compatibility-audit.json'), receipt['source_plan_sha256'])
    if receipt['target_method'] != target['methods']['prnet'] or receipt['target_shared'] != target['shared']:
        raise ValueError('import target identity changed')
    if receipt['source_runtime_commit'] != FROZEN_COMMIT:
        raise ValueError('wrong import runtime')
    for name, value in receipt['artifacts_sha256'].items():
        if Path(name).is_absolute() or '..' in Path(name).parts or matrix.digest(root/name) != value:
            raise ValueError('imported provenance/artifact changed: '+name)
    for name, value in original['code_sha256'].items():
        if matrix.digest(root/'import-provenance/runtime'/name) != value:
            raise ValueError('archived frozen runtime changed')
    if receipt['source_manifest_sha256']!=matrix.digest(root/'manifest.json') or receipt['source_config_sha256']!=matrix.digest(root/'config.json') or receipt['source_environment']!=read(root/'environment.json'):
        raise ValueError('original import metadata changed')
    replay = read(root/'import-provenance/fresh-process-replay.json')
    manifest = read(root/'manifest.json')
    if replay.get('status') != 'exact' or replay.get('prediction_sha256') != manifest['outputs']['prediction.npy']:
        raise ValueError('missing exact fresh-process import replay')
    if replay.get('checkpoint_sha256') != matrix.digest(root/'checkpoint.pt') or replay.get('input_manifest_sha256') != target['input_manifest_sha256']:
        raise ValueError('replay checkpoint/input identity changed')
    if replay.get('config_sha256') != matrix.digest(root/'config.json') or replay.get('environment_sha256') != matrix.digest(root/'environment.json'):
        raise ValueError('replay configuration/environment changed')
    matrix.verify_run(root,run,target['input_manifest_sha256'],target['device'])
    return receipt


def import_original(args, run, matrix):
    source_root = Path(args.import_original).resolve()
    runtime = Path(args.import_runtime).resolve()
    source = matrix.run_path(source_root,run)
    target = matrix.frozen_plan(args)
    plan_path = Path(args.candidate)/'plan.json'
    if plan_path.exists() and read(plan_path) != target:
        raise ValueError('import candidate plan identity changed')
    original_path = source_root/'plan.json'
    original = read(original_path)
    audit = read(args.compatibility_audit)
    validate_compatibility(original,target,audit,matrix.digest(original_path))
    commit = verify_runtime(runtime,original,matrix.digest)
    config = matrix.verify_run(source,run,target['input_manifest_sha256'],target['device'])
    # Native source comes from the original plan, never from an alternate current path.
    native_source = Path(original['sources']['prnet'])
    for name, value in config['upstream_files'].items():
        if matrix.digest(native_source/name) != value:
            raise ValueError('original native source changed')
    destination = matrix.run_path(args.candidate,run)
    if destination.exists():
        raise ValueError('import destination already exists')
    source_identity = matrix.tree_identity(source)
    destination.parent.mkdir(parents=True,exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix='.import-',dir=destination.parent))
    try:
        # Snapshot original bytes verbatim. Source/old candidate are read-only.
        shutil.copytree(source,temporary,dirs_exist_ok=True)
        provenance = temporary/'import-provenance'
        provenance.mkdir()
        shutil.copy2(original_path,provenance/'original-plan.json')
        shutil.copy2(args.compatibility_audit,provenance/'compatibility-audit.json')
        for name in original['code_sha256']:
            copied = provenance/'runtime'/name
            copied.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(runtime/STEP_RELATIVE/name,copied)
        command = [str(matrix.step/'execution/capped.sh'),'/usr/bin/time','-v','-o',str(provenance/'replay.time'),args.python,
                   str(matrix.step/'execution/import_replay.py'),'--runtime',str(runtime/STEP_RELATIVE/'output/api'),
                   '--inputs',str(args.inputs),'--run',str(source),'--source',str(native_source)]
        result = subprocess.run(command,text=True,capture_output=True)
        (provenance/'replay.stderr').write_text(result.stderr)
        if result.returncode:
            raise ValueError('frozen import replay failed: '+result.stderr[-4000:])
        replay = json.loads(result.stdout)
        (provenance/'fresh-process-replay.json').write_text(json.dumps(replay,indent=2)+'\n')
        if matrix.tree_identity(source) != source_identity:
            raise ValueError('original source evidence changed during import')
        artifacts = matrix.tree_identity(temporary)
        receipt = dict(schema='tenv5k-import-v1',run=run,origin='imported; not retrained',
            original_run_path=str(source),original_candidate=str(source_root),source_runtime_path=str(runtime),
            source_runtime_commit=commit,source_plan_sha256=matrix.digest(original_path),
            source_manifest_sha256=matrix.digest(source/'manifest.json'),
            source_config_sha256=matrix.digest(source/'config.json'),source_environment=read(source/'environment.json'),
            target_method=target['methods']['prnet'],target_shared=target['shared'],artifacts_sha256=artifacts,
            replay_command=command)
        (temporary/'import-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        verify_import(temporary,run,target,matrix)
        if not plan_path.exists():
            plan_path.parent.mkdir(parents=True,exist_ok=True)
            plan_path.write_text(json.dumps(target,indent=2)+'\n')
        os.rename(temporary,destination)
    finally:
        if temporary.exists():shutil.rmtree(temporary)
    return receipt
