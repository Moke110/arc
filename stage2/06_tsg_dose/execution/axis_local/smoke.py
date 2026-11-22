"""Exercise the dual-source seven-method fit, score and replay interfaces."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

EXECUTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXECUTION.parent / 'output/api'))
from tsg_round1.inputs import SOURCE_POPULATIONS, sha256
from tsg_round1.learner import verify_candidate
from tsg_round1.matrix import METHODS
from tsg_round1.protocol import SCENARIOS


def run(inputs, source, destination, device):
    """Run seven bounded methods in both sources serially and seal their replay proofs."""
    pair = json.loads((inputs / 'manifest.json').read_text())
    if pair.get('schema') != 'tsg-dose-round1-axis-local-pair-v1':
        raise ValueError('expected a paired axis-local input candidate')
    for name in SOURCE_POPULATIONS:
        if pair['sources'].get(name) != sha256(inputs / name / 'manifest.json'):
            raise ValueError(f'pair source digest mismatch: {name}')
    destination.mkdir(parents=True)
    records = []
    for name in SOURCE_POPULATIONS:
        for method_index, method in enumerate(METHODS):
            scenario = SCENARIOS[method_index % len(SCENARIOS)]
            candidate = destination / f'{name}-{scenario}-{method}'
            command = [sys.executable, str(EXECUTION / 'axis_local/transigen.py'),
                       '--inputs', str(inputs / name), '--source', str(source),
                       '--candidate', str(candidate), '--device', device]
            with (destination / f'{candidate.name}.txt').open('w') as log:
                for options in (['--scenario', scenario, '--method', method, '--smoke-rows', '8', '--max-epochs', '2'],
                                ['--evaluate-only']):
                    subprocess.run(command + options, stdout=log, stderr=subprocess.STDOUT, check=True)
            manifest = verify_candidate(candidate)
            terminal = json.loads((candidate / 'terminal.json').read_text())
            config = json.loads((candidate / 'config.json').read_text())
            if (not terminal['smoke'] or terminal['status'] != 'complete'
                    or config['source'] != name or config['scenario'] != scenario
                    or config['method'] != method
                    or manifest['input_manifest_sha256'] != pair['sources'][name]
                    or any(terminal[phase + '_replay']['process_token'] == terminal['training_process_token']
                           for phase in ('validation', 'test'))):
                raise ValueError(f'incomplete or mismatched smoke origin: {candidate.name}')
            records.append(dict(source=name, scenario=scenario, method=method, candidate=candidate.name,
                                configuration=config,
                                fit_sha256=sha256(candidate / 'manifest.json'), terminal=terminal))
            print(f'{name} {scenario} {method}: bounded exact validation and test replay', flush=True)
    (destination / 'acceptance.json').write_text(json.dumps(dict(
        schema='tsg-dose-round1-axis-local-smoke-v2', formal_fits=0,
        pair_sha256=sha256(inputs / 'manifest.json'), records=records), indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ('inputs', 'source', 'destination'):
        parser.add_argument('--' + argument, type=Path, required=True)
    parser.add_argument('--device', default='cpu')
    args = parser.parse_args()
    run(args.inputs, args.source, args.destination, args.device)
