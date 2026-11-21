"""Fit and verify a bounded fresh no-dose run in both tenv scenarios."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'output/api'))
from tsg_kpgt.inputs import sha256, WIDTH
from tsg_kpgt.learner import verify_candidate
from tsg_round1.protocol import SCENARIOS


def run(inputs, source, destination, device):
    destination.mkdir(parents=True)
    input_digest = sha256(inputs / 'manifest.json')
    records = []
    for scenario in SCENARIOS:
        candidate = destination / scenario
        command = [
            sys.executable,
            str(HERE / 'transigen.py'),
            '--inputs',
            str(inputs),
            '--source',
            str(source),
            '--candidate',
            str(candidate),
            '--device',
            device,
        ]
        with (destination / f'{scenario}.txt').open('w') as log:
            for args in (
                ['--scenario', scenario, '--smoke-rows', '8', '--max-epochs', '2'],
                ['--evaluate-only'],
            ):
                subprocess.run(
                    command + args, stdout=log, stderr=subprocess.STDOUT, check=True
                )
        manifest = verify_candidate(candidate)
        terminal = json.loads((candidate / 'terminal.json').read_text())
        config = json.loads((candidate / 'config.json').read_text())
        tokens = {
            terminal[key]
            for key in ('training_process_token', 'test_evaluation_process_token')
        }
        tokens.update(
            terminal[phase + '_replay']['process_token']
            for phase in ('validation', 'test')
        )
        if (
            not terminal['smoke']
            or terminal['status'] != 'complete'
            or not terminal['weights_updated']
            or config['source'] != 'tenv'
            or config['scenario'] != scenario
            or config['method'] != 'no-dose'
            or config['representation_width'] != WIDTH
            or manifest['input_manifest_sha256'] != input_digest
            or len(tokens) != 4
        ):
            raise ValueError(f'invalid KPGT smoke evidence: {scenario}')
        records.append(
            dict(
                scenario=scenario,
                configuration=config,
                fit_sha256=sha256(candidate / 'manifest.json'),
                terminal=terminal,
            )
        )
        print(f'{scenario}: exact independent validation and test replay', flush=True)
    (destination / 'acceptance.json').write_text(
        json.dumps(
            dict(
                schema='tsg-dose-tenv-kpgt-smoke-v1',
                input_manifest_sha256=input_digest,
                formal_fits=0,
                records=records,
            ),
            sort_keys=True,
            indent=2,
        )
        + '\n'
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ('inputs', 'source', 'destination'):
        parser.add_argument('--' + argument, type=Path, required=True)
    parser.add_argument('--device', default='cpu')
    args = parser.parse_args()
    run(args.inputs, args.source, args.destination, args.device)
