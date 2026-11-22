"""Future input generation with explicit ADR0077 provenance; old fits stay frozen."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys
STEP = Path(__file__).resolve().parents[2]
ROOT = STEP.parents[1]
sys.path.insert(0, str(ROOT / 'stage2/01_bases/output/api'))
from bases import view, release


def prepare(view_dir, destination):
    view_dir, destination = Path(view_dir), Path(destination)
    # Reject historical/mixed generations before creating or writing a candidate.
    view.verify(view_dir)
    manifest = release.read_json(view_dir / 'manifest.json')
    if destination.exists() and any(destination.iterdir()):
        raise ValueError('future inputs require a new empty candidate')
    spec = importlib.util.spec_from_file_location('_frozen_prepare', STEP / 'execution/prepare.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.prepare(view_dir, destination)
    config = release.read_json(destination / 'config.json')
    config['view_generation'] = {key: manifest[key] for key in
                               ('schema', 'recipe', 'normalization', 'axis')}
    config['view_generation']['preparation_code_sha256'] = release.sha256(Path(__file__))
    config['view_generation']['preserved_prepare_sha256'] = release.sha256(STEP / 'execution/prepare.py')
    config['view_generation']['manifest_sha256'] = release.sha256(view_dir / 'manifest.json')
    config['view_generation']['tenv_sha256'] = release.sha256(view_dir / 'tenv.h5ad')
    config['view_generation']['evidence_scope'] = 'new axis-local inputs; historical full-axis fits are not evidence for this generation'
    release.write_json(destination / 'config.json', config)
    frozen = release.read_json(destination / 'manifest.json')
    frozen['view_generation'] = config['view_generation']
    frozen['outputs']['config.json'] = release.sha256(destination / 'config.json')
    release.write_json(destination / 'manifest.json', frozen)
    return destination

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--view', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.view, args.candidate)
