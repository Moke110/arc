"""Build or independently replay evaluator-only tenv parent-dose truth."""
import argparse
import json
from pathlib import Path
from bases.truth_build import build, verify, smoke


def main():
    root = next(path for path in Path(__file__).resolve().parents if (path / '.git').exists())
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--training', type=Path, default=root / 'stage2/01_bases/output/training-v1')
    parser.add_argument('--native', type=Path, required=True,
                        help='explicit bases-native raw integer count corpus; never normalized training expression')
    parser.add_argument('--output', type=Path, default=root / 'stage2/01_bases/output/truth-tenv-v1')
    parser.add_argument('--rscript', type=Path, default=Path('/home/chang/miniconda3/envs/biomni_e1/bin/Rscript'))
    parser.add_argument('--batch-size', type=int, default=25)
    parser.add_argument('--workers', type=int, choices=(1, 2), default=2)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--smoke', action='store_true', help='actual ordinary/multi-plate R preflight only')
    parser.add_argument('--replay-groups', type=int, default=3)
    parser.add_argument('--record', type=Path)
    args = parser.parse_args()
    if args.smoke:
        result = smoke(args.training, args.native, args.output, args.rscript)
    elif args.verify:
        result = verify(args.training, args.native, args.output, args.rscript,
                        replay_groups=args.replay_groups)
    else:
        result = build(args.training, args.native, args.output, args.rscript,
                       batch_size=args.batch_size, workers=args.workers)
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(encoded)
    print(encoded, flush=True)


if __name__ == '__main__':
    main()
