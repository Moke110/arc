"""Create or fully replay an exact training-successor binding for accepted tenv truth."""
import argparse
import json
from pathlib import Path
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.truth_successor import create_successor, verify_successor, DEFAULT_RSCRIPT


def main():
    root = next(path for path in Path(__file__).resolve().parents if (path / '.git').exists())
    output = root / 'stage2/01_bases/output'
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--truth', type=Path, default=output / 'tenv-truth')
    parser.add_argument('--original-training', type=Path, default=output / 'historical-training-release')
    parser.add_argument('--training', type=Path, default=output / 'training-corpora')
    parser.add_argument('--original-native', type=Path, default=output / 'tenv-raw-counts/tenv.h5ad')
    parser.add_argument('--native', type=Path, default=output / 'tenv-raw-counts/tenv.h5ad')
    parser.add_argument('--accepted-verification', type=Path,
                        default=root / 'stage2/01_bases/log/records/truth_complete_verify.json',
                        help='complete v1 verification report with three representative raw R refits')
    parser.add_argument('--output', type=Path, default=output / 'tenv-truth-successor')
    parser.add_argument('--rscript', type=Path, default=DEFAULT_RSCRIPT)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--record', type=Path)
    args = parser.parse_args()
    if args.verify:
        result = verify_successor(args.output, args.training, rscript=args.rscript)
    else:
        if args.accepted_verification is None:
            parser.error('--accepted-verification is required to publish a successor')
        result = create_successor(args.truth, args.original_training, args.training,
                                  args.original_native, args.native, args.accepted_verification,
                                  args.output, rscript=args.rscript)
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(encoded)
    print(encoded, flush=True)


if __name__ == '__main__':
    main()
