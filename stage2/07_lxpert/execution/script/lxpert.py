"""Command-line entry point for Stage 2 Step 07."""
import argparse
import json
from pathlib import Path

from inputs import prepare, verify


REPO = next(path for path in Path(__file__).resolve().parents if (path / ".git").exists())
STEP = REPO / "stage2/07_lxpert"
STEP02 = REPO / "stage2/02_dataset_corpus"
CORPUS = STEP02 / "output/corpus/v2"
RELEASE = STEP / "output/input"


def source_pins():
    records = STEP02 / "execution/build_record"
    tahoe = json.loads((records / "tahoe_corpus_move.json").read_text())
    cm = json.loads((records / "corpus_build.json").read_text())
    return {"tahoe": tahoe["file_sha256"], "tenv": cm["corpora"]["tenv"]["file_sha256"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="command", required=True)
    build = subcommands.add_parser("prepare", help="build the expression input release")
    check = subcommands.add_parser("verify", help="verify an existing input release")
    for command in (build, check):
        command.add_argument("--release", type=Path, default=RELEASE)
    build.add_argument("--corpus-dir", type=Path, default=CORPUS)
    build.add_argument("--seed", type=int, default=20261109, help="development plate split seed")
    build.add_argument("--fixture", action="store_true", help="allow bounded synthetic corpora")
    check.add_argument("--corpus-dir", type=Path, help="also rehash the source files")
    args = parser.parse_args()
    try:
        if args.command == "prepare":
            manifest = prepare(args.corpus_dir, args.release, args.seed, args.fixture, source_pins())
        else:
            manifest = verify(args.release, args.corpus_dir, source_pins())
    except (KeyError, OSError, TypeError, ValueError) as error:
        parser.exit(1, f"{args.command}: {error}\n")
    print(json.dumps({"command": args.command, "release": str(args.release), **manifest},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
