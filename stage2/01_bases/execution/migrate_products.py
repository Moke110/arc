"""One-shot, idempotent rename migration to the ADR 0062 stable product paths (#156).

The migration is data, not a rebuild. Every pair is proved before anything is
touched: the source must exist, the destination must be absent, and both must
live on the same filesystem so the move is a rename. Each moved tree keeps its
inode, size and mtime, and every non-expression file keeps its sha256, so the
record proves bytes were not rewritten. If any step fails, all completed renames
are reversed and the previously accepted products stay where they were.

Run it from the checkout that actually owns the large payloads (the main
checkout has the 34 GB ignored output tree; a linked worktree does not)::

    python3 stage2/01_bases/execution/migrate_products.py \
        --step-root /home/chang/arc/stage2/01_bases \
        --record /path/to/log/156-migration-record.json

Exit status: 0 migrated or already migrated, 1 aborted/failed (products intact),
2 unusable step root or record. Reverse every recorded rename with::

    python3 stage2/01_bases/execution/migrate_products.py --rollback RECORD
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time

SCHEMA = "bases-156-migration-record-v1"
EXPRESSION_SUFFIXES = {".h5ad", ".npy", ".npz", ".bin"}

# The single source of truth for this ticket's layout change. Paths are relative
# to the step root. These are pre-migration inputs, never runtime defaults: the
# code in output/api and execution only names the stable paths.
MOVES = (
    # id, source, destination, kind
    ("gene-identity", "output/identity-v1", "output/gene-identity", "dir"),
    ("frozen-tables", "output/native-v1/tables", "output/gene-identity/tables", "dir"),
    ("tenv-raw-h5ad", "output/native-v1/tenv.h5ad", "output/tenv-raw-counts/tenv.h5ad", "file"),
    ("tenv-raw-manifest", "output/native-v1/tenv_manifest.json", "output/tenv-raw-counts/tenv_manifest.json", "file"),
    ("tenv-raw-probes", "output/native-v1/tenv_probe_mapping.csv", "output/tenv-raw-counts/tenv_probe_mapping.csv", "file"),
    # Resumption state leaves output/ before the release directory itself moves.
    ("tahoe-resumption-state", "output/tahoe-counts-v1/.checkpoints", "execution/checkpoints/tahoe-raw-counts", "dir"),
    ("tahoe-raw-counts", "output/tahoe-counts-v1", "output/tahoe-raw-counts", "dir"),
    ("shared-expression", "output/shared-v2", "output/shared-expression", "dir"),
    ("training-corpora", "output/training-v2", "output/training-corpora", "dir"),
    ("tenv-truth", "output/truth-tenv-v1", "output/tenv-truth", "dir"),
    ("tenv-truth-successor", "output/truth-tenv-v2", "output/tenv-truth-successor", "dir"),
    ("tenv-hvg", "output/hvg-tenv-v1", "output/tenv-hvg", "dir"),
    ("historical-training-release", "output/training-v1", "output/historical-training-release", "dir"),
    ("truth-fit-resumption-state", "output/truth-tenv-v1.building", "execution/checkpoints/truth-tenv-building", "dir"),
    ("shared-v1-mappings", "output/shared-v1", "log/historical/shared-v1-mappings", "dir"),
    ("diagnostic-evidence", "output/hvg-metrics-diagnostic-v1", "log/diagnostic", "dir"),
    ("repo-health", "execution/build_record/final_health", "log/health/final_health", "dir"),
    ("execution-records", "execution/build_record", "log/records", "dir"),
    ("figure-png", "diagnostic_tenv_pr.png", "figure/diagnostic_tenv_pr.png", "file"),
    ("figure-csv", "diagnostic_tenv_pr.csv", "figure/data/diagnostic_tenv_pr.csv", "file"),
)

EMPTY_SOURCE_DIRS = ("output/native-v1",)


def digest(path):
    value = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            value.update(chunk)
    return value.hexdigest()


def evidence(root):
    """Record every file's identity; non-expression bytes are also hashed."""
    root = Path(root)
    items = []
    for base, directories, files in os.walk(root):
        directories.sort()
        for name in sorted(files):
            path = Path(base) / name
            stat = path.stat()
            relative = str(path.relative_to(root))
            item = {"path": relative, "size": stat.st_size, "device": stat.st_dev,
                    "inode": stat.st_ino, "mtime_ns": stat.st_mtime_ns,
                    "sha256": None}
            if Path(name).suffix not in EXPRESSION_SUFFIXES and path.is_file():
                item["sha256"] = digest(path)
            items.append(item)
        for name in sorted(directories):
            path = Path(base) / name
            if path.is_symlink():
                raise ValueError(f"unsupported symlink in migrated tree: {path}")
    return items


def tree_evidence(root, kind):
    root = Path(root)
    if kind == "file":
        stat = root.stat()
        return [{"path": root.name, "size": stat.st_size, "device": stat.st_dev,
                 "inode": stat.st_ino, "mtime_ns": stat.st_mtime_ns,
                 "sha256": digest(root) if root.suffix not in EXPRESSION_SUFFIXES else None}]
    return evidence(root)


def compare(before, after):
    differences = []
    if [item["path"] for item in before] != [item["path"] for item in after]:
        differences.append("file set changed")
    for left, right in zip(before, after):
        for field in ("path", "size", "device", "inode", "mtime_ns", "sha256"):
            if left[field] != right[field]:
                differences.append(f"{left['path']}: {field} {left[field]!r} -> {right[field]!r}")
    return differences


def plan(step):
    entries = []
    for name, source, destination, kind in MOVES:
        old, new = step / source, step / destination
        if old.exists() and not new.exists():
            state = "pending"
        elif not old.exists() and new.exists():
            state = "migrated"
        elif old.exists() and new.exists():
            state = "conflict"
        else:
            state = "absent"
        entries.append({"id": name, "source": source, "destination": destination,
                        "kind": kind, "old": str(old), "new": str(new), "state": state})
    return entries


def nearest_existing(path):
    path = Path(path)
    while not path.exists():
        if path.parent == path:
            raise ValueError(f"no existing ancestor for {path}")
        path = path.parent
    return path


def validate(step, entries):
    problems = []
    for entry in entries:
        if entry["state"] != "pending":
            continue
        old, new = Path(entry["old"]), Path(entry["new"])
        if old.stat().st_dev != nearest_existing(new.parent).stat().st_dev:
            problems.append(f"{entry['id']}: source and destination are on different filesystems")
        if entry["kind"] == "file" and not old.is_file():
            problems.append(f"{entry['id']}: source is not a regular file")
        if entry["kind"] == "dir" and not old.is_dir():
            problems.append(f"{entry['id']}: source is not a directory")
    return problems


def rename(step, entry, *, dry_run):
    old, new = Path(entry["old"]), Path(entry["new"])
    new.parent.mkdir(parents=True, exist_ok=True)
    if dry_run:
        return None
    before = tree_evidence(old, entry["kind"])
    os.rename(old, new)
    after = tree_evidence(new, entry["kind"])
    return before, after


def rollback_moves(step, entries):
    reversed_entries = []
    for entry in reversed(entries):
        old, new = Path(entry["old"]), Path(entry["new"])
        if new.exists() and not old.exists():
            old.parent.mkdir(parents=True, exist_ok=True)
            os.rename(new, old)
            reversed_entries.append(entry["id"])
    return reversed_entries


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--step-root", type=Path,
                        default=Path(__file__).resolve().parents[1])
    parser.add_argument("--record", type=Path, help="write the migration record here")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--rollback", type=Path, metavar="RECORD",
                        help="reverse the renames a migration record lists")
    args = parser.parse_args()
    step = args.step_root.resolve()
    if not (step / "output").is_dir() or not (step / "execution").is_dir():
        print(f"ERROR: {step} is not a bases step root", file=sys.stderr)
        return 2

    if args.rollback:
        record = json.loads(args.rollback.read_text())
        moves = [{"id": entry["id"], "old": entry["old"], "new": entry["new"]}
                 for entry in record["moves"]]
        reversed_entries = rollback_moves(step, moves)
        print(json.dumps({"schema": SCHEMA, "result": "ROLLED_BACK",
                          "reversed": reversed_entries, "record": str(args.rollback)},
                         indent=2, sort_keys=True))
        return 0

    entries = plan(step)
    problems = validate(step, entries)
    blocked = [entry for entry in entries if entry["state"] in ("conflict", "absent")]
    for entry in blocked:
        problems.append(f"{entry['id']}: state={entry['state']} "
                        f"({entry['old']} -> {entry['new']})")
    if problems:
        print(json.dumps({"schema": SCHEMA, "result": "ABORTED", "problems": problems},
                         indent=2, sort_keys=True))
        print("No product was moved; the accepted layout is unchanged.", file=sys.stderr)
        return 1
    pending = [entry for entry in entries if entry["state"] == "pending"]
    if not pending:
        print(json.dumps({"schema": SCHEMA, "result": "ALREADY_MIGRATED",
                          "moves": len(entries)}, indent=2, sort_keys=True))
        return 0

    record = {"schema": SCHEMA, "ticket": 156, "parent": 154,
              "step_root": str(step), "started": time.time(),
              "command": [sys.executable, *sys.argv],
              "rollback": "python3 execution/migrate_products.py --rollback <record>",
              "method": "os.rename on one filesystem; non-expression sha256 and inode/size/mtime preserved",
              "moves": [], "cleanup": [], "problems": []}
    if args.dry_run:
        record["result"] = "DRY_RUN"
        record["moves"] = [{key: entry[key] for key in ("id", "source", "destination", "kind")}
                           for entry in pending]
        print(json.dumps(record, indent=2, sort_keys=True))
        return 0

    completed = []
    try:
        for entry in pending:
            before, after = rename(step, entry, dry_run=False)
            differences = compare(before, after)
            if differences:
                raise RuntimeError(f"{entry['id']}: post-move evidence differs: {differences[:5]}")
            completed.append(entry)
            record["moves"].append({"id": entry["id"], "source": entry["source"],
                                    "destination": entry["destination"], "kind": entry["kind"],
                                    "old": entry["old"], "new": entry["new"],
                                    "files": after, "verified": True,
                                    "state": "moved"})
            print(f"moved {entry['source']} -> {entry['destination']} "
                  f"({len(after)} file(s), sha256 preserved)", flush=True)
        for relative in EMPTY_SOURCE_DIRS:
            directory = step / relative
            if directory.is_dir() and not any(directory.iterdir()):
                directory.rmdir()
                record["cleanup"].append({"action": "rmdir", "path": str(directory),
                                          "reason": "source directory is empty after renames"})
                print(f"removed empty source directory {relative}", flush=True)
    except Exception as error:  # noqa: BLE001 - any failure must reverse the renames
        record["problems"].append(str(error))
        reversed_entries = rollback_moves(step, completed)
        record["result"] = "FAILED"
        record["reversed"] = reversed_entries
        record["finished"] = time.time()
        if args.record:
            args.record.parent.mkdir(parents=True, exist_ok=True)
            args.record.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        print(json.dumps({"result": "FAILED", "problem": str(error),
                          "reversed": reversed_entries}, indent=2, sort_keys=True))
        print("All completed renames were reversed; accepted products are intact.", file=sys.stderr)
        return 1

    stays = [entry for entry in entries if entry["state"] == "migrated"]
    if stays:
        record["skipped_already_migrated"] = [entry["id"] for entry in stays]
    record["result"] = "PASS"
    record["finished"] = time.time()
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        print(f"record: {args.record}")
    print(json.dumps({"result": "PASS", "moved": len(record["moves"]),
                      "skipped": len(stays), "seconds": round(record["finished"] - record["started"], 2)},
                     indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
