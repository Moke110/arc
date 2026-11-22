#!/usr/bin/env python3
"""Plan, move, and recover the ignored Stage 2 payloads across ADR 0083."""
import argparse
import json
import os
from pathlib import Path
import stat
import subprocess

# Longer prefixes first so formerly study-owned feature bytes move to util.
ROOTS = {
    "stage2/03_tsg_drug_emb/output/representations": "stage2/util/output/representations",
    "stage2/06_tsg_transfer": "stage2/08_tsg_transfer",
    "stage2/05_tsg_width_dose": "stage2/07_tsg_width_dose",
    "stage2/04_tsg_dose": "stage2/06_tsg_dose",
    "stage2/03_tsg_drug_emb": "stage2/05_tsg_drug_emb",
    "stage2/02_models_tenv-5k": "stage2/04_models_tenv-5k",
    "stage2/tsg_view": "stage2/03_view/03.1_hvg-deg_views/historical/tsg_view",
    "stage2/01_bases": "stage2/01_bases",
}
SCHEMA = "stage2-payload-migration-v1"


def destination(relative):
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"invalid repository-relative path: {relative}")
    for old, new in ROOTS.items():
        if path.is_relative_to(old):
            return (Path(new) / path.relative_to(old)).as_posix()
    raise ValueError(f"unmapped Stage 2 path: {relative}")


def checked_path(root, relative):
    path = root / relative
    for parent in path.relative_to(root).parents:
        ancestor = root / parent
        if ancestor.is_symlink() or (ancestor.exists() and not ancestor.is_dir()):
            raise ValueError(f"unsafe payload parent: {ancestor}")
    return path


def identity(path):
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode):
        raise ValueError(f"payload must be a regular file: {path}")
    return [info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns]


def present(path):
    return path.exists() or path.is_symlink()


def git_files(root, options, roots=ROOTS):
    raw = subprocess.check_output(
        ["git", "-C", str(root), "ls-files", "-z", *options, "--", *roots]
    )
    return sorted({os.fsdecode(name) for name in raw.split(b"\0") if name})


def plan(root):
    entries = [{"old": name, "identity": identity(checked_path(root, name))}
               for name in git_files(root, ["--others", "--ignored", "--exclude-standard"])]
    return {"schema": SCHEMA, "root": str(root), "entries": entries}


def paths(root, entries, reverse=False):
    moves = []
    for entry in entries:
        old = checked_path(root, entry["old"])
        new = checked_path(root, destination(entry["old"]))
        if old == new:
            if identity(old) != entry["identity"]:
                raise ValueError(f"source changed: {old}")
            continue
        source, target = (new, old) if reverse else (old, new)
        if present(target):
            if present(source) or identity(target) != entry["identity"]:
                raise ValueError(f"destination exists: {target}")
        elif not present(source) or identity(source) != entry["identity"]:
            raise ValueError(f"source missing or changed: {source}")
        else:
            parent = next(p for p in target.parents if p.exists())
            if parent.stat().st_dev != source.stat().st_dev:
                raise ValueError(f"cannot rename across filesystems: {target}")
            moves.append((source, target))
    return moves


def move(root, entries, reverse):
    retired_roots = {key for key, value in ROOTS.items() if key != value}
    if not reverse and git_files(root, [], retired_roots):
        raise ValueError("old Stage 2 paths still contain tracked files; apply after the Git merge")
    if not reverse:
        extras = set(git_files(root, ["--others"], retired_roots)).difference(
            entry["old"] for entry in entries
        )
        if extras:
            raise ValueError(f"unplanned retired payload: {min(extras)}")
    moves = paths(root, entries, reverse)
    moved = []
    try:
        for source, target in moves:
            target.parent.mkdir(parents=True, exist_ok=True)
            if present(target):
                raise ValueError(f"destination exists: {target}")
            source.rename(target)
            moved.append((source, target))
    except (OSError, ValueError):
        for source, target in reversed(moved):
            target.rename(source)
        raise
    retired = ({ROOTS[old] for old in retired_roots} if reverse else retired_roots)
    roots = [root / name for name in retired]
    old_directories = {parent for entry in entries
                       for parent in checked_path(
                           root, destination(entry["old"]) if reverse else entry["old"]
                       ).parents
                       if any(parent.is_relative_to(directory) for directory in roots)}
    for directory in sorted(old_directories, key=lambda path: -len(path.parts)):
        if directory.is_dir() and not directory.is_symlink() and not any(directory.iterdir()):
            directory.rmdir()


def link(root, entries, checkout, historical):
    if checkout == root or not (checkout / ".git").exists():
        raise ValueError("link target must be a separate Git checkout")
    links = []
    for entry in entries:
        old = entry["old"]
        if "__pycache__" in Path(old).parts or Path(old).suffix in {".pyc", ".pyo"}:
            continue
        new = destination(old)
        source = checked_path(root, new)
        if not present(source):
            source = checked_path(root, old)
        if identity(source) != entry["identity"]:
            raise ValueError(f"source changed: {source}")
        target = checked_path(checkout, old if historical else new)
        if present(target):
            if (historical and target.is_file() and not target.is_symlink()
                    and target.samefile(source)) or (not historical and target.is_symlink()
                                                     and target.resolve() == source):
                continue
            raise ValueError(f"destination exists: {target}")
        links.append((source, target))
    linked = []
    try:
        for source, target in links:
            target.parent.mkdir(parents=True, exist_ok=True)
            if historical:
                target.hardlink_to(source)
            else:
                target.symlink_to(source)
            linked.append(target)
    except OSError:
        for target in reversed(linked):
            target.unlink()
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("plan", "apply", "reverse", "link-snapshot", "link-active"))
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--target", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    try:
        if args.action == "plan":
            record = plan(root)
            with args.plan.open("x") as stream:
                json.dump(record, stream, indent=2)
                stream.write("\n")
        else:
            record = json.loads(args.plan.read_text())
            if record["schema"] != SCHEMA or record["root"] != str(root):
                raise ValueError("plan belongs to a different schema or checkout")
            entries = record["entries"]
            sources = [item["old"] for item in entries]
            if len(sources) != len(set(sources)) or len(sources) != len({destination(s) for s in sources}):
                raise ValueError("duplicate mapping in plan")
            if args.action in {"apply", "reverse"}:
                move(root, entries, args.action == "reverse")
            else:
                if args.target is None:
                    raise ValueError("--target is required for replay links")
                link(root, entries, args.target.resolve(), args.action == "link-snapshot")
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"{error}\n")
    print(f"{args.action}: {len(record['entries'])} Stage 2 payloads")


if __name__ == "__main__":
    main()
