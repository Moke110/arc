"""Fast-forward a materialized checkout across the ADR 0099 Stage 2 rename."""

import argparse
import fcntl
import json
import os
from pathlib import Path
import subprocess

from stage2_history import check_ignored, exists, file_record, git, paths_under, write_json

OLD = Path("stage2")
NEW = Path("stage2_CarPS-VAE")


def tree(root, revision, stage):
    result = {}
    for row in git(root, "ls-tree", "-r", "-z", revision, "--", str(stage)).split(b"\0"):
        if row:
            metadata, path = row.split(b"\t", 1)
            result[Path(os.fsdecode(path)).relative_to(stage).as_posix()] = metadata.decode()
    return result


def inventory(stage):
    if not stage.is_dir() or stage.is_symlink():
        raise ValueError(f"expected an ordinary Stage 2 directory: {stage}")
    files, directories = paths_under(stage, Path("."))
    result = {}
    for name in sorted(files | directories):
        path = stage / name
        if name in directories:
            info = path.lstat()
            record = {"directory": True, "device": info.st_dev, "inode": info.st_ino, "mode": info.st_mode}
        else:
            record = file_record(path)
            if "link" in record:
                target = path.stat()
                record["target"] = [target.st_dev, target.st_ino]
        result[name.as_posix()] = record
    return result


def unchanged(before, after, changed):
    for name, record in before.items():
        if name not in changed and after.get(name) != record:
            raise ValueError(f"materialized Stage 2 entry changed: {name}")
    for name in after.keys() - before.keys() - set(changed):
        if not after[name].get("directory"):
            raise ValueError(f"unexpected materialized Stage 2 entry: {name}")


def check_primary_branch(root):
    if git(root, "branch", "--show-current").strip() != b"main":
        raise ValueError("run the CarPS cutover from the checked-out main branch")


def prepare(root, target):
    check_primary_branch(root)
    if exists(root / NEW):
        raise ValueError(f"destination already exists: {root / NEW}")
    if git(root, "status", "--porcelain", "--untracked-files=no", "--", str(OLD)):
        raise ValueError(f"tracked {OLD} must match the predecessor")
    git(root, "diff", "--cached", "--quiet")
    before = git(root, "rev-parse", "HEAD").decode().strip()
    target = git(root, "rev-parse", "--verify", "--end-of-options", target + "^{commit}").decode().strip()
    if target != git(root, "rev-parse", "--verify", "refs/heads/dev").decode().strip():
        raise ValueError("the target must be the exact dev tip")
    git(root, "merge-base", "--is-ancestor", before, target)
    old, new = tree(root, before, OLD), tree(root, target, NEW)
    if not old or not old.keys() <= new.keys() or tree(root, target, OLD):
        raise ValueError("target must relocate every previous tracked Stage 2 file")
    entries = inventory(root / OLD)
    for name in new.keys() - old.keys():
        for path in (Path(name), *Path(name).parents):
            key = path.as_posix()
            if key in entries and (key == name or not entries[key].get("directory")):
                raise ValueError(f"tracked target overlaps materialized data: {key}")
    payload = sorted(name for name, record in entries.items()
                     if not record.get("directory") and name not in old)
    return {"before": before, "target": target, "entries": entries, "payload": payload,
            "changed": sorted(name for name, identity in new.items() if old.get(name) != identity)}


def git_head(root):
    return git(root, "rev-parse", "HEAD").decode().strip()


def finish(root, work, plan):
    check_primary_branch(root)
    old, new = root / OLD, root / NEW
    parked, checkout = work / "materialized", work / "target-tracked"
    if git_head(root) == plan["before"]:
        if not exists(parked):
            os.rename(old, parked)
        try:
            git(root, "merge", "--quiet", "--ff-only", plan["target"])
        except subprocess.CalledProcessError:
            if git_head(root) == plan["before"] and not exists(old):
                os.rename(parked, old)
            raise
    if git_head(root) != plan["target"] or exists(old):
        raise ValueError("checkout is outside the recorded cutover revisions or still has the predecessor root")

    check_ignored(root, {NEW / name for name in plan["payload"]}, "CarPS")
    if exists(parked):
        if not exists(checkout):
            os.rename(new, checkout)
        os.rename(parked, new)
    for name in plan["changed"]:
        source, target = checkout / name, new / name
        if exists(source):
            target.parent.mkdir(parents=True, exist_ok=True)
            os.replace(source, target)


def verify(root, plan):
    if exists(root / OLD):
        raise ValueError("the previous Stage 2 directory still exists")
    git(root, "merge-base", "--is-ancestor", plan["target"], "HEAD")
    check_ignored(root, {NEW / name for name in plan["payload"]}, "CarPS")
    if git(root, "status", "--porcelain", "--untracked-files=no", "--", str(NEW)):
        raise ValueError("tracked Stage 2 files differ from Git")
    current = inventory(root / NEW)
    unchanged(plan["entries"], current, plan["changed"])
    preserved = plan["entries"].keys() - set(plan["changed"])
    return {"status": "PASS", "before": plan["before"], "target": plan["target"],
            "materialized_entries_preserved": len(preserved),
            "archive_entries": sum(name == "history" or name.startswith("history/") for name in plan["entries"]),
            "external_links": sum("link" in record for record in plan["entries"].values()),
            "tracked_paths_updated": len(plan["changed"]), "ignored_payload_files": len(plan["payload"])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("action", choices=("deploy", "resume", "verify"))
    parser.add_argument("--target", help="fast-forward revision for deploy")
    args = parser.parse_args()
    root, work = args.root.resolve(), args.work.resolve()
    if work == root / "tmp" or not work.is_relative_to(root / "tmp"):
        parser.error("cutover work must use a subdirectory of checkout tmp/")
    try:
        if args.action == "deploy":
            if args.target is None:
                parser.error("deploy requires --target")
            work.mkdir(parents=True)
        if root.stat().st_dev != work.stat().st_dev:
            raise ValueError("the checkout and cutover work must share one filesystem")
        with (work / "lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            state = work / "before.json"
            if args.action == "deploy":
                plan = prepare(root, args.target)
                write_json(state, plan)
            else:
                plan = json.loads(state.read_text())
            if args.action != "verify":
                finish(root, work, plan)
            result = verify(root, plan)
            write_json(work / "after.json", result)
        print(json.dumps(result, sort_keys=True))
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"CarPS Stage 2 cutover: {error}\n")


if __name__ == "__main__":
    main()
