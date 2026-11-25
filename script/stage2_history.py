"""Install the ADR 0087 local history from a materialized Stage 2 checkout."""
import argparse
import fcntl
from hashlib import file_digest
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
from tempfile import TemporaryDirectory

HISTORY = Path("stage2/history")
TRAINING = Path("stage2/03_data_train")
STUDY = Path("stage2/04_tsg_axis")
ARCHIVE_ROOTS = tuple(Path(name) for name in (
    "stage2/01_data_raw/historical",
    "stage2/05_tsg_loss", "stage2/06_tsg_drug_emb", "stage2/07_tsg_dose",
    "stage2/08_tsg_width_dose", "stage2/09_tsg_transfer",
    "tmp/stage2-297-rollback",
    "tmp/stage2-305-deploy/canonical-run/predecessor-roots",
    "tmp/stage2-305-deploy/predecessor-payloads.json",
    "tmp/stage2-305-deploy/cutover-notes.md",
    "tmp/stage2-305-deploy/relocate_payloads.py",
    "tmp/stage2-302-deploy/predecessor-history",
    "tmp/stage2-302-deploy/predecessor-raw-first-v1",
    "tmp/stage2-302-deploy/retired-step",
    "tmp/stage2-302-deploy/migration.json",
    "tmp/stage2-302-deploy/predecessor-tracked-history.sha256",
)) + tuple(TRAINING / name for name in (
    "03.1_hvg-deg-5k", ".output.publication.lock", ".output.previous",
    "execution/build_training.py", "execution/tests/test_training.py",
    "output/api/training_data.py", "output/manifest.json", "output/tenv", "output/tahoe-bulk",
    "log/305-review-promotion.json", "log/307-tenv-candidate.json",
    "log/308-tahoe-paired-candidate.json", "log/309-cutover.json",
)) + tuple(STUDY / name for name in (
    "README.md", "historical", "output/formal", "output/report", "output/.report.publication.lock",
    "figure/data/delta-metrics.csv", "figure/data/disease-f1.csv",
    "figure/delta-metrics.png", "figure/disease-f1.png", "figure/script/render.py",
    "log/316-runs", "log/310-runtime-replay.json", "log/314-smoke-and-migration.json",
    "log/315-metrics-smoke.json", "log/316-formal.json", "log/delta-metrics-style.json",
    "execution/verify_historical.py", "execution/replay_pinned.py",
))
RELEASES = {
    TRAINING / name: TRAINING / name.replace("-deg-", "-")
    for name in (
        "03.1_tenv_hvg-deg-5k", "03.2_tahoe_hvg-deg-5k", "03.3_tenv_hvg-deg-3k",
        "03.4_tenv_hvg-deg-4k", "03.5_tenv_denoise-deg-3k", "03.6_tenv_denoise-deg-4k",
        "03.7_tenv_denoise-deg-5k",
    )
}
REPLAY = {
    "retired_bases": "130afdd8b96188fc714f586d65461956e21cdc79",
    "full_view": "f15cb32174c291800eb7017fbd19d460bb321e26",
    "paired_training": "b28a35fc03370ce45c06155a4c472366eb7aebd8",
    "axis_wide_control": "986f8eacb2cda5b9e0793c03e12a27761aaa1c36",
    "axis_wide_control_figures": "85bc014",
    "matched_control_report": "3f4cbb6569b739c6c3761efb1f9c0dde254f9be3",
}
SMALL_BYTES = 1024 * 1024


def git(root, *args):
    return subprocess.check_output(("git", *args), cwd=root)


def names(raw):
    return {Path(os.fsdecode(name)) for name in raw.split(b"\0") if name}


def exists(path):
    return os.path.lexists(path)


def check_parents(root, relative):
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"invalid repository-relative archive path: {relative}")
    for parent in relative.parents:
        if (root / parent).is_symlink():
            raise ValueError(f"archive path has a symlink ancestor: {root / parent}")


def sha256(path):
    with path.open("rb") as stream:
        return file_digest(stream, "sha256").hexdigest()


def file_record(path, *, force_digest=False):
    info = path.lstat()
    result = {"device": info.st_dev, "inode": info.st_ino, "mode": info.st_mode,
              "size": info.st_size, "mtime_ns": info.st_mtime_ns}
    if stat.S_ISLNK(info.st_mode):
        result["link"] = os.readlink(path)
    elif stat.S_ISREG(info.st_mode):
        if force_digest or info.st_size <= SMALL_BYTES:
            result["sha256"] = sha256(path)
    else:
        raise ValueError(f"unsupported archive payload type: {path}")
    return result


def paths_under(root, source):
    path = root / source
    if not exists(path):
        return set(), set()
    if path.is_symlink() or not path.is_dir():
        return {source}, set()
    files, dirs = set(), {source}
    for directory, subdirs, leaves in os.walk(path):
        current = Path(directory)
        for name in subdirs.copy():
            child = current / name
            if child.is_symlink():
                subdirs.remove(name)
                leaves.append(name)
            else:
                dirs.add(child.relative_to(root))
        files.update((current / name).relative_to(root) for name in leaves)
    return files, dirs


def write_json(path, record):
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")


def prepare(root, state, revision):
    if exists(state) or exists(root / HISTORY):
        raise ValueError("prepare requires an unused state and absent Stage 2 history")
    if git(root, "status", "--porcelain", "--untracked-files=no", "--", "stage2", ".gitignore"):
        raise ValueError("prepare requires a clean tracked Stage 2 checkout")
    before = git(root, "rev-parse", "HEAD").decode().strip()
    target = git(root, "rev-parse", "--verify", "--end-of-options", revision + "^{commit}").decode().strip()
    tracked = names(git(root, "ls-files", "-z"))
    successor = names(git(root, "ls-tree", "-r", "-z", "--name-only", target))
    if any(path.is_relative_to(HISTORY) for path in successor):
        raise ValueError("the successor's Git tree contains Stage 2 history")
    if b"/stage2/history/" not in git(root, "show", f"{target}:.gitignore").splitlines():
        raise ValueError("the successor must ignore the Stage 2 history directory")
    state.parent.mkdir(parents=True, exist_ok=True)
    if root.stat().st_dev != state.parent.stat().st_dev:
        raise ValueError("cutover state must use the checkout filesystem")

    entries, directories, groups = {}, set(), {}
    with TemporaryDirectory(prefix="stage2-history-", dir=state.parent) as temporary:
        candidate = Path(temporary) / "state"
        archive = candidate / "archive"
        archive.mkdir(parents=True)
        for source in (*ARCHIVE_ROOTS, *RELEASES):
            check_parents(root, source)
            if (root / source).is_symlink():
                raise ValueError(f"archive root must not be a symlink: {source}")
            files, dirs = paths_under(root, source)
            if not files and not dirs:
                continue
            release = RELEASES.get(source)
            groups[source.as_posix()] = (release or HISTORY / source).as_posix()
            if release is None:
                directories.update(dirs)
            for path in sorted(files):
                original = root / path
                record = file_record(original, force_digest=path in tracked)
                if record["device"] != root.stat().st_dev:
                    raise ValueError(f"source is on another filesystem: {path}")
                if release is not None:
                    new_path = release / path.relative_to(source)
                    if path in tracked:
                        if new_path not in successor or path in successor:
                            raise ValueError(f"unverified tracked release rename: {path}")
                        if path.name == "manifest.json":
                            old_blob = git(root, "rev-parse", f"{before}:{path}").strip()
                            new_blob = git(root, "rev-parse", f"{target}:{new_path}").strip()
                            if old_blob != new_blob:
                                raise ValueError(f"release manifest changes in the successor: {path}")
                        kind = "release-tracked"
                    else:
                        kind = "release-payload"
                elif path in tracked:
                    kind = "retained-snapshot" if path == source and path in successor else "snapshot"
                    if kind == "snapshot" and path in successor:
                        raise ValueError(f"successor still tracks an archival path: {path}")
                    new_path = HISTORY / path
                    snapshot = archive / path
                    snapshot.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(original, snapshot, follow_symlinks=False)
                    check_snapshot(snapshot, record)
                else:
                    kind, new_path = "archive-payload", HISTORY / path
                entries[path.as_posix()] = {"kind": kind, "destination": new_path.as_posix(), **record}
        if not entries:
            raise ValueError("there are no materialized archive or release paths")
        for directory in directories:
            (archive / directory).mkdir(parents=True, exist_ok=True)
        record = {
            "schema": "stage2-local-history-v1",
            "predecessor_git_revision": before,
            "target_git_revision": target,
            "pinned_replay_revisions": REPLAY,
            "path_map": groups,
            "archive_directories": sorted(path.as_posix() for path in directories),
            "entries": entries,
        }
        write_json(candidate / "plan.json", record)
        candidate.rename(state)
    return record


def check_snapshot(path, record):
    current = file_record(path, force_digest=True)
    for key in {"mode", "size", "mtime_ns", "sha256", "link"} & record.keys():
        if current.get(key) != record[key]:
            raise ValueError(f"archived tracked snapshot differs: {path}")


def check_payload(path, record):
    current = file_record(path)
    for key, value in current.items():
        if record.get(key) != value:
            raise ValueError(f"payload identity or bytes differ: {path}")


def check_archive_boundary(root, tracked):
    if any(path.is_relative_to(HISTORY) for path in tracked):
        raise ValueError("the archive or one of its entries is tracked")
    ignore = subprocess.run(("git", "check-ignore", "-q", "--", str(HISTORY / "archive-manifest.json")), cwd=root)
    if ignore.returncode:
        raise ValueError(f"Git ignore does not cover {root / HISTORY}")


def checked_locations(root, state, plan, tracked):
    live = root / HISTORY
    staged = state / "archive"
    check_parents(root, HISTORY)
    if (any(exists(path) and (path.is_symlink() or not path.is_dir()) for path in (live, staged))
            or live.is_dir() == staged.is_dir()):
        raise ValueError("exactly one plain Stage 2 archive directory must exist")
    archive = live if live.is_dir() else staged
    expected_files = {"archive-manifest.json"} if exists(archive / "archive-manifest.json") else set()
    locations = []
    for raw_path, record in plan["entries"].items():
        source = Path(raw_path)
        destination = Path(record["destination"])
        kind = record["kind"]
        check_parents(root, source)
        original = root / source
        if destination.is_relative_to(HISTORY):
            relative = destination.relative_to(HISTORY)
            check_parents(archive, relative)
            replacement = archive / relative
            if kind != "archive-payload" or exists(replacement):
                expected_files.add(relative.as_posix())
        else:
            check_parents(root, destination)
            replacement = root / destination
        if kind in {"snapshot", "retained-snapshot"}:
            check_snapshot(replacement, record)
            if (kind == "retained-snapshot") != (source in tracked):
                raise ValueError(f"unexpected tracked snapshot state: {source}")
        elif kind == "release-tracked":
            if destination not in tracked or not replacement.is_file():
                raise ValueError(f"canonical tracked release file is missing: {destination}")
            if source.name == "manifest.json" and sha256(replacement) != record["sha256"]:
                raise ValueError(f"canonical release manifest changed: {destination}")
        else:
            if not exists(original) and not exists(replacement):
                raise ValueError(f"both original and replacement payload are missing: {source}")
            if exists(original):
                check_payload(original, record)
            if exists(replacement):
                check_payload(replacement, record)
            locations.append((original, replacement))
    archive_files, _ = paths_under(archive, Path("."))
    if {path.as_posix() for path in archive_files} != expected_files:
        raise ValueError("unexpected file in the local archive")
    for name in plan["path_map"]:
        current, _ = paths_under(root, Path(name))
        for path in current:
            record = plan["entries"].get(path.as_posix())
            if not record or record["kind"] not in {
                "retained-snapshot", "archive-payload", "release-payload"
            }:
                raise ValueError(f"unexpected or residual file after Git cutover: {path}")
    return archive, locations


def move(source, destination):
    if not exists(source):
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not exists(destination):
        os.link(source, destination, follow_symlinks=False)
    if (source.lstat().st_dev, source.lstat().st_ino) != (destination.lstat().st_dev, destination.lstat().st_ino):
        raise ValueError(f"payload destination is occupied: {destination}")
    source.unlink()


def prune_sources(root, plan):
    for name in plan["path_map"]:
        _, dirs = paths_under(root, Path(name))
        for directory in sorted(dirs, key=lambda item: len(item.parts), reverse=True):
            path = root / directory
            if not any(path.iterdir()):
                path.rmdir()


def apply(root, state, action):
    plan = json.loads((state / "plan.json").read_text())
    if plan["schema"] != "stage2-local-history-v1":
        raise ValueError("unrecognized Stage 2 archive plan")
    git(root, "merge-base", "--is-ancestor", plan["target_git_revision"], "HEAD")
    if git(root, "status", "--porcelain", "--untracked-files=no", "--", "stage2", ".gitignore"):
        raise ValueError("archive install requires a clean tracked Stage 2 checkout")
    tracked = names(git(root, "ls-files", "-z"))
    check_archive_boundary(root, tracked)
    archive, locations = checked_locations(root, state, plan, tracked)
    live = root / HISTORY
    if action == "undo":
        for original, replacement in reversed(locations):
            move(replacement, original)
        if archive == live:
            archive.rename(state / "archive")
        return plan
    if action == "install":
        for original, replacement in locations:
            move(original, replacement)
        prune_sources(root, plan)
        if archive != live:
            write_json(archive / "archive-manifest.json", plan)
            archive.rename(live)
    archive, locations = checked_locations(root, state, plan, tracked)
    if archive != live or any(exists(original) or not exists(replacement) for original, replacement in locations):
        raise ValueError("Stage 2 archive installation is incomplete")
    if json.loads((archive / "archive-manifest.json").read_text()) != plan:
        raise ValueError("installed archive manifest differs from the cutover inventory")
    if any(not (archive / name).is_dir() for name in plan["archive_directories"]):
        raise ValueError("installed archive is missing a predecessor directory")
    return plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "install", "verify", "undo"))
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True, help="cutover state under the checkout's tmp/")
    parser.add_argument("--target", help="successor Git revision for prepare")
    args = parser.parse_args()
    root, state = args.root.resolve(), args.state.resolve()
    if state == root / "tmp" or not state.is_relative_to(root / "tmp"):
        parser.error("cutover state must be in a subdirectory of the checkout's tmp/")
    try:
        if args.action == "prepare":
            if args.target is None:
                parser.error("--target is required for prepare")
            plan = prepare(root, state, args.target)
        else:
            with (state / "lock").open("a") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                plan = apply(root, state, args.action)
        print(json.dumps({"result": "PASS", "action": args.action, "files": len(plan["entries"]),
                          "roots": plan["path_map"], "target": plan["target_git_revision"]}, sort_keys=True))
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"stage2 history: {error}\n")


if __name__ == "__main__":
    main()
