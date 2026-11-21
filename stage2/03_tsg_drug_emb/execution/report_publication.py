"""Verify an isolated report candidate or publish the reviewed report into this worktree."""
import argparse
import fcntl
import json
from pathlib import Path
from shutil import copyfile
import subprocess
import sys
from tempfile import TemporaryDirectory

from publish_report import FIGURES, RENDERER, STEP, sha256, write_json

RECEIPT = "log/acceptance.json"
DRAFT_README_SHA256 = "f48ca0035407f9921a513a76d7e156edf17b72f36ea522ff00628ef2541242b7"
TABLES = (
    "runs", "summary", "dose", "curves", "parent_pcc", "pooled",
    "paired_parents", "paired", "axis_segments",
)
REPORT_FILES = {
    "README.md", "log/run_inventory.json",
    *(f"figure/{name}" for name in FIGURES),
    *(f"figure/data/{table}.csv" for table in TABLES),
}


def payload_inventory(candidate):
    files = {}
    for path in candidate.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"report payload contains a symlink: {path}")
        relative = path.relative_to(candidate).as_posix()
        if path.is_file() and relative != RECEIPT:
            files[relative] = sha256(path)
    return files


def verify_report(report, *, published=False, origins=None):
    receipt_path = report / RECEIPT
    if receipt_path.is_symlink():
        raise ValueError("report receipt must be a regular file")
    receipt = json.loads(receipt_path.read_text())
    schema = "tsg-drug-emb-report-acceptance-v1" if published else "tsg-drug-emb-report-candidate-v1"
    status = "PASS" if published else "candidate"
    if (receipt.get("schema") != schema or receipt.get("status") != status
            or receipt.get("renderer_sha256") != sha256(RENDERER)
            or set(receipt.get("files", ())) != REPORT_FILES):
        raise ValueError("report receipt schema, inventory or renderer differs")
    if origins is not None:
        for key, path in origins.items():
            if sha256(path) != receipt[key]:
                raise ValueError(f"report origin differs: {key}")
    if published:
        for name, digest in receipt["files"].items():
            path = report / name
            parents = path.relative_to(report).parents
            if (path.is_symlink() or any((report / part).is_symlink() for part in parents)
                    or not path.is_file() or sha256(path) != digest):
                raise ValueError(f"published report payload differs: {name}")
    elif payload_inventory(report) != receipt["files"]:
        raise ValueError("candidate report payload differs from its sealed inventory")
    with TemporaryDirectory(prefix="tsg-drug-emb-report-replay-") as directory:
        replay = Path(directory)
        subprocess.run([sys.executable, str(RENDERER), "--data", str(report / "figure/data"),
                        "--output", directory], check=True)
        if {path.name for path in replay.iterdir()} != set(FIGURES):
            raise ValueError("CSV-only replay produced an unexpected figure inventory")
        for figure in FIGURES:
            if sha256(replay / figure) != receipt["files"][f"figure/{figure}"]:
                raise ValueError(f"fresh CSV-only replay differs: {figure}")
    return receipt


def install_report(candidate, step, receipt, visual_review):
    if not visual_review.strip():
        raise ValueError("record the visual review of the three report figures before publication")
    if not step.is_dir() or step.is_symlink():
        raise ValueError("report publication requires an existing Stage 2 step directory")
    lock_path = step / "execution/checkpoints/.report-publication.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        _install_under_lock(candidate, step, receipt, visual_review)


def _install_under_lock(candidate, step, receipt, visual_review):
    readme = step / "README.md"
    targets = (RECEIPT, *sorted(REPORT_FILES - {"README.md"}))
    for name in targets:
        target = step / name
        if (target.exists() or target.is_symlink()
                or any((step / part).is_symlink() for part in target.relative_to(step).parents)):
            raise FileExistsError(f"report target already exists or is symlinked: {target}")
    if readme.is_symlink() or sha256(readme) != DRAFT_README_SHA256:
        raise ValueError("Stage 2 README differs from the unpublished draft")
    with TemporaryDirectory(prefix=".report-staging-", dir=step) as directory:
        staging = Path(directory)
        previous = staging / "draft.md"
        copyfile(readme, previous)
        for name in REPORT_FILES:
            target = staging / name
            target.parent.mkdir(parents=True, exist_ok=True)
            copyfile(candidate / name, target)
            if sha256(target) != receipt["files"][name]:
                raise ValueError(f"candidate changed during copy: {name}")
        write_json(staging / RECEIPT, dict(
            receipt, schema="tsg-drug-emb-report-acceptance-v1", status="PASS",
            candidate_receipt_sha256=sha256(candidate / RECEIPT),
            visual_review=visual_review,
        ))
        installed = []
        try:
            for name in targets[1:]:
                target = step / name
                target.parent.mkdir(parents=True, exist_ok=True)
                (staging / name).rename(target)
                installed.append(target)
            (staging / "README.md").replace(readme)
            (staging / RECEIPT).rename(step / RECEIPT)
            installed.append(step / RECEIPT)
            verify_report(step, published=True)
        except BaseException:
            for path in reversed(installed):
                path.unlink()
            previous.replace(readme)
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for command in ("verify-candidate", "verify-published", "publish"):
        group = commands.add_parser(command)
        if command != "verify-published":
            group.add_argument("--candidate", required=True, type=Path)
        for name in ("inputs", "representations", "matrix"):
            group.add_argument(f"--{name}", type=Path, required=True)
        if command == "publish":
            group.add_argument("--visual-review", required=True)
    args = parser.parse_args()
    origins = {
        "formal_acceptance_sha256": args.matrix / "acceptance.json",
        "inputs_sha256": args.inputs / "manifest.json",
        "representations_sha256": args.representations / "manifest.json",
    }
    published = args.command == "verify-published"
    report = STEP if published else args.candidate.resolve()
    receipt = verify_report(report, published=published, origins=origins)
    if args.command == "publish":
        install_report(report, STEP, receipt, args.visual_review)
    print(f"{args.command}: verified report digests, origins and fresh CSV-only replay", flush=True)


if __name__ == "__main__":
    main()
