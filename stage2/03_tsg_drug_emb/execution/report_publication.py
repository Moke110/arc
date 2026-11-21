"""Verify an isolated report candidate or publish the reviewed report into this worktree."""
import argparse
import fcntl
import json
from pathlib import Path
from shutil import copyfile
import subprocess
import sys
from tempfile import TemporaryDirectory

from publish_report import FIGURES, RENDERER, STEP, figure_names, sha256, write_json

RECEIPT = "log/acceptance.json"
TAHOE_RECEIPT = "log/tahoe/acceptance.json"
TENV_README = Path("log/tenv/README.md")
TENV_PLOT_MOVES = {Path("figure") / name: Path("figure") / f"tenv_{name}" for name in FIGURES}
TENV_RENDERER_SHA256 = "223feccf854ba40d4e1ba303fbd13c9e57f84c849d60b658e4b18f9e437287b4"
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


def tahoe_files():
    return {
        "log/tahoe/README.md", "log/tahoe/run_inventory.json",
        *(f"figure/{name}" for name in figure_names("tahoe_bulk")),
        *(f"figure/data/tahoe/{table}.csv" for table in TABLES),
    }


def payload_inventory(candidate, receipt=RECEIPT):
    files = {}
    for path in candidate.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"report payload contains a symlink: {path}")
        relative = path.relative_to(candidate).as_posix()
        if path.is_file() and relative != receipt:
            files[relative] = sha256(path)
    return files


def published_path(report, name, source):
    path = report / name
    if source == "tenv":
        if name == "README.md" and ((report / TENV_README).exists() or (report / TENV_README).is_symlink()):
            return report / TENV_README
        moved = TENV_PLOT_MOVES.get(Path(name))
        if moved is not None and not path.exists() and not path.is_symlink():
            return report / moved
    return path


def verify_report(report, *, published=False, origins=None, source="tenv"):
    tahoe = source == "tahoe_bulk"
    receipt_name = TAHOE_RECEIPT if tahoe else RECEIPT
    receipt_path = report / receipt_name
    if receipt_path.is_symlink():
        raise ValueError("report receipt must be a regular file")
    receipt = json.loads(receipt_path.read_text())
    schema = "tsg-drug-emb-report-acceptance-v1" if published else "tsg-drug-emb-report-candidate-v1"
    status = "PASS" if published else "candidate"
    allowed_renderers = {sha256(RENDERER)}
    if published and source == "tenv":
        allowed_renderers.add(TENV_RENDERER_SHA256)
    if (receipt.get("schema") != schema or receipt.get("status") != status
            or receipt.get("source", "tenv") != source
            or receipt.get("renderer_sha256") not in allowed_renderers
            or set(receipt.get("files", ())) != (tahoe_files() if tahoe else REPORT_FILES)):
        raise ValueError("report receipt schema, inventory or renderer differs")
    if origins is not None:
        for key, path in origins.items():
            if sha256(path) != receipt[key]:
                raise ValueError(f"report origin differs: {key}")
    if published:
        for name, digest in receipt["files"].items():
            path = published_path(report, name, source)
            parents = path.relative_to(report).parents
            if (path.is_symlink() or any((report / part).is_symlink() for part in parents)
                    or not path.is_file() or sha256(path) != digest):
                raise ValueError(f"published report payload differs: {name}")
    elif payload_inventory(report, receipt_name) != receipt["files"]:
        raise ValueError("candidate report payload differs from its sealed inventory")
    with TemporaryDirectory(prefix="tsg-drug-emb-report-replay-") as directory:
        replay = Path(directory)
        data = report / "figure/data"
        if tahoe:
            data /= "tahoe"
        command = [sys.executable, str(RENDERER), "--data", str(data), "--output", directory,
                   "--source", source]
        if tahoe:
            command.append("--prefixed")
        subprocess.run(command, check=True)
        figures = figure_names(source)
        if {path.name for path in replay.iterdir()} != set(figures):
            raise ValueError("CSV-only replay produced an unexpected figure inventory")
        for figure in figures:
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
    source = receipt.get("source", "tenv")
    tahoe = source == "tahoe_bulk"
    receipt_name = TAHOE_RECEIPT if tahoe else RECEIPT
    files = tahoe_files() if tahoe else REPORT_FILES
    readme = step / "README.md"
    targets = (*sorted(files - {"README.md"}), receipt_name)
    for name in targets:
        target = step / name
        if (target.exists() or target.is_symlink()
                or any((step / parent).is_symlink() for parent in Path(name).parents)):
            raise FileExistsError(f"report target already exists or is symlinked: {target}")
    if not tahoe and (readme.is_symlink() or sha256(readme) != DRAFT_README_SHA256):
        raise ValueError("Stage 2 README differs from the unpublished draft")
    if tahoe:
        verify_report(step, published=True)
        for new in (TENV_README, *TENV_PLOT_MOVES.values()):
            target = step / new
            if (target.exists() or target.is_symlink()
                    or any((step / parent).is_symlink() for parent in new.parents)):
                raise FileExistsError(f"tenv migration target exists or is symlinked: {target}")
    with TemporaryDirectory(prefix=".report-staging-", dir=step) as directory:
        staging = Path(directory)
        previous = staging / "draft.md"
        copyfile(readme, previous)
        if tahoe:
            active = readme.read_text()
            for old, new in TENV_PLOT_MOVES.items():
                active = active.replace(f"]({old.as_posix()})", f"]({new.as_posix()})")
            active += "\nThe [Tahoe comparison](log/tahoe/README.md) uses the three frozen Tahoe cold-drug folds.\n"
            (staging / "README.md").write_text(active)
        for name in files:
            target = staging / name
            target.parent.mkdir(parents=True, exist_ok=True)
            copyfile(candidate / name, target)
            if sha256(target) != receipt["files"][name]:
                raise ValueError(f"candidate changed during copy: {name}")
        write_json(staging / receipt_name, dict(
            receipt, schema="tsg-drug-emb-report-acceptance-v1", status="PASS",
            candidate_receipt_sha256=sha256(candidate / receipt_name),
            visual_review=visual_review,
        ))
        installed = []
        moved = []
        try:
            if tahoe:
                for old, new in {Path("README.md"): TENV_README, **TENV_PLOT_MOVES}.items():
                    (step / new).parent.mkdir(parents=True, exist_ok=True)
                    (step / old).rename(step / new)
                    moved.append((old, new))
            for name in targets:
                target = step / name
                target.parent.mkdir(parents=True, exist_ok=True)
                (staging / name).rename(target)
                installed.append(target)
            if not tahoe:
                (staging / "README.md").replace(readme)
            else:
                (staging / "README.md").rename(readme)
                verify_report(step, published=True)
                expected = {*figure_names("tahoe_bulk"), *(path.name for path in TENV_PLOT_MOVES.values())}
                if {path.name for path in (step / "figure").glob("*.png")} != expected:
                    raise ValueError("paired report must contain exactly six source-prefixed PNGs")
            verify_report(step, published=True, source=source)
        except BaseException:
            for path in reversed(installed):
                path.unlink()
            for old, new in reversed(moved):
                (step / new).replace(step / old)
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
    source = json.loads((args.inputs / "config.json").read_text())["source"]
    receipt = verify_report(report, published=published, origins=origins, source=source)
    if args.command == "publish":
        install_report(report, STEP, receipt, args.visual_review)
    print(f"{args.command}: verified report digests, origins and fresh CSV-only replay", flush=True)


if __name__ == "__main__":
    main()
