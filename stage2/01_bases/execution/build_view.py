#!/usr/bin/env python3
"""Publish or independently replay the paired 5k-local expression release (ADR 0077).

Select unchanged G_5000 integer counts before CP10K and natural log1p. Each
vehicle is normalized separately before matched averaging. Build and verify the
complete paired candidate outside output before replacement. --output supports
candidate publication without modifying the live release. --repo-root binds
provenance when reading the owning materialized checkout from a worktree.
"""
import argparse
import json
import sys
from pathlib import Path

import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import view

STEP = Path(__file__).resolve().parents[1]
REPO = next(path for path in Path(__file__).resolve().parents if (path / ".git").exists())
VIEW_OUTPUTS = {size: STEP / f"output/view-{size // 1000}k" for size in view.VIEW_NAMES}


def view_name(size):
    return f"view-{int(size) // 1000}k"


def _relative(path):
    try:
        return Path(path).absolute().relative_to(REPO.absolute()).as_posix()
    except ValueError:
        return Path(path).name


def summarize(output, size, peak_rss_kib):
    manifest = json.loads((Path(output) / "manifest.json").read_text())
    return {
        "name": view_name(size),
        "path": _relative(output),
        "size": int(size),
        "manifest_sha256": view.sha256(Path(output) / "manifest.json"),
        "outputs": {
            name: {
                "sha256": record["sha256"],
                "bytes": record["bytes"],
                "rows": record.get("rows"),
                "columns": record.get("columns"),
            }
            for name, record in manifest["outputs"].items()
        },
        "sources": {
            name: {
                "rows": summary_["rows"],
                "selected_axis_library_total": summary_["selected_axis_library_total"],
            }
            for name, summary_ in manifest["sources"].items()
        },
        "schema": manifest["schema"],
        "normalization": manifest["normalization"],
        "axis": manifest["axis"],
        "peak_rss_kib": int(peak_rss_kib),
    }


def run(action, corpus, axis_dir, staging, output=None, repo_root=REPO):
    results = {}
    outputs = {5000: Path(output)} if output is not None else VIEW_OUTPUTS
    for size, output in sorted(outputs.items()):
        kwargs = {"corpus_dir": corpus, "axis_dir": axis_dir, "size": int(size),
                  "name": view_name(size), "repo_root": repo_root}
        if action == "publish":
            view.publish_release(output, root=staging, **kwargs)
        else:
            view.replay_release(output, root=staging, **kwargs)
        results[int(size)] = output
    import resource
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {"action": action,
            "views": {view_name(size): summarize(output, size, peak)
                      for size, output in sorted(results.items())}}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=("publish", "verify"))
    parser.add_argument("--corpus", type=Path, default=STEP / "output/corpus")
    parser.add_argument("--axis", type=Path, default=STEP / "output/gene-axis")
    parser.add_argument("--staging", type=Path,
                        default=STEP / "execution/checkpoints")
    parser.add_argument("--output", type=Path, help="paired 5k release destination")
    parser.add_argument("--repo-root", type=Path, default=REPO)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    result = run(args.action, args.corpus, args.axis, args.staging, args.output, args.repo_root)
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(text)
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
