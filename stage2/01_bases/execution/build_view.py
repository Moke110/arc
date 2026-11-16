#!/usr/bin/env python3
"""Publish or replay the four Stage 2 model-facing views (ADR 0073).

Each view is the corpus normalized once on the **full 18,164-gene library total**
(CP10K, then ``log1p``), then subset to one released axis's genes, with
``control`` = the mean of the row's matched corpus controls normalized the same
way and ``delta = X - control``. The denominator is never the view's own gene
subset.

This driver only binds paths; every scientific definition and the public reader
live once in ``output/api/bases/view.py``.

    python stage2/01_bases/execution/build_view.py publish
    python stage2/01_bases/execution/build_view.py verify

Publication is candidate-then-replace (ADR 0062): each view is built in a
candidate under ``execution/checkpoints/``, verified, and only then swapped in;
a failed build leaves the accepted release usable. ``verify`` replays all four views
from the pinned corpus and axis inputs and compares the manifest and every output
digest. One view at a time; the corpus is streamed in bounded batches.
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
                "library_total": summary_["full_axis_library_total"],
                "basis_share": summary_["axis_share_of_full_axis_cp10k"],
                "basis_share_treated": summary_["treated_axis_share_of_full_axis_cp10k"],
            }
            for name, summary_ in manifest["sources"].items()
        },
        "basis": manifest["basis"]["sources"],
        "peak_rss_kib": int(peak_rss_kib),
    }


def run(action, corpus, axis_dir, staging):
    results = {}
    for size, output in sorted(VIEW_OUTPUTS.items()):
        kwargs = {"corpus_dir": corpus, "axis_dir": axis_dir, "size": int(size),
                  "name": view_name(size), "repo_root": REPO}
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
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    result = run(args.action, args.corpus, args.axis, args.staging)
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(text)
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
