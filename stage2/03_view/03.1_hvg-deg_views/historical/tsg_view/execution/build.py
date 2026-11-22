#!/usr/bin/env python3
"""Build one ADR 0082 expression-view product.

The CLI resolves the pinned inputs and a staged destination, then calls
``views/<view>.py::build(context)``. A variant module owns only its recipe; the shared
harness owns the product schema, rows, payload layout, companion tables and digests, and
this driver owns staging, installation and the build receipt.

Only ``full`` is implemented in ticket 01; tickets 02-04 add ``merge``, ``solo``,
``solo_deepest`` and ``full_pooled_label`` modules without touching this file or the
harness.
"""
import argparse
import importlib
import json
import resource
import shutil
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEP = HERE.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from views import (MODES, SOURCES, VIEWS, Context, HarnessError,  # noqa: E402
                   environment, relative, sha256)

RECEIPT_SCHEMA = "tsg-view-build-receipt-v1"
INTERPRETER_WARNING = (
    "use /home/chang/miniconda3/envs/arc/bin/python: other interpreters do not reproduce "
    "the accepted digests")


def variant_module(view):
    """Import ``views/<view>.py`` (hyphens map to underscores); fail closed if absent."""
    name = view.replace("-", "_")
    path = HERE / "views" / f"{name}.py"
    if not path.is_file():
        raise HarnessError(
            f"view {view!r} has no builder module: expected {path}; ticket 01 implements "
            "only 'full'")
    module = importlib.import_module(f"views.{name}")
    if not callable(getattr(module, "build", None)):
        raise HarnessError(f"{path} does not define build(context)")
    return module


def build_product(*, view, source, output, staging, repo_root_path=None, size=5000,
                  fixture=False, corpus=None, axis=None, accepted=None,
                  provenance_root="repo", keep_staging=False):
    """Build one source's product into ``output/<view>/`` and return its receipt."""
    if view not in VIEWS:
        raise HarnessError(f"unknown view: {view}")
    if source not in SOURCES:
        raise HarnessError(f"unknown source: {source}")
    output, staging = Path(output), Path(staging)
    staging.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    work = Path(tempfile.mkdtemp(prefix=f".{view}-{source}-", dir=staging))
    context = Context(view=view, source=source, output_dir=output, staging_dir=work,
                      repo_root_path=repo_root_path, size=size, fixture=fixture,
                      corpus_dir=corpus, axis_dir=axis, accepted_dir=accepted,
                      provenance_root=provenance_root)
    module = variant_module(view)
    try:
        manifest = module.build(context)
        installed = context.install()
    except BaseException:
        if not keep_staging:
            shutil.rmtree(work, ignore_errors=True)
        raise
    if not keep_staging:
        shutil.rmtree(work, ignore_errors=True)
    product_dir = output / view
    manifest_path = product_dir / f"{source}.json"
    receipt = {
        "schema": RECEIPT_SCHEMA,
        "view": view,
        "source": source,
        "product": manifest["product"],
        "mode": MODES[view],
        "rows": manifest["rows"],
        "treated": manifest["treated"],
        "conditions": manifest["conditions"],
        "manifest": {"file": relative(manifest_path, context.root),
                     "sha256": sha256(manifest_path)},
        "files": {path.name: {"sha256": sha256(path), "bytes": path.stat().st_size}
                  for path in installed},
        "evidence": context.evidence,
        "environment": environment(),
        "interpreter": sys.executable,
        "repo_root": str(context.root),
        "staging_removed": not keep_staging,
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "peak_rss_kib": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        "interpreter_warning": INTERPRETER_WARNING,
    }
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--view", required=True, choices=VIEWS)
    parser.add_argument("--source", required=True, choices=SOURCES)
    parser.add_argument("--output", type=Path, default=STEP / "output/views",
                        help="views root; the product lands in <output>/<view>/")
    parser.add_argument("--staging", type=Path, default=STEP / "execution/checkpoints",
                        help="disposable build directory outside output/")
    parser.add_argument("--repo-root", type=Path, default=None,
                        help="materialized checkout owning the pinned payloads")
    parser.add_argument("--record", type=Path, help="write the build receipt JSON here")
    parser.add_argument("--keep-staging", action="store_true",
                        help="keep the staged bytes for failure evidence")
    # Hermetic test hooks; production always uses the released size and non-fixture inputs.
    parser.add_argument("--size", type=int, default=5000, help=argparse.SUPPRESS)
    parser.add_argument("--fixture", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--corpus", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--axis", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--accepted", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--provenance-root", default="repo", help=argparse.SUPPRESS)
    args = parser.parse_args()
    receipt = build_product(view=args.view, source=args.source, output=args.output,
                            staging=args.staging, repo_root_path=args.repo_root,
                            size=args.size, fixture=args.fixture, corpus=args.corpus,
                            axis=args.axis, accepted=args.accepted,
                            provenance_root=args.provenance_root,
                            keep_staging=args.keep_staging)
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(text)
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
