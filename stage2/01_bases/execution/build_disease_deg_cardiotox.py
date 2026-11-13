#!/usr/bin/env python3
"""Publish or replay the cardiotoxicity class DEG products (ticket #189).

``output/disease-deg/`` gains three ticket-owned products:

* ``cardiotox/`` — the pinned code->drug->mechanism crosswalk, the class
  selection/membership/excluded-code evidence and the per-source eligibility
  record, with its manifest;
* ``cardiotox_TOP2/`` — one pooled class DEG table and within-dataset rank per
  source (GSE243674, DToxS-iPSC/GSE217421, GSE244740), ``catalog.json`` and
  ``manifest.json``;
* ``cardiotox_other/`` — one pooled class DEG table and within-dataset rank per
  entered DToxS non-TOP2 class, ``catalog.json`` and ``manifest.json``.

Every row is produced by the single public implementation in
``output/api/bases/cardiotox.py``; this driver only binds paths and publishes.

    python stage2/01_bases/execution/build_disease_deg_cardiotox.py publish
    python stage2/01_bases/execution/build_disease_deg_cardiotox.py verify
"""
import argparse
import json
import sys
from pathlib import Path

import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import cardiotox

STEP = Path(__file__).resolve().parents[1]
REPO = next(path for path in Path(__file__).resolve().parents
            if (path / ".git").exists())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("publish", "verify"))
    parser.add_argument("--staging", type=Path,
                        default=STEP / "execution/checkpoints")
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    if args.action == "publish":
        manifests, resources = cardiotox.publish_release(REPO, args.staging)
    else:
        manifests, resources = cardiotox.replay_release(REPO, args.staging)
    result = {
        "action": args.action,
        "output_root": "stage2/01_bases/output/disease-deg",
        "products": {
            name: {
                "manifest_sha256": cardiotox.sha256(
                    STEP / "output/disease-deg" / name / "manifest.json"),
                "counts": manifest.get("counts"),
            }
            for name, manifest in (manifests.items()
                                   if args.action == "publish" else {}.items())
        },
        "resources": resources,
        "replay": manifests if args.action == "verify" else None,
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(text)
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
