"""Regenerate the tracked ticket-03 discarded-member evidence table.

Ticket 03's item 4 asks each member-choice product (``solo``, ``solo-deepest``) to record
how many members its per-condition choice discarded. The payload's ``members.csv`` is
gitignored with the rest of the payload, so the counts were only derivable. This module
turns the published ``members.csv`` companion of each member-choice product into one
reviewable, tracked table::

    python evidence.py --views output/views --output log/ticket-03-discarded-members.csv

A row is one condition of one member-choice product:

* ``source`` -- ``tenv`` or ``tahoe_bulk``;
* ``profile_id`` -- the published product row (the chosen member's own id);
* ``chosen_member_id`` -- the member the rule kept (equal to ``profile_id`` by design,
  recorded separately so the equality is explicit in the evidence);
* ``member_count`` -- the condition's released member count, from the ``members.csv``
  group that the product emitted;
* ``discarded_count`` -- ``member_count - 1`` (the members the choice dropped);
* ``selection_rule`` -- ``draw`` for ``solo``, ``deepest`` for ``solo-deepest``.

The module reads only companion tables (never a payload matrix); the test recomputes the
same counts from the independent ``merge`` member table and reconciles them.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import common  # noqa: E402

COLUMNS = ("source", "profile_id", "chosen_member_id", "member_count", "discarded_count",
           "selection_rule")
#: The member-choice products, in publication order, with their selection rule.
MEMBERS = (("solo", "draw"), ("solo-deepest", "deepest"))


def _members(path):
    return pd.read_csv(path, dtype={"profile_id": str, "member_profile_id": str})


def discarded_members(views_root):
    """One row per condition of every published member-choice product.

    Fails closed when a product has no ``members.csv``: a member-choice view without its
    member table cannot support the discarded-count claim.
    """
    views_root = Path(views_root)
    rows = []
    for view, rule in MEMBERS:
        for source in common.SOURCE_ORDER:
            path = views_root / view / f"{source}.members.csv"
            if not path.is_file():
                if view == "solo-deepest" and source == "tahoe_bulk":
                    continue  # the sensitivity arm is tenv-only by design
                raise ValueError(f"{view}/{source}: the member table is missing: {path}")
            table = _members(path)
            grouped = table.groupby("profile_id", sort=True)
            for identity, block in grouped:
                members = [str(member) for member in block["member_profile_id"]]
                count = len(members)
                # The harness pairs the chosen identity with every member of its
                # condition, so the chosen member must be one of the block's members.
                if str(identity) not in members:
                    raise ValueError(
                        f"{view}/{source}: the row identity {identity} is not a member of "
                        "its own condition")
                rows.append({
                    "source": source,
                    "profile_id": str(identity),
                    "chosen_member_id": str(identity),
                    "member_count": count,
                    "discarded_count": count - 1,
                    "selection_rule": rule,
                })
    frame = pd.DataFrame(rows, columns=list(COLUMNS))
    frame["_source"] = frame["source"].map(
        {source: index for index, source in enumerate(common.SOURCE_ORDER)})
    frame["_rule"] = frame["selection_rule"].map(
        {rule: index for index, (_, rule) in enumerate(MEMBERS)})
    return frame.sort_values(["_source", "_rule", "profile_id"], kind="mergesort").drop(
        columns=["_source", "_rule"]).reset_index(drop=True)


def write_evidence(views_root, output):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame = discarded_members(views_root)
    frame.to_csv(output, index=False, lineterminator="\n")
    return output, frame


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--views", type=Path, required=True,
                        help="the published views root (carries <view>/<source>.members.csv)")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args(argv)
    path, frame = write_evidence(arguments.views, arguments.output)
    counts = frame.groupby(["source", "selection_rule"]).size().to_dict()
    print(f"{path}: {len(frame)} rows {counts}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
