"""Shared identity-table helper for the Stage 2 corpus tables (ticket #16).

The retired v1 corpus writer (``freeze_corpus`` / ``load_frozen_corpus`` /
``check_drift``) lived here until the 02/04 merge left it with no caller.  The
frozen identity tables now live read-only in
``stage2/_frozen/cm-data-v1/corpus/v1/`` and ADR 0033 retired the writer; what
remains is the one helper the v2 table build still uses.
"""

from __future__ import annotations

import pandas as pd


def source_gene_set_frame(index, source_maps: dict) -> pd.DataFrame:
    """One row per (source, resolved HGNC gene).

    Produces ``output/corpus/v2/tables/source_gene_sets.csv``: the same content
    the frozen ``corpus/v1/source_gene_sets.csv`` recorded (104,149 rows: the
    sum of the four native axes).
    """
    source_rows = []
    for dataset_id, mapping in source_maps.items():
        seen = set()
        for hgnc_id, match, verified in zip(mapping.hgnc_ids, mapping.match_types,
                                            mapping.verified):
            if hgnc_id is None or hgnc_id in seen:
                continue
            seen.add(hgnc_id)
            source_rows.append({
                "dataset_id": dataset_id,
                "hgnc_id": hgnc_id,
                "symbol": index.ref(hgnc_id).symbol,
                "match_type": match,
                "entrez_verified": bool(verified),
            })
    return pd.DataFrame(source_rows).sort_values(
        ["dataset_id", "symbol"], kind="stable").reset_index(drop=True)
