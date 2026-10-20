"""Gene identity: source identifiers -> one current HGNC symbol.

Rules (ADR 0026 decision 1)
---------------------------
* The universe is the four-source intersection, named by the **latest HGNC
  symbol**.
* ``dtoxs_primary`` / ``dtoxs_ipsc`` are symbol spaces; the DToxS package's
  own ``UCSC-hg38-Gene-HGNC-Entrez-Map`` supplies an Entrez ID used to verify
  (or correct) a symbol that HGNC has since reassigned.
* ``gse244740`` is an Ensembl space (version suffix stripped).
* ``gse262419`` probes are ``SYMBOL_probeid`` rows.  The symbol is the identity
  key (duplicate probes collapse by symbol); the numeric suffix is a probe id,
  used only as an Entrez fallback when the symbol does not resolve.  It is not
  treated as an Entrez ID (e.g. ``A1BG_25586``).
* Unresolvable or ambiguous identifiers are dropped and counted, never guessed.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .annotation import HgncIndex

MATCH_TYPES = ["direct", "alias", "previous", "entrez", "ensembl",
               "unmapped", "ambiguous"]


@dataclass
class SourceGeneMap:
    """Resolution of one source's raw gene ids."""

    dataset_id: str
    raw_ids: list[str]
    hgnc_ids: np.ndarray            # object array, None where unresolved
    match_types: np.ndarray         # object array of MATCH_TYPES
    verified: np.ndarray            # bool: source Entrez agreed with HGNC
    counts: Counter = field(default_factory=Counter)

    @property
    def resolved_ids(self) -> set[str]:
        return {h for h in self.hgnc_ids if h is not None}

    def targets(self, universe: "GeneUniverse") -> np.ndarray:
        """Position of each raw gene in the universe, or -1."""
        return np.array([universe.position.get(h, -1) if h is not None else -1
                         for h in self.hgnc_ids], dtype=np.int32)


class GeneUniverse:
    """Fixed-order gene axis (one row per current HGNC symbol)."""

    def __init__(self, genes: pd.DataFrame):
        self.genes = genes.reset_index(drop=True)
        self.symbols = self.genes["symbol"].astype(str).tolist()
        self.position = {h: i for i, h in enumerate(self.genes["hgnc_id"])}

    @property
    def var(self) -> pd.DataFrame:
        var = self.genes[["symbol", "entrez_id", "ensembl_id",
                          "match_type"]].set_index("symbol")
        return var

    def __len__(self) -> int:
        return len(self.genes)


def _resolve_symbol_space(index: HgncIndex, raw_id: str,
                          entrez: str | None) -> tuple[str | None, str, bool]:
    """Resolve a symbol-space identifier, Entrez-first when available."""
    if entrez:
        hgnc_id, _ = index.resolve_entrez(entrez)
        if hgnc_id is not None:
            ref = index.ref(hgnc_id)
            return hgnc_id, ("direct" if ref.symbol == raw_id else "entrez"), True
    hgnc_id, match = index.resolve_symbol(raw_id)
    return hgnc_id, match, False


def map_symbol_space(dataset_id: str, raw_ids: list[str], index: HgncIndex,
                     entrez_map: dict[str, str] | None = None) -> SourceGeneMap:
    hgnc_ids, match_types, verified = [], [], []
    for raw_id in raw_ids:
        entrez = (entrez_map or {}).get(raw_id, "")
        hgnc_id, match, ok = _resolve_symbol_space(index, raw_id, entrez)
        hgnc_ids.append(hgnc_id)
        match_types.append(match)
        verified.append(ok)
    mapping = SourceGeneMap(dataset_id, raw_ids,
                            np.array(hgnc_ids, dtype=object),
                            np.array(match_types, dtype=object),
                            np.array(verified, dtype=bool))
    mapping.counts = Counter(match_types)
    return mapping


def map_ensembl_space(dataset_id: str, raw_ids: list[str],
                      index: HgncIndex) -> SourceGeneMap:
    hgnc_ids, match_types = [], []
    for raw_id in raw_ids:
        hgnc_id, _ = index.resolve_ensembl(raw_id)
        hgnc_ids.append(hgnc_id)
        match_types.append("ensembl" if hgnc_id is not None else "unmapped")
    mapping = SourceGeneMap(dataset_id, raw_ids,
                            np.array(hgnc_ids, dtype=object),
                            np.array(match_types, dtype=object),
                            np.zeros(len(raw_ids), dtype=bool))
    mapping.counts = Counter(match_types)
    return mapping


def split_probe(probe: str) -> tuple[str, str]:
    """``A1BG_25586`` -> (``A1BG``, ``25586``); symbol may contain underscores."""
    parts = str(probe).rsplit("_", 1)
    return (parts[0], parts[1]) if len(parts) == 2 else (str(probe), "")


def map_probe_space(dataset_id: str, raw_ids: list[str],
                    index: HgncIndex) -> SourceGeneMap:
    hgnc_ids, match_types, verified = [], [], []
    for probe in raw_ids:
        symbol, suffix = split_probe(probe)
        hgnc_id, match = index.resolve_symbol(symbol)
        ok = False
        if hgnc_id is not None:
            # The suffix is a probe id, not an Entrez ID.  It may coincidentally
            # equal the gene's Entrez; record agreement but never use it to
            # override a resolved symbol.
            if suffix.isdigit():
                by_entrez, _ = index.resolve_entrez(suffix)
                ok = by_entrez == hgnc_id
        elif suffix.isdigit():
            hgnc_id, _ = index.resolve_entrez(suffix)
            match = "entrez" if hgnc_id is not None else "unmapped"
        hgnc_ids.append(hgnc_id)
        match_types.append(match)
        verified.append(ok)
    mapping = SourceGeneMap(dataset_id, raw_ids,
                            np.array(hgnc_ids, dtype=object),
                            np.array(match_types, dtype=object),
                            np.array(verified, dtype=bool))
    mapping.counts = Counter(match_types)
    return mapping


def _aggregate_match_type(types: set[str]) -> str:
    """Strongest evidence needed to recover a gene in any source."""
    if "ambiguous" in types:
        return "ambiguous"
    if "alias" in types or "previous" in types:
        return "alias"
    if "entrez" in types:
        return "entrez"
    return "direct"


def compute_universe(index: HgncIndex,
                     mappings: list[SourceGeneMap]) -> GeneUniverse:
    """Four-source intersection in fixed order (HGNC symbol, lexicographic)."""
    resolved = [mapping.resolved_ids for mapping in mappings]
    shared = set.intersection(*resolved) if resolved else set()
    rows = []
    for hgnc_id in shared:
        ref = index.ref(hgnc_id)
        types = set()
        for mapping in mappings:
            positions = np.flatnonzero(mapping.hgnc_ids == hgnc_id)
            types.update(mapping.match_types[positions])
        rows.append({
            "symbol": ref.symbol,
            "hgnc_id": ref.hgnc_id,
            "entrez_id": ref.entrez_id,
            "ensembl_id": ref.ensembl_id,
            "match_type": _aggregate_match_type(types),
        })
    genes = pd.DataFrame(rows).sort_values("symbol", kind="stable").reset_index(drop=True)
    return GeneUniverse(genes)


def universe_from_genes_csv(path) -> GeneUniverse:
    genes = pd.read_csv(path, dtype={"entrez_id": str, "ensembl_id": str})
    required = {"symbol", "hgnc_id", "entrez_id", "ensembl_id", "match_type"}
    missing = required - set(genes.columns)
    if missing:
        raise ValueError(f"frozen genes.csv missing columns: {sorted(missing)}")
    if genes["symbol"].duplicated().any() or genes["symbol"].isna().any():
        raise ValueError("frozen genes.csv has duplicate or empty symbols")
    if genes["hgnc_id"].duplicated().any():
        raise ValueError("frozen genes.csv has duplicate HGNC ids")
    return GeneUniverse(genes)
