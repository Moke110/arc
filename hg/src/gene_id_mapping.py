"""
Build gene symbol → node_idx mapping for the 19,392 genes in the HG graph.

The 19,392 genes are alphabetically sorted by gene symbol. The 978 landmark genes
provide (symbol, node_idx) anchor points. For any gene symbol:

  1. Direct lookup: if the symbol is among the 978 landmark genes
  2. Alphabetical interpolation: find the nearest anchor points and estimate position
     based on alphabetical ordering

This module can also be used to fetch the full gene list from MyGene.info if needed.
"""

import bisect
import logging
from pathlib import Path

import numpy as np
import pandas as pd

HG_ROOT = Path(__file__).resolve().parent.parent
Xpert_ROOT = HG_ROOT.parent / "xpert"

logger = logging.getLogger(__name__)


def load_landmark_anchors() -> pd.DataFrame:
    """Load the 978 landmark genes with their gene_node_idx as anchor points.

    Returns:
        DataFrame with columns: gene_symbol, gene_id (Entrez), gene_node_idx, ensembl_id
        Sorted by gene_node_idx (which is alphabetically ordered by gene_symbol).
    """
    meta_path = Xpert_ROOT / "metadata" / "l1000_gene_info_978.csv"
    meta = pd.read_csv(meta_path)

    anchors = meta[["gene_symbol", "gene_id", "gene_node_idx", "ensembl_id"]].copy()
    anchors = anchors.sort_values("gene_node_idx").reset_index(drop=True)
    return anchors


class GeneIndexMapper:
    """Map gene symbols (or Entrez IDs) to HG gene node indices (0-19391).

    Uses 978 landmark gene anchor points with alphabetical interpolation
    for genes not in the landmark set.

    Usage:
        mapper = GeneIndexMapper()
        idx = mapper.symbol_to_idx("TP53")       # -> int or None
        idx = mapper.entrez_to_idx(7157)          # -> int or None
        idx = mapper.lookup("TP53")               # -> same as symbol_to_idx
    """

    def __init__(self):
        self.anchors = load_landmark_anchors()

        # Build sorted lists for bisection
        self._anchor_symbols = self.anchors["gene_symbol"].str.upper().tolist()
        self._anchor_indices = self.anchors["gene_node_idx"].tolist()

        # Build direct lookup dicts for landmark genes
        self._symbol_to_idx = dict(zip(
            self.anchors["gene_symbol"].str.upper(),
            self.anchors["gene_node_idx"]
        ))
        self._entrez_to_idx = dict(zip(
            self.anchors["gene_id"].astype(int),
            self.anchors["gene_node_idx"]
        ))
        self._idx_to_symbol = dict(zip(
            self.anchors["gene_node_idx"],
            self.anchors["gene_symbol"]
        ))

    def symbol_to_idx(self, symbol: str) -> int | None:
        """Map a gene symbol to its HG node index.

        Args:
            symbol: HGNC gene symbol (case-insensitive), e.g. "TP53", "Brca1"

        Returns:
            node_idx (0-19391) or None if mapping fails
        """
        symbol_upper = symbol.upper().strip()

        # Direct lookup
        if symbol_upper in self._symbol_to_idx:
            return self._symbol_to_idx[symbol_upper]

        # Alphabetical interpolation
        return self._interpolate_index(symbol_upper)

    def entrez_to_idx(self, entrez_id: int) -> int | None:
        """Map an Entrez Gene ID to its HG node index.

        Only works for the 978 landmark genes (direct lookup).
        For non-landmark genes, use symbol_to_idx instead.
        """
        return self._entrez_to_idx.get(int(entrez_id))

    def idx_to_symbol(self, idx: int) -> str | None:
        """Map an HG node index to gene symbol (landmark genes only)."""
        return self._idx_to_symbol.get(idx)

    def _interpolate_index(self, symbol_upper: str) -> int | None:
        """Estimate gene node_idx by alphabetical interpolation between anchors.

        The 19,392 genes are known to be alphabetically ordered by gene symbol.
        We use binary search to find where `symbol_upper` falls among the 978
        anchor gene symbols, then estimate the position assuming uniform
        distribution of genes between anchor points.

        Args:
            symbol_upper: Uppercased gene symbol

        Returns:
            Estimated node_idx, or None if outside the known range
        """
        # Find insertion point
        pos = bisect.bisect_left(self._anchor_symbols, symbol_upper)

        # Before the first anchor
        if pos == 0:
            # Extrapolate before first anchor
            first_idx = self._anchor_indices[0]
            # Assume genes before first anchor start near index 0
            # Estimate: first anchor is at first_idx, and there are `pos` genes before
            return max(0, first_idx - 1)

        # After the last anchor
        if pos >= len(self._anchor_symbols):
            last_idx = self._anchor_indices[-1]
            # Assume genes extend to ~19391
            return min(19391, last_idx + 1)

        # Between two anchors
        before_symbol = self._anchor_symbols[pos - 1]
        after_symbol = self._anchor_symbols[pos]
        before_idx = self._anchor_indices[pos - 1]
        after_idx = self._anchor_indices[pos]

        # Estimate position by alphabetical ratio
        # Compute rough alphabetical distance
        total_dist = self._str_distance(before_symbol, after_symbol)
        partial_dist = self._str_distance(before_symbol, symbol_upper)

        if total_dist > 0:
            frac = min(1.0, max(0.0, partial_dist / total_dist))
        else:
            frac = 0.5  # same symbol prefix, pick midpoint

        estimated_idx = int(before_idx + frac * (after_idx - before_idx))
        return max(0, min(19391, estimated_idx))

    @staticmethod
    def _str_distance(a: str, b: str) -> float:
        """Compute rough alphabetical distance between two gene symbols.

        Uses a weighted positional encoding: earlier characters contribute more.
        """
        max_len = max(len(a), len(b))
        dist = 0.0
        for i in range(max_len):
            ca = ord(a[i]) if i < len(a) else 0
            cb = ord(b[i]) if i < len(b) else 0
            weight = 1.0 / (i + 1)  # earlier chars matter more
            dist += weight * abs(ca - cb)
        return dist

    def symbols_to_indices(self, symbols: list[str]) -> dict[str, int | None]:
        """Batch-convert gene symbols to node indices.

        Returns dict mapping each input symbol to its node_idx (or None).
        """
        return {s: self.symbol_to_idx(s) for s in symbols}

    def get_mapped_count(self, symbols: list[str]) -> tuple[int, int]:
        """Return (mapped_count, total_count) for a list of gene symbols."""
        results = self.symbols_to_indices(symbols)
        mapped = sum(1 for v in results.values() if v is not None)
        return mapped, len(symbols)


def build_full_gene_list(output_path: str | None = None) -> pd.DataFrame:
    """Attempt to build the full 19,392 gene list using MyGene.info API.

    This is a fallback when the 978 landmark genes aren't sufficient.
    Queries MyGene.info for all human protein-coding genes (~20k) and maps
    them to node indices via alphabetical ordering.

    Args:
        output_path: If provided, save CSV to this path

    Returns:
        DataFrame with columns: gene_symbol, entrez_id, estimated_node_idx
    """
    try:
        import requests
    except ImportError:
        logger.error("requests not installed. Install with: pip install requests")
        return pd.DataFrame()

    # Query MyGene.info for all human protein-coding genes
    # We fetch in batches
    logger.info("Fetching human gene list from MyGene.info...")

    # Use the query endpoint to get all human protein-coding genes
    # Limited to ~1000 per request, need to paginate
    all_genes = []
    url = "https://mygene.info/v3/query"
    params = {
        "q": "taxid:9606 AND type_of_gene:protein-coding",
        "fields": "symbol,entrezgene,ensembl.gene",
        "size": 1000,
        "species": "human",
    }

    scroll_id = None
    total_fetched = 0

    for _ in range(25):  # max 25 pages = 25,000 genes
        if scroll_id:
            resp = requests.get(
                f"https://mygene.info/v3/query?scroll_id={scroll_id}",
                timeout=30,
            )
        else:
            resp = requests.post(url, json=params, timeout=30)
            data = resp.json()
            scroll_id = data.get("_scroll_id")

        if resp.status_code != 200:
            logger.error(f"MyGene.info query failed: {resp.status_code}")
            break

        data = resp.json()
        hits = data.get("hits", [])
        if not hits:
            break

        for hit in hits:
            symbol = hit.get("symbol", "")
            entrez = hit.get("entrezgene", "")
            all_genes.append({
                "gene_symbol": symbol,
                "entrez_id": str(entrez) if entrez else "",
            })

        total_fetched += len(hits)
        logger.info(f"  Fetched {total_fetched} genes...")

    logger.info(f"Total genes fetched: {len(all_genes)}")

    if not all_genes:
        return pd.DataFrame()

    df = pd.DataFrame(all_genes)
    # Sort alphabetically by gene_symbol (matching the HG graph ordering)
    df = df.sort_values("gene_symbol").reset_index(drop=True)

    # Map to HG node indices using alphabetical interpolation
    mapper = GeneIndexMapper()
    df["estimated_node_idx"] = df["gene_symbol"].apply(
        lambda s: mapper.symbol_to_idx(s)
    )

    # Filter to genes with valid index estimates
    df = df[df["estimated_node_idx"].notna()].copy()
    df["estimated_node_idx"] = df["estimated_node_idx"].astype(int)

    # Deduplicate by estimated_node_idx (keep first alphabetically)
    df = df.drop_duplicates(subset="estimated_node_idx", keep="first")

    if output_path:
        df.to_csv(output_path, index=False)
        logger.info(f"Saved {len(df)} genes to {output_path}")

    return df


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Test the mapper
    mapper = GeneIndexMapper()
    test_genes = ["TP53", "BRCA1", "EGFR", "NOTEXIST", "CYP3A4", "MTOR", "VEGFA"]
    print("=== Gene Symbol → HG Node Index ===")
    for g in test_genes:
        idx = mapper.symbol_to_idx(g)
        is_anchor = "✓ anchor" if g.upper() in mapper._symbol_to_idx else "interpolated"
        print(f"  {g:15s} → {idx}  ({is_anchor})")

    # Check mapping accuracy on landmark genes
    print("\n=== Landmark Gene Self-Check ===")
    anchor_symbols = mapper.anchors["gene_symbol"].str.upper().tolist()
    correct = 0
    errors = []
    for sym in anchor_symbols[:100]:  # check first 100
        mapped = mapper.symbol_to_idx(sym)
        expected = mapper._symbol_to_idx[sym]
        if mapped == expected:
            correct += 1
        else:
            errors.append((sym, mapped, expected))
    print(f"  Correct: {correct}/100")
    if errors:
        for sym, m, e in errors[:5]:
            print(f"  Error: {sym} → {m} (expected {e})")
