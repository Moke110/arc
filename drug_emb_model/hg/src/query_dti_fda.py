"""
Query drug-target interaction (DTI) data for FDA drugs from ChEMBL API.

Two-step approach:
  1. Query activity endpoint → get target_chembl_id for each drug's protein targets
  2. Resolve target_chembl_id → gene symbol via target endpoint (with caching)

Also supports batch mode and resumption via local cache files.

Output: numpy dict with edge_index for FDA drug→gene connections.
  drug indices are LOCAL (0..N_fda-1), gene indices are HG global (0-19391).
"""

import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HG_ROOT = Path(__file__).resolve().parent.parent
Project_ROOT = HG_ROOT.parent

sys.path.insert(0, str(HG_ROOT))
from src.gene_id_mapping import GeneIndexMapper

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"
CACHE_DIR = HG_ROOT / "data" / "dti_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# Global target→gene_symbol cache
TARGET_CACHE_FILE = CACHE_DIR / "target_gene_symbols.json"


# ---------------------------------------------------------------------------
# API helpers
# ---------------------------------------------------------------------------
def _chembl_request(endpoint: str, params: dict | None = None,
                    timeout: int = 30, max_retries: int = 3) -> dict | None:
    """Make a ChEMBL API request with retry and rate-limit handling."""
    import requests
    url = f"{CHEMBL_API_BASE}/{endpoint}"
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, params=params, timeout=timeout)
            if resp.status_code == 200:
                return resp.json()
            elif resp.status_code == 404:
                return None
            elif resp.status_code == 429:
                wait = 2 ** attempt
                logger.warning(f"Rate limited, waiting {wait}s...")
                time.sleep(wait)
            else:
                logger.warning(f"API error {resp.status_code}: {url}")
                time.sleep(1)
        except requests.RequestException as e:
            logger.warning(f"Request failed (attempt {attempt+1}): {e}")
            time.sleep(2 ** attempt)
    return None


# ---------------------------------------------------------------------------
# Target → Gene Symbol resolution
# ---------------------------------------------------------------------------
def _load_target_cache() -> dict[str, str | None]:
    """Load the target_chembl_id → gene_symbol cache from disk."""
    if TARGET_CACHE_FILE.exists():
        with open(TARGET_CACHE_FILE) as f:
            return json.load(f)
    return {}


def _save_target_cache(cache: dict):
    """Save the target_chembl_id → gene_symbol cache."""
    with open(TARGET_CACHE_FILE, "w") as f:
        json.dump(cache, f, indent=2)


def resolve_target_to_gene(target_chembl_id: str,
                           cache: dict[str, str | None]) -> str | None:
    """Resolve a ChEMBL target ID to a gene symbol.

    Queries the ChEMBL target endpoint to extract the GENE_SYMBOL synonym
    from target_component_synonyms. Results are cached persistently.

    Args:
        target_chembl_id: e.g. "CHEMBL319"
        cache: mutable dict, updated in-place

    Returns:
        Gene symbol (e.g. "ADRA1A") or None
    """
    if target_chembl_id in cache:
        return cache[target_chembl_id]

    data = _chembl_request("target.json",
                           params={"target_chembl_id": target_chembl_id, "limit": 1})
    if not data:
        cache[target_chembl_id] = None
        return None

    targets = data.get("targets", [])
    if not targets:
        cache[target_chembl_id] = None
        return None

    # Search target_component_synonyms for GENE_SYMBOL type
    for tc in targets[0].get("target_components", []):
        for syn in tc.get("target_component_synonyms", []):
            if syn.get("syn_type") == "GENE_SYMBOL":
                gene_symbol = syn.get("component_synonym", "")
                if gene_symbol:
                    cache[target_chembl_id] = gene_symbol
                    return gene_symbol

    cache[target_chembl_id] = None
    return None


def batch_resolve_targets(target_ids: list[str]) -> dict[str, str | None]:
    """Resolve a batch of target_chembl_ids to gene symbols.

    Uses persistent cache; only queries API for uncached targets.
    Saves cache after batch resolution.
    """
    cache = _load_target_cache()
    logger.info(f"Target cache: {len(cache)} entries loaded, "
                f"{sum(1 for v in cache.values() if v is not None)} resolved")

    uncached = [t for t in target_ids if t not in cache]
    if uncached:
        logger.info(f"Resolving {len(uncached)} new target IDs...")
        for i, tid in enumerate(uncached):
            resolve_target_to_gene(tid, cache)
            if (i + 1) % 100 == 0:
                logger.info(f"  Resolved {i+1}/{len(uncached)} targets")
                _save_target_cache(cache)
            if i % 10 == 0 and i > 0:
                time.sleep(0.15)  # rate limiting
        _save_target_cache(cache)
        logger.info(f"Target cache updated: {len(cache)} entries")

    return {tid: cache.get(tid) for tid in target_ids}


# ---------------------------------------------------------------------------
# Drug target query
# ---------------------------------------------------------------------------
def query_drug_activities(chembl_id: str,
                          min_pchembl: float = 5.0) -> list[dict]:
    """Get all protein target activities for a drug from ChEMBL.

    Returns unique (target_chembl_id, pchembl_value, standard_type) tuples.

    pChEMBL thresholds:
      5.0 ≈ 10µM (broad)
      6.0 ≈ 1µM  (standard)
      7.0 ≈ 100nM (high affinity)
    """
    cache_file = CACHE_DIR / f"{chembl_id}_activities.json"
    if cache_file.exists():
        with open(cache_file) as f:
            return json.load(f)

    data = _chembl_request("activity.json", params={
        "molecule_chembl_id": chembl_id,
        "standard_type__in": "IC50,Ki,Kd,EC50,AC50,Potency,Inhibition,Activity",
        "pchembl_value__gte": min_pchembl,
        "target_type": "SINGLE PROTEIN",
        "limit": 500,
    })

    if not data:
        return []

    targets = {}
    for act in data.get("activities", []):
        tid = act.get("target_chembl_id", "")
        if not tid:
            continue
        pchembl = act.get("pchembl_value") or 0
        if tid not in targets or pchembl > targets[tid]["pchembl_value"]:
            targets[tid] = {
                "target_chembl_id": tid,
                "target_pref_name": act.get("target_pref_name", ""),
                "standard_type": act.get("standard_type", ""),
                "pchembl_value": pchembl,
                "standard_value": act.get("standard_value"),
                "standard_units": act.get("standard_units"),
            }

    result = sorted(targets.values(), key=lambda x: x["pchembl_value"], reverse=True)

    with open(cache_file, "w") as f:
        json.dump(result, f, indent=2)

    return result


# ---------------------------------------------------------------------------
# Main DTI collection
# ---------------------------------------------------------------------------
def collect_fda_dti(fda_csv_path: str,
                    output_path: str | None = None,
                    min_pchembl: float = 5.0,
                    start_idx: int = 0,
                    end_idx: int | None = None,
                    resume: bool = True) -> dict:
    """Collect DTI data for FDA drugs from ChEMBL.

    Two-phase approach:
      Phase 1: Query all drug activities, collect unique target_chembl_ids
      Phase 2: Batch-resolve target_chembl_ids → gene symbols
      Phase 3: Map gene symbols → HG node indices, build edges

    Args:
        fda_csv_path: path to fda_drugs.csv
        output_path: save edge data as .npy
        min_pchembl: minimum pChEMBL (default 5.0 = 10µM)
        start_idx, end_idx: slice of FDA drug list to process
        resume: skip drugs with existing activity cache

    Returns:
        dict with edge_index, drug_fda_local_indices, gene_hg_indices, stats
    """
    df = pd.read_csv(fda_csv_path)
    if end_idx is None:
        end_idx = len(df)
    df = df.iloc[start_idx:end_idx].copy()
    n_drugs = len(df)
    logger.info(f"Processing {n_drugs} FDA drugs (indices {start_idx}-{end_idx-1})")
    logger.info(f"pChEMBL threshold: {min_pchembl}")

    # ------------------------------------------------------------------
    # Phase 1: Query activities for all drugs
    # ------------------------------------------------------------------
    drug_targets_list = []  # list[list[dict]] — per drug
    all_target_ids = set()

    for i, (_, row) in enumerate(df.iterrows()):
        cid = row["chembl_id"]
        name = row.get("pref_name", cid)

        if (i + 1) % 100 == 0:
            logger.info(f"  Phase 1: {i+1}/{n_drugs} drugs queried, "
                        f"{len(all_target_ids)} unique target IDs found")

        cache_file = CACHE_DIR / f"{cid}_activities.json"
        if resume and cache_file.exists():
            with open(cache_file) as f:
                activities = json.load(f)
        else:
            activities = query_drug_activities(cid, min_pchembl=min_pchembl)

        drug_targets_list.append(activities)
        for act in activities:
            all_target_ids.add(act["target_chembl_id"])

        if i % 5 == 0 and i > 0:
            time.sleep(0.15)  # rate limiting

    logger.info(f"Phase 1 done: {len(all_target_ids)} unique target IDs across {n_drugs} drugs")

    # ------------------------------------------------------------------
    # Phase 2: Resolve target IDs → gene symbols
    # ------------------------------------------------------------------
    target_id_list = sorted(all_target_ids)
    target_to_gene = batch_resolve_targets(target_id_list)

    n_resolved = sum(1 for v in target_to_gene.values() if v is not None)
    logger.info(f"Phase 2 done: {n_resolved}/{len(target_id_list)} targets resolved to gene symbols")

    # ------------------------------------------------------------------
    # Phase 3: Map gene symbols → HG indices, build edges
    # ------------------------------------------------------------------
    mapper = GeneIndexMapper()
    all_edges = []
    drug_stats = []

    for i, activities in enumerate(drug_targets_list):
        row = df.iloc[i]
        drug_targets = set()

        for act in activities:
            gene_sym = target_to_gene.get(act["target_chembl_id"])
            if gene_sym:
                hg_idx = mapper.symbol_to_idx(gene_sym)
                if hg_idx is not None:
                    pval = act.get("pchembl_value", 0) or 0
                    drug_targets.add((
                        gene_sym.upper(),
                        hg_idx,
                        f"chembl:{act['standard_type']}:p{float(pval):.1f}",
                    ))

        for gs, gene_idx, source in drug_targets:
            all_edges.append({
                "fda_local_idx": i,
                "gene_hg_idx": gene_idx,
                "gene_symbol": gs,
                "source": source,
            })

        drug_stats.append({
            "fda_local_idx": i,
            "chembl_id": row["chembl_id"],
            "drug_name": row.get("pref_name", ""),
            "n_targets": len(drug_targets),
        })

    # ------------------------------------------------------------------
    # Build result arrays
    # ------------------------------------------------------------------
    n_edges = len(all_edges)
    edge_index = np.zeros((2, max(n_edges, 1)), dtype=np.int64)
    gene_indices = np.zeros(max(n_edges, 1), dtype=np.int64)
    sources = []

    for i, e in enumerate(all_edges):
        edge_index[0, i] = e["fda_local_idx"]
        edge_index[1, i] = e["gene_hg_idx"]
        gene_indices[i] = e["gene_hg_idx"]
        sources.append(e["source"])

    if n_edges == 0:
        edge_index = np.zeros((2, 0), dtype=np.int64)
        gene_indices = np.zeros(0, dtype=np.int64)

    drugs_with_targets = sum(1 for s in drug_stats if s["n_targets"] > 0)
    unique_genes = len(set(gene_indices[gene_indices >= 0])) if n_edges > 0 else 0
    anchor_genes = sum(1 for g in gene_indices if g >= 0 and mapper.idx_to_symbol(int(g)) is not None)

    result = {
        "edge_index": edge_index,
        "drug_fda_local_indices": edge_index[0].copy(),
        "gene_hg_indices": gene_indices.copy(),
        "edge_sources": sources,
        "drug_stats": drug_stats,
        "stats": {
            "n_fda_drugs": n_drugs,
            "n_dti_edges": n_edges,
            "drugs_with_targets": drugs_with_targets,
            "drugs_without_targets": n_drugs - drugs_with_targets,
            "unique_target_genes": unique_genes,
            "targets_are_landmark": anchor_genes,
        },
    }

    logger.info(f"\n{'='*50}")
    logger.info(f"Results: {n_edges} DTI edges for {drugs_with_targets}/{n_drugs} drugs")
    logger.info(f"  Unique target genes: {unique_genes} "
                f"({anchor_genes} landmark, {unique_genes - anchor_genes} interpolated)")
    if n_drugs > 0:
        logger.info(f"  Avg targets per drug (all): {n_edges/n_drugs:.1f}")
    if drugs_with_targets > 0:
        logger.info(f"  Avg targets per drug (with targets): {n_edges/drugs_with_targets:.1f}")

    # Top targeted genes
    if n_edges > 0:
        from collections import Counter
        gene_counts = Counter(e["gene_symbol"] for e in all_edges)
        logger.info(f"  Top 10 target genes: {gene_counts.most_common(10)}")

    if output_path:
        save_dict = {
            "edge_index": edge_index,
            "gene_hg_indices": gene_indices,
            "edge_sources": np.array(sources, dtype=object),
            "stats": result["stats"],
        }
        np.save(output_path, save_dict, allow_pickle=True)
        logger.info(f"\nSaved DTI edges to {output_path}")

        stats_df = pd.DataFrame(drug_stats)
        stats_path = Path(output_path).with_suffix(".csv")
        stats_df.to_csv(stats_path, index=False)
        logger.info(f"Saved per-drug stats to {stats_path}")

    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    import argparse
    ap = argparse.ArgumentParser(description="Query DTI data for FDA drugs from ChEMBL")
    ap.add_argument("--fda-csv", default=str(Project_ROOT / "xpert/predict/fda-cm/fda_drugs.csv"))
    ap.add_argument("--output", default=str(HG_ROOT / "data/DTI_edge_index_fda_drugs.npy"))
    ap.add_argument("--min-pchembl", type=float, default=5.0,
                    help="Minimum pChEMBL (5.0=10µM, 6.0=1µM, 7.0=100nM)")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=None)
    ap.add_argument("--test", type=int, default=0,
                    help="Test with first N drugs")
    ap.add_argument("--no-resume", action="store_true",
                    help="Re-query even if cache exists")
    args = ap.parse_args()

    if args.test > 0:
        collect_fda_dti(
            args.fda_csv, output_path=None,
            min_pchembl=args.min_pchembl,
            start_idx=0, end_idx=args.test,
            resume=not args.no_resume,
        )
    else:
        collect_fda_dti(
            args.fda_csv, output_path=args.output,
            min_pchembl=args.min_pchembl,
            start_idx=args.start, end_idx=args.end,
            resume=not args.no_resume,
        )
