#!/usr/bin/env python3
"""
Rebuild MIRI pathway gene sets for drug efficacy scoring.

Pipeline:
  1. Fetch KEGG/GO pathway gene sets from Enrichr (cache to .claude/)
  2. Run preranked GSEA on GEO DEGs against KEGG/GO gene sets
  3. Run gpt-researcher for MIRI literature pathway survey
  4. Combine GEO + literature pathways, intersect with L1000 landmark genes
  5. Apply dual filters: n_genes > 5  AND  coverage >= 978/20000
  6. Merge identical gene sets (KEGG+GO with same genes+direction → KEGG&GO-)
  7. Select top ~30 UP, ~30 DOWN (literature-supported pathways prioritized)
  8. Output MIRI_pathways.json

Usage:
    python scripts/rebuild_pathways.py [--target-up 30] [--target-down 30]
                                      [--skip-literature] [--refresh-cache]

Author: ARC project, 2026-06-11
"""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

# ── Paths ───────────────────────────────────────────────────────────
PROJECT_DIR = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = PROJECT_DIR / "scripts"
GEO_DIR = PROJECT_DIR / "geo"
REF_DIR = PROJECT_DIR / "reference-genesets"
CACHE_DIR = PROJECT_DIR / ".claude"

CACHE_DIR.mkdir(parents=True, exist_ok=True)
REF_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(SCRIPTS_DIR))
from utils import (
    enrichment_score,
    load_l1000_genes,
    passes_coverage_filter,
)

L1000_CSV = "/home/chang/arc/xpert/metadata/l1000_gene_info_978.csv"

# ── Constants ────────────────────────────────────────────────────────
TOTAL_HUMAN_GENES = 20000
N_LANDMARK = 978
COVERAGE_FACTOR = 1.0
MIN_L1000_GENES = 6  # >5
COVERAGE_THRESHOLD = (N_LANDMARK / TOTAL_HUMAN_GENES) * COVERAGE_FACTOR


# ══════════════════════════════════════════════════════════════════════
# 1. Enrichr API — KEGG & GO gene set fetching
# ══════════════════════════════════════════════════════════════════════

def fetch_enrichr_library(library_name, cache_dir=CACHE_DIR):
    """
    Fetch a gene set library from Enrichr, with caching.

    Parameters
    ----------
    library_name : str
        e.g. 'KEGG_2021_Human', 'GO_Biological_Process_2023'

    Returns
    -------
    dict[str, list[str]]
        {pathway_name: [gene_symbol, ...]}
    """
    import urllib.request
    import urllib.error

    cache_file = cache_dir / f"{library_name}.json"

    # Load from cache if available
    if cache_file.exists():
        with open(cache_file) as f:
            cached = json.load(f)
        age_hours = (time.time() - cached.get("_fetched_at", 0)) / 3600
        if age_hours < 168:  # 1 week cache
            terms = cached.get("terms", {})
            print(f"  Loaded {len(terms)} pathways from cache ({cache_file.name}, "
                  f"{age_hours:.0f}h old)")
            return terms

    # Fetch from Enrichr
    url = ("https://maayanlab.cloud/Enrichr/geneSetLibrary"
           f"?mode=json&libraryName={library_name}")
    print(f"  Fetching {library_name} from Enrichr...")

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "ARC-MIRI/1.0"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read().decode())
    except urllib.error.URLError as e:
        print(f"  WARNING: Enrichr fetch failed: {e}")
        if cache_file.exists():
            print(f"  Falling back to cached data")
            with open(cache_file) as f:
                cached = json.load(f)
            return cached.get("terms", {})
        raise

    lib = data.get(library_name, {})
    terms = lib.get("terms", {})

    # Convert {gene: weight} dicts to gene lists
    result = {}
    for name, gene_dict in terms.items():
        result[name] = sorted(gene_dict.keys())

    # Cache
    cache_data = {
        "_fetched_at": time.time(),
        "_library_name": library_name,
        "_n_terms": len(result),
        "terms": result,
    }
    with open(cache_file, 'w') as f:
        json.dump(cache_data, f, indent=2, ensure_ascii=False)
    print(f"  Fetched {len(result)} pathways, cached to {cache_file.name}")

    return result


# ══════════════════════════════════════════════════════════════════════
# 2. GEO DEG loading + ortholog mapping
# ══════════════════════════════════════════════════════════════════════

def load_geo_deg(deg_csv_path, gene_col="Gene", logfc_col="logFC",
                 pval_col="P.Value"):
    """
    Load full DEG results from a GEO dataset.

    Returns DataFrame with columns: gene, logFC, p_value
    """
    df = pd.read_csv(deg_csv_path)
    # Standardize column names
    result = pd.DataFrame({
        "gene": df[gene_col].astype(str).tolist(),
        "logFC": df[logfc_col].astype(float).tolist(),
        "p_value": df[pval_col].astype(float).tolist() if pval_col in df.columns
                   else np.ones(len(df)),
    })
    # Drop rows with missing values
    result = result.dropna(subset=["gene", "logFC"])
    return result


def map_to_human_mygene(species_genes, species, cache_dir=CACHE_DIR):
    """
    Map mouse/rat/pig gene symbols to human using mygene + homologene.

    Approach:
      1. Query species genes → get homologene group IDs
      2. Extract human Entrez gene IDs from homologene groups (taxid=9606)
      3. Batch query human Entrez IDs → get human gene symbols

    Parameters
    ----------
    species_genes : list[str]
    species : str — 'mouse', 'rat', 'pig'

    Returns
    -------
    dict[str, str or None]
        {species_gene: human_gene_symbol or None}
    """
    import mygene

    species_to_taxid = {
        "mouse": 10090,
        "rat": 10116,
        "pig": 9823,
    }

    taxid = species_to_taxid.get(species)
    if taxid is None:
        raise ValueError(f"Unknown species: {species}")

    cache_file = cache_dir / f"ortholog_{species}_to_human.json"

    # Load cache
    cached = {}
    if cache_file.exists():
        with open(cache_file) as f:
            cached = json.load(f)

    # Only query genes not in cache
    to_query = [g for g in species_genes if g not in cached]
    if not to_query:
        return {g: cached.get(g) for g in species_genes}

    print(f"  Mapping {len(to_query)} {species} genes to human via mygene+homologene...")

    mg = mygene.MyGeneInfo()
    batch_size = 1000
    new_mappings = {}

    # ── Phase 1: Query species genes for homologene data ────────────
    gene_to_human_entrez = {}  # species_gene → human_entrez_str

    for start in tqdm(range(0, len(to_query), batch_size),
                      desc=f"  {species}→homologene", unit="batch"):
        batch = to_query[start:start + batch_size]
        try:
            results = mg.querymany(batch, scopes="symbol", species=taxid,
                                   fields="homologene")
        except Exception as e:
            print(f"    mygene query error: {e}")
            for g in batch:
                new_mappings[g] = None
            continue

        for item in results:
            query = item.get("query")
            hg = item.get("homologene", {})
            if hg and "genes" in hg:
                for hg_taxid, hg_entrez in hg["genes"]:
                    if hg_taxid == 9606:  # human
                        gene_to_human_entrez[query] = str(hg_entrez)
                        break
            if query not in gene_to_human_entrez:
                new_mappings[query] = None  # no human ortholog

    # ── Phase 2: Query human Entrez IDs for symbols ─────────────────
    unique_human_entrez = list(set(gene_to_human_entrez.values()))
    entrez_to_symbol = {}

    for start in range(0, len(unique_human_entrez), batch_size):
        batch = unique_human_entrez[start:start + batch_size]
        try:
            results = mg.querymany(batch, scopes="entrezgene", species="human",
                                   fields="symbol")
            for item in results:
                q = item.get("query")
                symbol = item.get("symbol")
                if symbol:
                    entrez_to_symbol[q] = symbol
        except Exception as e:
            print(f"    Human entrez lookup error: {e}")

    # ── Phase 3: Build final mapping ────────────────────────────────
    for species_gene, human_entrez in gene_to_human_entrez.items():
        human_symbol = entrez_to_symbol.get(human_entrez)
        new_mappings[species_gene] = human_symbol

    # Fill in any remaining genes
    for g in to_query:
        if g not in new_mappings:
            new_mappings[g] = None

    # Update cache
    cached.update(new_mappings)
    with open(cache_file, 'w') as f:
        json.dump(cached, f, indent=2, ensure_ascii=False)

    # Ensure all requested genes have a mapping
    for g in species_genes:
        if g not in cached:
            cached[g] = None

    n_mapped = sum(1 for v in cached.values() if v is not None)
    print(f"    Mapped {n_mapped}/{len(species_genes)} {species} genes to human")
    return {g: cached.get(g) for g in species_genes}


def build_human_ranked_list(deg_df, species, cache_dir=CACHE_DIR):
    """
    Map all genes in a DEG dataset to human, producing a ranked list
    sorted by logFC (descending) for GSEA input.

    Returns
    -------
    tuple[list[str], np.ndarray]
        (human_gene_symbols, logFC values) — sorted descending by logFC.
        Genes without human orthologs are dropped.
        Duplicate human genes take the max |logFC| entry.
    """
    species_genes = deg_df["gene"].tolist()
    mapping = map_to_human_mygene(species_genes, species, cache_dir)

    # Build human dataframe
    df = deg_df.copy()
    df["human_gene"] = df["gene"].map(mapping)
    df = df.dropna(subset=["human_gene"])
    df = df[df["human_gene"] != "None"]

    if df.empty:
        print(f"  WARNING: No human orthologs found for {species}")
        return [], np.array([])

    # Resolve duplicates: keep max |logFC|
    df["abs_logFC"] = df["logFC"].abs()
    df = df.sort_values("abs_logFC", ascending=False)
    df = df.drop_duplicates(subset="human_gene", keep="first")
    df = df.drop(columns=["abs_logFC"])

    # Sort descending by logFC for GSEA
    df = df.sort_values("logFC", ascending=False)

    return df["human_gene"].tolist(), df["logFC"].values


# ══════════════════════════════════════════════════════════════════════
# 3. Preranked GSEA against KEGG/GO
# ══════════════════════════════════════════════════════════════════════

def compute_pathway_es(ranked_genes, ranked_logfc, pathway_genes):
    """
    Compute enrichment score for a single pathway against a ranked gene list.

    Parameters
    ----------
    ranked_genes : list[str]
        Human gene symbols, ranked descending by logFC.
    ranked_logfc : np.ndarray
        logFC values in same order.
    pathway_genes : set[str]
        Genes in the pathway (human symbols).

    Returns
    -------
    dict with: ES, direction ("UP"/"DOWN"), n_pathway, n_overlap, n_overlap_ranked
    """
    n_in_ranked = 0
    gene_to_rank = {}
    for i, g in enumerate(ranked_genes):
        gene_to_rank[g] = i
        if g in pathway_genes:
            n_in_ranked += 1

    if n_in_ranked == 0:
        return {"ES": 0.0, "direction": None, "n_pathway": len(pathway_genes),
                "n_overlap": 0, "n_overlap_ranked": 0}

    # Build hit mask
    N = len(ranked_genes)
    hit_mask = np.zeros(N, dtype=bool)
    for i, g in enumerate(ranked_genes):
        if g in pathway_genes:
            hit_mask[i] = True

    es = enrichment_score(ranked_logfc, hit_mask, weight=1.0)
    direction = "UP" if es > 0 else "DOWN"

    return {
        "ES": float(es),
        "direction": direction,
        "n_pathway": len(pathway_genes),
        "n_overlap": int(hit_mask.sum()),
        "n_overlap_ranked": n_in_ranked,
    }


def run_gsea_on_dataset(ranked_genes, ranked_logfc, pathway_sets,
                        dataset_name, direction_label):
    """
    Compute ES for all pathways against one GEO dataset.

    Parameters
    ----------
    ranked_genes, ranked_logfc : GSEA-ranked input
    pathway_sets : dict[str, list[str]]
        {pathway_name: [gene_symbols]}
    dataset_name : str
    direction_label : str — 'kegg' or 'go'

    Returns
    -------
    list[dict] — one per pathway, sorted by |ES| descending
    """
    results = []
    for pw_name, pw_genes in tqdm(pathway_sets.items(),
                                   desc=f"  {dataset_name}/{direction_label}",
                                   unit="pathway"):
        res = compute_pathway_es(ranked_genes, ranked_logfc, set(pw_genes))
        res["pathway"] = pw_name
        res["source"] = f"GEO_{dataset_name}"
        res["database"] = direction_label.upper()
        results.append(res)

    results.sort(key=lambda r: abs(r["ES"]), reverse=True)
    return results


# ══════════════════════════════════════════════════════════════════════
# 4. Literature research via gpt-researcher
# ══════════════════════════════════════════════════════════════════════

def run_gpt_researcher(output_dir=CACHE_DIR):
    """
    Run gpt-researcher for MIRI pathway literature survey.
    Returns cached report if available, otherwise runs the research.

    Returns
    -------
    str — path to the generated report, or None if failed.
    """
    # Check for cached report first
    cache_report = output_dir / "literature_pathways_report.md"
    if cache_report.exists():
        age_hours = (time.time() - cache_report.stat().st_mtime) / 3600
        if age_hours < 168:  # 1 week
            print(f"  Using cached literature report ({age_hours:.0f}h old)")
            return str(cache_report)

    ref_dir = Path("/home/chang/arc/references")
    for f in sorted(ref_dir.glob("*miri_pathways_literature*"), reverse=True):
        age_hours = (time.time() - f.stat().st_mtime) / 3600
        if age_hours < 168:
            print(f"  Using existing literature report: {f}")
            import shutil
            shutil.copy(f, cache_report)
            return str(f)

    skill_dir = Path("/home/chang/arc/.agents/skills/gpt-researcher").resolve()
    runner = skill_dir / "run_research.py"

    if not runner.exists():
        print("  WARNING: gpt-researcher not found, skipping literature research")
        return None

    query = (
        "Comprehensive survey of all KEGG pathways and GO biological processes "
        "reported as significantly up-regulated or down-regulated in myocardial "
        "ischemia-reperfusion injury (MIRI). For each pathway specify: "
        "(1) exact pathway name with KEGG pathway ID (e.g. hsa04210) or GO term ID "
        "(e.g. GO:0006915), "
        "(2) regulation direction (UP = activated in MIRI, DOWN = suppressed), "
        "(3) supporting evidence from peer-reviewed literature with PMIDs. "
        "List at least 30 UP-regulated and 30 DOWN-regulated pathways, "
        "covering both KEGG and GO categories. "
        "Organize by regulation direction and database source."
    )

    print("  Running gpt-researcher for MIRI pathway literature survey...")
    print("  This may take several minutes...")

    try:
        result = subprocess.run(
            [sys.executable, str(runner),
             "--report-type", "deep",
             "--output", "miri_pathways_literature",
             query],
            cwd=str(output_dir),
            capture_output=True,
            text=True,
            timeout=600,
        )
        if result.returncode != 0:
            print(f"  WARNING: gpt-researcher failed: {result.stderr[:500]}")
            return None

        # Find the generated report
        report_files = sorted(output_dir.glob("*miri_pathways_literature*"))
        if not report_files:
            # Also check project references directory
            ref_dir = Path("/home/chang/arc/references")
            report_files = sorted(ref_dir.glob("*miri_pathways_literature*"))
        if report_files:
            report_path = report_files[0]
            print(f"  Literature report: {report_path}")
            # Copy to cache for future reference
            import shutil
            shutil.copy(report_path, output_dir / "literature_pathways_report.md")
            return str(report_path)
        else:
            print("  WARNING: No report file found after gpt-researcher run")
            return None

    except subprocess.TimeoutExpired:
        print("  WARNING: gpt-researcher timed out (600s)")
        return None
    except Exception as e:
        print(f"  WARNING: gpt-researcher error: {e}")
        return None


def parse_literature_report(report_path):
    """
    Parse the gpt-researcher markdown report to extract structured pathway info.

    Expects markdown tables with columns:
    | Type | Pathway/GO Term (ID) | Full Name | Direction | Evidence |

    Returns
    -------
    list[dict] — [{name, source (KEGG/GO), direction (UP/DOWN), pathway_id}, ...]
    """
    import re

    if report_path is None or not Path(report_path).exists():
        return []

    with open(report_path) as f:
        text = f.read()

    pathways = []
    seen = set()

    # Pattern: match table rows like:
    # | KEGG | hsa04668 | TNF signaling pathway | UP | ... |
    # | GO:BP | GO:0006954 | Inflammatory response | UP | ... |
    table_row_re = re.compile(
        r'^\|\s*(KEGG|GO:BP|GO)\s*\|\s*([^|]+)\s*\|\s*([^|]+)\s*\|\s*(UP|DOWN)\s*\|',
        re.MULTILINE
    )

    for match in table_row_re.finditer(text):
        source_type = match.group(1).strip()
        pw_id = match.group(2).strip()
        full_name = match.group(3).strip()
        direction = match.group(4).strip()

        # Normalize source type
        if source_type in ("GO:BP", "GO"):
            source = "GO"
        else:
            source = "KEGG"

        # Use full name as the primary name
        name = full_name

        # Store ID for matching purposes
        key = (name.lower(), source, direction)
        if key in seen:
            continue
        seen.add(key)

        pathways.append({
            "name": name,
            "source": source,
            "direction": direction,
            "pathway_id": pw_id,
            "evidence": "literature",
        })

    print(f"  Parsed {len(pathways)} pathways from literature report "
          f"({sum(1 for p in pathways if p['direction']=='UP')} UP, "
          f"{sum(1 for p in pathways if p['direction']=='DOWN')} DOWN)")
    return pathways


# ══════════════════════════════════════════════════════════════════════
# 5. Name normalization
# ══════════════════════════════════════════════════════════════════════

def normalize_pathway_name(name, source):
    """
    Normalize a pathway name with source prefix.

    - KEGG: "KEGG-Apoptosis" or "KEGG-TNF signaling pathway"
    - GO: "GO-Apoptotic Process"
    - Merged: "KEGG&GO-Apoptosis"

    Also strips GO ID suffixes like "(GO:0006915)".
    """
    import re

    # Strip GO ID
    name = re.sub(r'\s*\(GO:\d+\)\s*', '', name).strip()
    # Strip trailing/leading quotes
    name = name.strip("'\" ")

    # Already has a prefix
    if name.startswith(("KEGG-", "GO-", "KEGG&GO-")):
        return name

    return f"{source}-{name}"


# ══════════════════════════════════════════════════════════════════════
# 6. Combine, filter, merge, select
# ══════════════════════════════════════════════════════════════════════

def build_candidate_pool(geo_results_kegg, geo_results_go,
                         literature_pathways,
                         kegg_sets, go_sets,
                         l1000_set):
    """
    Build a pool of all candidate pathways with L1000 intersection data.

    Returns
    -------
    list[dict] — each with full metadata for filtering
    """
    candidates = {}
    # key = (normalized_name, direction, frozenset(l1000_genes))

    # From GEO KEGG
    for res in geo_results_kegg:
        pw_name = res["pathway"]
        full_genes = set(kegg_sets.get(pw_name, []))
        if not full_genes:
            continue
        l1000_genes = sorted(full_genes & l1000_set)
        geo_support = {
            "dataset": res["source"],
            "ES": res["ES"],
            "n_pathway_total": res["n_pathway"],
            "n_overlap_with_ranked": res["n_overlap_ranked"],
        }
        key = (normalize_pathway_name(pw_name, "KEGG"), res["direction"])
        if key in candidates:
            candidates[key]["geo_support"].append(geo_support)
        else:
            candidates[key] = {
                "name": normalize_pathway_name(pw_name, "KEGG"),
                "source_type": "KEGG",
                "direction": res["direction"],
                "full_genes": full_genes,
                "l1000_genes": l1000_genes,
                "n_total": len(full_genes),
                "n_l1000": len(l1000_genes),
                "geo_support": [geo_support],
                "literature_support": False,
                "pmids": [],
            }

    # From GEO GO
    for res in geo_results_go:
        pw_name = res["pathway"]
        full_genes = set(go_sets.get(pw_name, []))
        if not full_genes:
            continue
        l1000_genes = sorted(full_genes & l1000_set)
        geo_support = {
            "dataset": res["source"],
            "ES": res["ES"],
            "n_pathway_total": res["n_pathway"],
            "n_overlap_with_ranked": res["n_overlap_ranked"],
        }
        key = (normalize_pathway_name(pw_name, "GO"), res["direction"])
        if key in candidates:
            candidates[key]["geo_support"].append(geo_support)
        else:
            candidates[key] = {
                "name": normalize_pathway_name(pw_name, "GO"),
                "source_type": "GO",
                "direction": res["direction"],
                "full_genes": full_genes,
                "l1000_genes": l1000_genes,
                "n_total": len(full_genes),
                "n_l1000": len(l1000_genes),
                "geo_support": [geo_support],
                "literature_support": False,
                "pmids": [],
            }

    # From literature — look up genes from KEGG/GO cache
    for lit_pw in literature_pathways:
        lit_name = lit_pw["name"]
        lit_source = lit_pw["source"]
        lit_dir = lit_pw["direction"]
        lit_id = lit_pw.get("pathway_id", "")

        # Try to match against cached gene sets
        matched_genes = None
        if lit_source == "KEGG":
            # Try exact match first
            matched_genes = kegg_sets.get(lit_name, None)
            # Fuzzy match if not found
            if matched_genes is None:
                for k, v in kegg_sets.items():
                    if lit_name.lower() in k.lower() or k.lower() in lit_name.lower():
                        matched_genes = v
                        break
        elif lit_source == "GO":
            # Try matching by GO ID first (e.g., "GO:0006954" in "Inflammatory Response (GO:0006954)")
            if lit_id:
                for k, v in go_sets.items():
                    if lit_id in k:
                        matched_genes = v
                        break
            # Fall back to name matching
            if matched_genes is None:
                matched_genes = go_sets.get(lit_name, None)
            if matched_genes is None:
                for k, v in go_sets.items():
                    if lit_name.lower() in k.lower() or k.lower() in lit_name.lower():
                        matched_genes = v
                        break

        if matched_genes is None:
            continue

        full_genes = set(matched_genes)
        l1000_genes = sorted(full_genes & l1000_set)

        key = (normalize_pathway_name(lit_name, lit_source), lit_dir)
        if key in candidates:
            candidates[key]["literature_support"] = True
        else:
            candidates[key] = {
                "name": normalize_pathway_name(lit_name, lit_source),
                "source_type": lit_source,
                "direction": lit_dir,
                "full_genes": full_genes,
                "l1000_genes": l1000_genes,
                "n_total": len(full_genes),
                "n_l1000": len(l1000_genes),
                "geo_support": [],
                "literature_support": True,
                "pmids": [],
            }

    return list(candidates.values())


def filter_candidates(candidates):
    """
    Apply dual filters: n_l1000 > 5 AND coverage >= threshold.

    Returns (passing, rejected)
    """
    passing = []
    rejected = []
    for c in candidates:
        n_l1000 = c["n_l1000"]
        n_total = c["n_total"]
        passes_count = n_l1000 > 5
        passes_cov = passes_coverage_filter(n_l1000, n_total, TOTAL_HUMAN_GENES,
                                            N_LANDMARK, COVERAGE_FACTOR)
        if passes_count and passes_cov:
            passing.append(c)
        else:
            reasons = []
            if not passes_count:
                reasons.append(f"n_l1000={n_l1000}≤5")
            if not passes_cov:
                cov = n_l1000 / n_total if n_total > 0 else 0
                reasons.append(f"coverage={cov:.4f}<{COVERAGE_THRESHOLD:.4f}")
            c["reject_reason"] = "; ".join(reasons)
            rejected.append(c)
    return passing, rejected


def merge_identical_pathways(candidates):
    """
    Merge pathways with identical L1000 gene sets AND same direction.

    If a KEGG and a GO pathway merge, the new prefix is 'KEGG&GO-'.
    """
    merged = []
    used = set()
    n = len(candidates)

    for i in range(n):
        if i in used:
            continue
        pw_i = candidates[i]
        genes_i = set(pw_i["l1000_genes"])
        merged_pw = dict(pw_i)  # shallow copy
        merged_pw["merged_from"] = []

        for j in range(i + 1, n):
            if j in used:
                continue
            pw_j = candidates[j]
            genes_j = set(pw_j["l1000_genes"])

            if genes_i == genes_j and pw_i["direction"] == pw_j["direction"]:
                # Merge
                merged_pw["merged_from"].append(pw_j["name"])
                # Combine source types
                sources = {pw_i["source_type"], pw_j["source_type"]}
                if len(sources) > 1:
                    merged_pw["source_type"] = "KEGG&GO"
                    # Rebuild name
                    core = pw_i["name"]
                    for prefix in ["KEGG-", "GO-"]:
                        core = core.replace(prefix, "")
                    merged_pw["name"] = f"KEGG&GO-{core}"
                # Combine geo support
                merged_pw["geo_support"].extend(pw_j["geo_support"])
                # Combine literature
                merged_pw["literature_support"] = (
                    merged_pw["literature_support"] or pw_j["literature_support"]
                )
                used.add(j)

        merged.append(merged_pw)

    return merged


def score_and_rank(candidates):
    """
    Score each pathway for ranking within direction group.

    Score = GEO significance + literature support + ES magnitude + gene count bonus
    """
    for c in candidates:
        score = 0.0

        # GEO significance (primary driver)
        for gs in c.get("geo_support", []):
            es = gs.get("ES", 0.0)
            dataset = gs.get("dataset", "")
            if "GSE160516" in dataset:
                score += 10.0 * min(abs(es) / 0.3, 2.0)  # cap at 2x
            elif "GSE108940" in dataset:
                score += 5.0 * min(abs(es) / 0.3, 1.0)
            elif "GSE137293" in dataset:
                score += 3.0 * min(abs(es) / 0.3, 1.0)

        # Literature support
        if c.get("literature_support"):
            score += 8.0

        # Gene count bonus (more L1000 genes = more robust scoring)
        score += min(c["n_l1000"] / 20.0, 1.5)

        c["_score"] = score

    return candidates


def select_top_by_direction(candidates, target_up=30, target_down=30):
    """
    Select top pathways per direction using two-phase prioritization.

    Phase 1: Select ALL literature-supported pathways first, sorted by score,
             up to the target count. Literature-validated pathways are always
             prioritized over GEO-only pathways regardless of enrichment signal.
    Phase 2: Fill remaining slots with top-scored non-literature pathways.

    If ``--skip-literature`` is used, all candidates have literature_support=False,
    so Phase 1 selects nothing and this degrades to pure score-based selection.
    """
    up = [c for c in candidates if c["direction"] == "UP"]
    down = [c for c in candidates if c["direction"] == "DOWN"]

    def _select_group(group, target):
        lit = [c for c in group if c.get("literature_support")]
        non_lit = [c for c in group if not c.get("literature_support")]

        lit.sort(key=lambda c: c.get("_score", 0), reverse=True)
        non_lit.sort(key=lambda c: c.get("_score", 0), reverse=True)

        # Phase 1: literature pathways first (up to target)
        selected = lit[:min(target, len(lit))]
        n_lit = len(selected)
        remaining_slots = target - n_lit

        # Phase 2: fill remaining with top non-literature
        if remaining_slots > 0:
            n_fill = min(remaining_slots, len(non_lit))
            selected.extend(non_lit[:n_fill])

        if n_lit > 0:
            print(f"    Literature: {n_lit}, GEO-fill: {remaining_slots}/{len(non_lit)} available")
        return selected

    up_selected = _select_group(up, target_up)
    down_selected = _select_group(down, target_down)

    # Remove internal _score from output
    for c in up_selected + down_selected:
        c.pop("_score", None)

    return up_selected + down_selected, up, down


# ══════════════════════════════════════════════════════════════════════
# 7. Output
# ══════════════════════════════════════════════════════════════════════

def build_output_json(pathways, n_total, n_valid, kegg_sets, go_sets):
    """
    Build the final MIRI_pathways.json structure.

    Parameters
    ----------
    pathways : list[dict]
        Selected pathways with full metadata.
    n_total : int
        Total candidate pathways before filtering.
    n_valid : int
        Pathways after filtering + merging (before top-N selection).
    """
    output = {
        "name": "MIRI_pathways",
        "source": "Enrichr (KEGG_2021_Human + GO_Biological_Process_2023) + GEO GSEA + literature",
        "description": (
            f"Programmatically rebuilt MIRI pathway gene sets. "
            f"Sources: GEO DEG enrichment (GSE160516 mouse primary, GSE108940 rat, "
            f"GSE137293 pig) + gpt-researcher literature survey. "
            f"KEGG/GO gene sets from Enrichr, intersected with L1000 978 landmark genes. "
            f"Filters: n_l1000 > 5 AND n_l1000/n_total >= {N_LANDMARK}/{TOTAL_HUMAN_GENES} "
            f"(={COVERAGE_THRESHOLD:.4f}). "
            f"UP = pathway activated in MIRI (drug should DOWN-regulate). "
            f"DOWN = pathway suppressed in MIRI (drug should UP-regulate). "
        ),
        "generated_at": pd.Timestamp.now().isoformat(),
        "parameters": {
            "total_human_genome_genes": TOTAL_HUMAN_GENES,
            "n_landmark": N_LANDMARK,
            "coverage_factor": COVERAGE_FACTOR,
            "coverage_threshold": COVERAGE_THRESHOLD,
            "min_l1000_genes": MIN_L1000_GENES,
        },
        "n_pathways_total": n_total,
        "n_pathways_before_selection": n_valid,
        "n_pathways_final": len(pathways),
        "n_up": sum(1 for p in pathways if p["direction"] == "UP"),
        "n_down": sum(1 for p in pathways if p["direction"] == "DOWN"),
        "kegg_gene_sets_available": len(kegg_sets),
        "go_gene_sets_available": len(go_sets),
        "pathways": [],
    }

    for pw in pathways:
        entry = {
            "name": pw["name"],
            "direction": pw["direction"],
            "description": _build_description(pw),
            "genes": pw["l1000_genes"],
            "n_genes": pw["n_l1000"],
            "n_total_genes_before_l1000": pw["n_total"],
            "sources": _build_sources(pw),
            "merged_from": pw.get("merged_from", []) or None,
        }
        output["pathways"].append(entry)

    return output


def _build_description(pw):
    """Build a human-readable description for a pathway entry."""
    parts = []

    source = pw.get("source_type", "Unknown")
    name = pw["name"]
    direction = pw["direction"]

    if direction == "UP":
        dir_desc = "UP-regulated (activated) in MIRI. Drug should DOWN-regulate these genes."
    else:
        dir_desc = "DOWN-regulated (suppressed) in MIRI. Drug should UP-regulate these genes."

    parts.append(f"{source} pathway: {name}. {dir_desc}")

    # GEO evidence
    for gs in pw.get("geo_support", []):
        parts.append(f"{gs['dataset']}: ES={gs['ES']:.3f}")

    if pw.get("literature_support"):
        parts.append("Supported by MIRI literature.")

    if pw.get("merged_from"):
        parts.append(f"Merged from: {', '.join(pw['merged_from'])}")

    return " ".join(parts)


def _build_sources(pw):
    """Build the sources list for a pathway entry."""
    sources = []
    for gs in pw.get("geo_support", []):
        sources.append(gs["dataset"])
    if pw.get("literature_support"):
        sources.append("literature")
    return sorted(set(sources)) if sources else ["unknown"]


def save_intermediate(data, filename, output_dir=CACHE_DIR):
    """Save intermediate data to .claude/ directory."""
    path = output_dir / filename
    with open(path, 'w') as f:
        json.dump(data, f, indent=2, ensure_ascii=False, default=str)
    print(f"  → Saved intermediate: {path}")
    return path


# ══════════════════════════════════════════════════════════════════════
# 8. Main pipeline
# ══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Rebuild MIRI pathway gene sets from GEO + literature")
    parser.add_argument("--target-up", type=int, default=30,
                        help="Target number of UP pathways (default: 30)")
    parser.add_argument("--target-down", type=int, default=30,
                        help="Target number of DOWN pathways (default: 30)")
    parser.add_argument("--refresh-cache", action="store_true",
                        help="Force re-fetch from Enrichr")
    parser.add_argument("--skip-literature", action="store_true",
                        help="Skip gpt-researcher literature survey")
    parser.add_argument("--total-human-genes", type=int, default=20000,
                        help="Total human protein-coding genes (default: 20000)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output JSON path (default: reference-genesets/MIRI_pathways.json)")
    args = parser.parse_args()

    global TOTAL_HUMAN_GENES, COVERAGE_THRESHOLD
    TOTAL_HUMAN_GENES = args.total_human_genes
    COVERAGE_THRESHOLD = (N_LANDMARK / TOTAL_HUMAN_GENES) * COVERAGE_FACTOR

    print("=" * 65)
    print("Rebuilding MIRI Pathway Gene Sets")
    print(f"  Target: ~{args.target_up} UP, ~{args.target_down} DOWN")
    print(f"  Filters: n_genes > 5, coverage >= {COVERAGE_THRESHOLD:.4f}")
    print(f"  Total human genes: {TOTAL_HUMAN_GENES}")
    print("=" * 65)

    # ── 1. Load L1000 ────────────────────────────────────────────────
    print("\n[1/7] Loading L1000 landmark genes...")
    l1000_set = load_l1000_genes(L1000_CSV)
    print(f"  Loaded {len(l1000_set)} landmark genes")

    # ── 2. Fetch KEGG/GO gene sets ───────────────────────────────────
    print("\n[2/7] Fetching KEGG and GO gene sets from Enrichr...")

    if args.refresh_cache:
        for f in CACHE_DIR.glob("KEGG_2021_Human.json"):
            f.unlink()
        for f in CACHE_DIR.glob("GO_Biological_Process_2023.json"):
            f.unlink()

    kegg_sets = fetch_enrichr_library("KEGG_2021_Human")
    go_sets = fetch_enrichr_library("GO_Biological_Process_2023")

    print(f"  KEGG pathways: {len(kegg_sets)}")
    print(f"  GO BP terms:   {len(go_sets)}")

    # ── 3. GEO enrichment ────────────────────────────────────────────
    print("\n[3/7] Running GEO-based pathway enrichment...")

    # Primary: GSE160516 (mouse)
    geo_configs = [
        {
            "name": "GSE160516",
            "csv": GEO_DIR / "GSE160516_deg_results.csv",
            "species": "mouse",
            "gene_col": "Gene",
            "logfc_col": "logFC",
            "pval_col": "P.Value",
            "weight": 1.0,  # primary
        },
        {
            "name": "GSE108940",
            "csv": GEO_DIR / "GSE108940_WT_deg_results.csv",
            "species": "rat",
            "gene_col": "Gene",
            "logfc_col": "logFC",
            "pval_col": "P.Value",
            "weight": 0.5,
        },
        {
            "name": "GSE137293",
            "csv": GEO_DIR / "GSE137293_deg_results.csv",
            "species": "pig",
            "gene_col": "gene",
            "logfc_col": "log2FC",
            "pval_col": "P_Value",
            "weight": 0.3,
        },
    ]

    all_geo_kegg = []
    all_geo_go = []

    for cfg in geo_configs:
        deg_path = cfg["csv"]
        if not deg_path.exists():
            print(f"  SKIP {cfg['name']}: {deg_path} not found")
            continue

        print(f"\n  ── {cfg['name']} ({cfg['species']}) ──")
        deg_df = load_geo_deg(deg_path, cfg["gene_col"],
                              cfg["logfc_col"], cfg["pval_col"])
        print(f"    Loaded {len(deg_df)} DEGs")

        ranked_genes, ranked_logfc = build_human_ranked_list(
            deg_df, cfg["species"])
        print(f"    Mapped to {len(ranked_genes)} human genes")

        if len(ranked_genes) == 0:
            print(f"    WARNING: No human-mapped genes, skipping")
            continue

        # Run GSEA against KEGG
        kegg_results = run_gsea_on_dataset(
            ranked_genes, ranked_logfc, kegg_sets, cfg["name"], "kegg")
        # Run GSEA against GO
        go_results = run_gsea_on_dataset(
            ranked_genes, ranked_logfc, go_sets, cfg["name"], "go")

        all_geo_kegg.extend(kegg_results)
        all_geo_go.extend(go_results)

        print(f"    KEGG top-5 UP: "
              + ", ".join(f"{r['pathway'][:30]}({r['ES']:.3f})"
                          for r in kegg_results[:5] if r["direction"] == "UP"))
        print(f"    KEGG top-5 DOWN: "
              + ", ".join(f"{r['pathway'][:30]}({r['ES']:.3f})"
                          for r in kegg_results[:5] if r["direction"] == "DOWN"))

    # Save intermediate
    save_intermediate(
        {"kegg": all_geo_kegg[:100], "go": all_geo_go[:100]},
        "geo_enrichment_top100.json")

    # ── 4. Literature research ───────────────────────────────────────
    print("\n[4/7] Literature research...")

    literature_pathways = []
    if args.skip_literature:
        print("  Skipped (--skip-literature)")
    else:
        try:
            report_path = run_gpt_researcher()
            if report_path:
                literature_pathways = parse_literature_report(report_path)
                save_intermediate(literature_pathways,
                                 "literature_pathways_parsed.json")
            else:
                print("  No literature data available, proceeding with GEO only")
        except Exception as e:
            print(f"  Literature research failed ({e}), proceeding with GEO only")

    # ── 5. Build candidate pool ──────────────────────────────────────
    print("\n[5/7] Building candidate pool...")
    candidates = build_candidate_pool(
        all_geo_kegg, all_geo_go,
        literature_pathways,
        kegg_sets, go_sets,
        l1000_set,
    )
    n_total = len(candidates)
    print(f"  Total candidates (before filtering): {n_total}")

    # ── 6. Filter + Merge ────────────────────────────────────────────
    print("\n[6/7] Filtering and merging...")
    passing, rejected = filter_candidates(candidates)
    print(f"  Passed filters: {len(passing)}")
    print(f"  Rejected: {len(rejected)}")
    if rejected:
        reasons = {}
        for r in rejected:
            reason = r.get("reject_reason", "unknown")
            reasons[reason] = reasons.get(reason, 0) + 1
        for reason, count in sorted(reasons.items(), key=lambda x: -x[1]):
            print(f"    {reason}: {count}")

    save_intermediate([{"name": r["name"], "direction": r["direction"],
                        "reason": r.get("reject_reason")} for r in rejected],
                      "rejected_pathways.json")

    # Merge identical
    merged = merge_identical_pathways(passing)
    n_merges = len(passing) - len(merged)
    if n_merges > 0:
        print(f"  Merged {n_merges} pathways into existing ones")
    print(f"  After merging: {len(merged)} pathways")

    # Score and select
    scored = score_and_rank(merged)
    selected, all_up, all_down = select_top_by_direction(
        scored, args.target_up, args.target_down)
    print(f"  Selected: {sum(1 for p in selected if p['direction']=='UP')} UP "
          f"+ {sum(1 for p in selected if p['direction']=='DOWN')} DOWN")

    if len(all_up) > args.target_up:
        print(f"  ({(len(all_up))} UP available, selected top {args.target_up})")
    if len(all_down) > args.target_down:
        print(f"  ({(len(all_down))} DOWN available, selected top {args.target_down})")

    # ── 7. Output ────────────────────────────────────────────────────
    print("\n[7/7] Writing output...")

    output = build_output_json(selected, n_total, len(passing), kegg_sets, go_sets)

    output_path = args.output or str(REF_DIR / "MIRI_pathways.json")
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, ensure_ascii=False)
    print(f"  → Saved final output: {output_path}")

    # Summary
    print("\n" + "=" * 65)
    print("Rebuild complete!")
    print(f"  Final pathways: {output['n_pathways_final']} "
          f"({output['n_up']} UP + {output['n_down']} DOWN)")
    print(f"  Sources: {output['n_pathways_before_selection']} passed filters "
          f"from {output['n_pathways_total']} total candidates")
    print(f"  Intermediate files: {CACHE_DIR}")
    print("=" * 65)


if __name__ == "__main__":
    main()
