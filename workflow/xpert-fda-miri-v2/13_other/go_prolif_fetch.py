#!/usr/bin/env python3
"""
13_other — cache proliferation-related GO (biological process) data.

Pulls two small, reproducible reference pieces from the Gene Ontology:

  1. QuickGO term search  -> every BP term whose *name* contains
     "proliferation" (this is the "proliferation GO pathway" vocabulary).
  2. QuickGO annotation search rooted at GO:0008283 (cell population
     proliferation). The default relationship expansion returns human protein
     annotations for the whole proliferation subtree, i.e. the proliferation
     processes themselves plus their regulation / positive-regulation /
     negative-regulation descendants.

Raw responses are cached under output/cache/ so all downstream scoring is
offline and reproducible. Human taxon only (9606).

Outputs (all under output/cache/):
  go_prolif_terms.json    {go_id: name, ...}
  go_prolif_annotations.csv  one row per (annotation id, go_id, symbol, ...)

Usage:
    python 13_other/go_prolif_fetch.py            # download + cache
    python 13_other/go_prolif_fetch.py --offline  # reuse cache (must exist)
"""

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import requests

OUT = Path(__file__).resolve().parent / "output" / "cache"
QGO = "https://www.ebi.ac.uk/QuickGO/services"
ROOT_GO = "GO:0008283"  # cell population proliferation
TAXON = 9606
UA = {"Accept": "application/json"}


def _get(url, params, tries=5):
    for attempt in range(tries):
        r = requests.get(url, params=params, headers=UA, timeout=120)
        if r.status_code == 200:
            return r.json()
        time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"QuickGO {url} params={params} -> HTTP {r.status_code}: {r.text[:400]}")


def fetch_term_search():
    """All BP terms matching 'proliferation' (name/synonym/definition)."""
    terms, page, total_pages = {}, 1, 1
    while page <= total_pages:
        d = _get(f"{QGO}/ontology/go/search",
                 {"query": "proliferation", "aspect": "biological_process",
                  "page": page, "pageSize": 100})
        info = d.get("pageInfo", {})
        if info:
            total_pages = info.get("total", 1)
        rows = d.get("results", [])
        for t in rows:
            tid = t.get("id")
            if tid and tid.startswith("GO:"):
                terms.setdefault(tid, t.get("name"))
        if not rows:
            break
        page += 1
        time.sleep(0.2)
    return terms


def fetch_annotations(start_page=1):
    """Human protein annotations under the cell-population-proliferation subtree."""
    page, total_pages = start_page, start_page
    rows = []
    while page <= total_pages:
        d = _get(f"{QGO}/annotation/search",
                 {"goId": ROOT_GO, "taxonId": TAXON,
                  "includeFields": "goName", "page": page, "pageSize": 100})
        info = d.get("pageInfo", {})
        total_pages = info.get("total", page) if info else page
        results = d.get("results", [])
        rows.extend(results)
        if page % 25 == 0 or not results:
            print(f"  annotations page {page}/{total_pages} "
                  f"({len(rows)} rows collected)", flush=True)
        page += 1
        time.sleep(0.2)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true")
    ap.add_argument("--resume", action="store_true", help="append rows if a partial cache exists")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    term_file = OUT / "go_prolif_terms.json"
    ann_file = OUT / "go_prolif_annotations.csv"

    if args.offline:
        if not (term_file.exists() and ann_file.exists()):
            sys.exit("offline cache missing; run without --offline first")
        terms = json.loads(term_file.read_text())
        with open(ann_file, newline="") as f:
            n_rows = sum(1 for _ in f) - 1
        print(f"offline cache: {len(terms)} terms, {n_rows} annotation rows")
        return

    terms = fetch_term_search()
    name_matched = {g: n for g, n in terms.items()
                    if n and "proliferation" in n.lower()}
    term_file.write_text(json.dumps(name_matched, indent=1, sort_keys=True))
    print(f"term cache written: {len(name_matched)} terms", flush=True)

    start_page = 1
    if ann_file.exists() and args.resume:
        with open(ann_file, newline="") as f:
            n_existing = sum(1 for _ in f) - 1
        start_page = 1 + (n_existing // 25)
        print(f"resume: {n_existing} cached rows, restarting at page {start_page}", flush=True)

    print("annotation pull under GO:0008283 ...", flush=True)
    anns = fetch_annotations(start_page=start_page)
    keep = ["id", "geneProductId", "goId", "goName", "qualifier", "symbol",
            "goEvidence", "evidenceCode", "assignedBy"]
    mode = "a" if start_page > 1 else "w"
    with open(ann_file, mode, newline="") as f:
        w = csv.DictWriter(f, fieldnames=keep, extrasaction="ignore")
        if start_page == 1:
            w.writeheader()
        w.writerows(anns)
    print("annotation cache written", flush=True)

    goids = {a.get("goId") for a in anns}
    overlap = name_matched.keys() & goids
    print(f"term search: {len(terms)} BP hits, {len(name_matched)} name-contains-'proliferation'")
    print(f"annotation pull under {ROOT_GO}: {len(anns)} rows, {len(goids)} distinct go_ids")
    print(f"name-matched terms present in the annotation subtree: {len(overlap)}/{len(name_matched)}")
    print(f"cached -> {term_file}")
    print(f"cached -> {ann_file}")


if __name__ == "__main__":
    main()
