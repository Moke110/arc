#!/usr/bin/env python3
"""Build the Stage 2 CM corpus and one dataset bundle.

Single entry point (ADR 0026): sources -> frozen corpus -> filtered views.

    python stage2/02_cm_data/build_cm_datasets.py --bundle datasets-1

The build is deterministic for a fixed corpus: the gene axis order is frozen,
row order follows the source files, and a rerun produces the same content
hashes (recorded in ``manifest.csv``).  Any new gene or drug outside
``corpus/v1/`` is a hard failure; promoting the corpus is a deliberate
``--rebuild-corpus`` edit.

Other modes:
    --check             print per-source well/control/compound counts
    --corpus-only       freeze or verify the corpus, build no bundle
    --rebuild-corpus    deliberately replace corpus/v1 (drift reset)
    --skip-stage1       skip the read-only Stage 1 regression gate
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from cmdata import corpus, paths, qa, sources, views  # noqa: E402
from cmdata.annotation import AnnotationError, HgncIndex, load_assets  # noqa: E402
from cmdata.drugs import (CompoundError, CompoundTable,  # noqa: E402
                          annotate_source_meta, build_compound_table,
                          tier_counts)
from cmdata.genes import compute_universe  # noqa: E402


class BuildError(RuntimeError):
    pass


def _content_digest(adata) -> str:
    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(adata.X).tobytes())
    digest.update(np.ascontiguousarray(adata.layers["log1p_cp10k"]).tobytes())
    digest.update(adata.obs.to_csv(index=False, lineterminator="\n").encode())
    digest.update(adata.var.to_csv(lineterminator="\n").encode())
    means = adata.uns["control_mean"]
    digest.update(means.to_csv(lineterminator="\n").encode())
    return digest.hexdigest()


def _load_sources(index: HgncIndex) -> dict:
    entrez = sources.load_dtoxs_entrez_map()
    loaded = {}
    for dataset_id in paths.SOURCE_DATASETS:
        started = time.time()
        source = sources.load_source(dataset_id, index, entrez)
        loaded[dataset_id] = source
        print(f"  loaded {dataset_id}: {source.meta.shape[0]:,} wells x "
              f"{len(source.gene_ids):,} source genes "
              f"({time.time() - started:.1f}s)")
    return loaded


def _crosswalk(sources: dict, index: HgncIndex, universe) -> pd.DataFrame:
    rows = []
    for dataset_id in paths.SOURCE_DATASETS:
        gene_map = sources[dataset_id].gene_map
        for raw_id, hgnc_id, match, verified in zip(
                gene_map.raw_ids, gene_map.hgnc_ids, gene_map.match_types,
                gene_map.verified):
            symbol = index.ref(hgnc_id).symbol if hgnc_id is not None else ""
            rows.append({
                "dataset_id": dataset_id,
                "source_gene_id": raw_id,
                "symbol": symbol,
                "hgnc_id": hgnc_id or "",
                "match_type": match,
                "entrez_verified": bool(verified),
                "in_corpus": hgnc_id in universe.position if hgnc_id else False,
            })
    return pd.DataFrame(rows)


def _compounds_with_membership(table: CompoundTable,
                               built_views: dict) -> pd.DataFrame:
    frame = table.compounds.copy()
    for view_name, view in built_views.items():
        stem = view_name.removesuffix(".h5ad")
        used = set(view.adata.obs.loc[~view.adata.obs["is_control"], "drug_key"])
        frame[stem] = frame["drug_key"].isin(used)
    return frame


def _write_manifest(bundle: Path, built_views: dict, source_commit: str) -> pd.DataFrame:
    rows = []
    for view_name, view in built_views.items():
        path = bundle / view_name
        adata = view.adata
        rows.append({
            "view": view_name,
            "file": view_name,
            "n_obs": adata.n_obs,
            "n_perturbed": view.n_perturbed,
            "n_controls": view.n_controls,
            "n_drugs": view.n_drugs,
            "n_genes": adata.n_vars,
            "n_control_groups": adata.uns["control_mean"].shape[0],
            "corpus_version": paths.CORPUS_VERSION,
            "source_commit": source_commit,
            "content_sha256": _content_digest(adata),
            "bytes": path.stat().st_size,
        })
    return pd.DataFrame(rows)


def build(args) -> None:
    started = time.time()
    if args.bundle == paths.DEFAULT_BUNDLE:
        bundle = paths.BUNDLE_DIR
    else:
        bundle = paths.CM_DATA_DIR / args.bundle
    temp = bundle.with_name(bundle.name + ".tmp")

    print("Loading pinned annotation assets")
    index, annotation = load_assets()
    print(f"  HGNC genes: {len(index.genes):,}; "
          f"compound annotation rows: {len(annotation.table):,}")

    print("Loading sources")
    loaded = _load_sources(index)
    mappings = [loaded[dataset_id].gene_map for dataset_id in paths.SOURCE_DATASETS]
    computed_universe = compute_universe(index, mappings)
    print(f"  four-source intersection: {len(computed_universe):,} genes")

    print("Building compound table")
    source_metas = {dataset_id: loaded[dataset_id].meta
                    for dataset_id in paths.SOURCE_DATASETS}
    table = build_compound_table(source_metas, annotation)
    print(f"  {table.compounds['drug_key'].nunique():,} drug keys, "
          f"{len(table.compounds):,} key/name rows; "
          f"tiers {tier_counts(table)}")
    print(f"  exclusions (not in corpus): {len(table.exclusions)}")

    source_maps = {dataset_id: loaded[dataset_id].gene_map
                   for dataset_id in paths.SOURCE_DATASETS}
    source_stats = {}
    universe_ids = set(computed_universe.genes["hgnc_id"])
    for dataset_id in paths.SOURCE_DATASETS:
        gene_map = loaded[dataset_id].gene_map
        resolved = gene_map.resolved_ids
        source_stats[dataset_id] = {
            "counts": dict(gene_map.counts),
            "resolved_unique": len(resolved),
            "entrez_verified": int(gene_map.verified.sum()),
            "in_intersection": len(resolved & universe_ids),
        }

    if args.rebuild_corpus and paths.CORPUS_DIR.exists():
        shutil.rmtree(paths.CORPUS_DIR)
    if paths.GENES_CSV.exists():
        frozen_universe, frozen_compounds, frozen_source_sets = (
            corpus.load_frozen_corpus())
        corpus.check_drift(computed_universe, table, source_maps, index,
                           frozen_universe, frozen_compounds,
                           frozen_source_sets)
        universe = frozen_universe
        drift_note = ("Rerun matched `corpus/v1/` exactly: no gene or drug "
                      "drift in any source.")
        print("  corpus/v1 verified (no drift)")
    else:
        corpus.freeze_corpus(computed_universe, table, source_stats,
                             source_maps, index)
        universe = computed_universe
        drift_note = ("`corpus/v1/` frozen by this build; there was no previous "
                      "version to drift against.")
        print(f"  froze corpus/v1 ({len(universe):,} genes)")

    if args.corpus_only:
        print("Corpus only; no bundle written.")
        return

    print("Mapping source counts to the corpus gene axis")
    for dataset_id, source in loaded.items():
        source.mapped_counts = views.map_counts(source.counts, source.gene_map,
                                                universe)
        source.meta = annotate_source_meta(source.meta, table)
        unresolved = source.meta.loc[
            (~source.meta["is_control"]) & (~source.meta["is_blank"])
            & (source.meta["entity_type"] == "small_molecule")
            & (source.meta["drug_key"].astype(str) == "")]
        if len(unresolved):
            raise BuildError(
                f"{dataset_id}: small-molecule rows without a drug_key: "
                f"{unresolved[['compound_name', 'well']].head(5).to_dict('records')}")

    if temp.exists():
        shutil.rmtree(temp)
    temp.mkdir(parents=True)
    try:
        built_views = {}
        for view_name, datasets, tier in paths.VIEWS:
            print(f"Building {view_name} ({tier}; {', '.join(datasets)})")
            view = views.build_view(view_name, datasets, tier, loaded,
                                    universe, table)
            built_views[view_name] = view
            print(f"  {view.n_perturbed:,} perturbed + {view.n_controls:,} "
                  f"controls; {view.n_drugs:,} drugs")

        expected = paths.EXPECTED_VIEWS
        for view_name, counts in expected.items():
            view = built_views[view_name]
            actual = {"perturbed": view.n_perturbed,
                      "controls": view.n_controls, "drugs": view.n_drugs}
            if actual != counts:
                raise BuildError(f"{view_name}: expected {counts}, got {actual}")

        tier_results = []
        for group, expected_counts in paths.EXPECTED_TIER_COUNTS.items():
            actual_counts = tier_counts(table, paths.TIER_DATASETS[group])
            ok = actual_counts == expected_counts
            tier_results.append(qa.GateResult(
                f"tier counts {group}", ok,
                f"{actual_counts}" if ok
                else f"expected {expected_counts}, got {actual_counts}"))

        gate_results = (qa.cross_view_gates(built_views, universe)
                        + tier_results
                        + qa.source_value_check(loaded))
        if not paths.STAGE1_OUTPUT.exists():
            stage1_results = [qa.GateResult(
                "stage1 regression", True,
                f"SKIPPED: baseline missing ({paths.STAGE1_OUTPUT})",
                skipped=True)]
        else:
            stage1_results = qa.stage1_regression(loaded, universe)

        failed = [result for result in gate_results + stage1_results
                  if not result.passed and not result.skipped]
        if failed:
            raise BuildError("QA gates failed: " + "; ".join(
                f"{result.name}: {result.detail}" for result in failed))

        print("Writing bundle")
        for view_name, view in built_views.items():
            view.adata.write_h5ad(temp / view_name, compression="gzip")

        source_commit = paths.SOURCE_COMMIT
        compounds = _compounds_with_membership(table, built_views)
        compounds.to_csv(temp / "compounds.csv", index=False, lineterminator="\n")
        crosswalk = _crosswalk(loaded, index, universe)
        crosswalk.to_csv(temp / "gene_crosswalk.csv", index=False,
                         lineterminator="\n")
        manifest = _write_manifest(temp, built_views, source_commit)
        manifest.to_csv(temp / "manifest.csv", index=False, lineterminator="\n")
        qa.write_report(
            temp / "qa_report.md", views=built_views, sources=loaded,
            universe=universe, table=table, gate_results=gate_results,
            stage1_results=stage1_results, manifest=manifest,
            drift_note=drift_note)

        if bundle.exists():
            shutil.rmtree(bundle)
        temp.rename(bundle)
    except Exception:
        if temp.exists():
            shutil.rmtree(temp)
        raise

    print(f"\nBundle written: {bundle}")
    print(manifest.to_string(index=False))
    print(f"Total build time: {time.time() - started:.1f}s")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bundle", default=paths.DEFAULT_BUNDLE,
                        help="bundle directory name under stage2/02_cm_data")
    parser.add_argument("--check", action="store_true",
                        help="print per-source counts and exit")
    parser.add_argument("--corpus-only", action="store_true",
                        help="freeze/verify the corpus and stop")
    parser.add_argument("--rebuild-corpus", action="store_true",
                        help="deliberately replace corpus/v1 (drift reset)")
    args = parser.parse_args(argv)

    if args.check:
        return sources.check_sources()
    try:
        build(args)
    except (corpus.DriftError, views.ViewError, sources.SourceError,
            CompoundError, AnnotationError, BuildError) as exc:
        print(f"\nBUILD FAILED: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
