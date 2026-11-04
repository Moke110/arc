#!/usr/bin/env python3
"""Build the four CM native-axis corpus h5ad files (ADR 0033, ticket #15).

One row per raw well; ``X`` = raw counts summed per resolved HGNC id on the
source's native gene axis (dense ``int32``); one ``log1p_cp10k`` layer
(``float32``); the source's own row order; vehicle controls and technical blanks
kept, pharmacology tiers not applied.  The 18,164-gene CM axis is a *projection
target* owned by ``tables/genes.csv`` (ticket #16), not a corpus universe.

    python stage2/02_dataset_corpus/execution/build_corpus.py --corpus all
    python stage2/02_dataset_corpus/execution/build_corpus.py --corpus tenv

Each invocation writes one corpus at a time, atomically (``.tmp`` then rename)
with ``track_times=False``, and records counts / content digest / file hash /
bytes / wall time / peak RSS in ``execution/build_record/corpus_build.json``.
``execution/build_record/excluded_compounds.csv`` records the non-small-molecule
perturbagens (8 keys, 143 wells) and the 32 retained technical blanks.

Other modes:
    --check         print the four raw-source counts and exit
    --repeat-check  build the selected corpus twice, assert identical file
                    sha256 (gate G5), record ``corpus_reproducibility.json``
    --verify        G1-G6 gates: NOT IMPLEMENTED until ticket #17
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import platform
import resource
import shutil
import sys
import tempfile
import time
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from cmdata import paths, qa, sources  # noqa: E402
from cmdata.annotation import AnnotationError, load_assets  # noqa: E402
from cmdata.corpus_h5ad import (CorpusError, build_corpus, content_digest,  # noqa: E402
                                corpus_stats, write_corpus)
from cmdata.drugs import (CompoundError, CompoundTable,  # noqa: E402
                          build_compound_table, tier_counts)

MIN_FREE_GIB = 12.0


class BuildError(RuntimeError):
    pass


def _iso_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _available_gib() -> float:
    with open("/proc/meminfo") as handle:
        for line in handle:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1024 ** 2
    raise BuildError("cannot read MemAvailable from /proc/meminfo")


def _environment() -> dict:
    return {
        "python": platform.python_version(),
        "anndata": version("anndata"),
        "h5py": version("h5py"),
        "numpy": version("numpy"),
        "pandas": version("pandas"),
        "platform": platform.platform(),
    }


def _load_sources(index) -> dict:
    """Load the four raw sources exactly as the v1 writer did.

    The compound table is global (identity must not depend on which corpus is
    being built), so every build loads all four sources even when it writes one
    corpus; the other three are used only for their row metadata.
    """
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


def _excluded_compounds(table: CompoundTable) -> pd.DataFrame:
    """The 8 non-small-molecule keys (143 wells) + the 32 blank wells."""
    exclusions = table.exclusions
    rows = []
    compounds = exclusions[exclusions["entity_type"] != "technical_blank"]
    for (name, entity), group in compounds.groupby(
            ["compound_name", "entity_type"], sort=True):
        rows.append({
            "dataset_id": ";".join(sorted(group["dataset_id"].astype(str))),
            "compound_name": name,
            "compound_code": ";".join(sorted({
                str(code) for code in group["compound_code"] if str(code)})),
            "entity_type": entity,
            "n_wells": int(group["n_wells"].sum()),
            "reason": f"non-small-molecule perturbagen ({entity})",
        })
    blanks = exclusions[exclusions["entity_type"] == "technical_blank"]
    if len(blanks):
        rows.append({
            "dataset_id": ";".join(sorted(blanks["dataset_id"].astype(str))),
            "compound_name": "Blank",
            "compound_code": "",
            "entity_type": "technical_blank",
            "n_wells": int(blanks["n_wells"].sum()),
            "reason": ("technical blank (concentration 0; kept as a corpus row, "
                       "neither control nor perturbation)"),
        })
    return pd.DataFrame(rows, columns=[
        "dataset_id", "compound_name", "compound_code", "entity_type",
        "n_wells", "reason"])


def _write_frame_atomic(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    frame.to_csv(temp, index=False, lineterminator="\n")
    os.replace(temp, path)


def _write_json_atomic(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temp, path)


def _load_record() -> dict:
    if paths.CORPUS_BUILD_RECORD.exists():
        return json.loads(paths.CORPUS_BUILD_RECORD.read_text())
    return {"corpora": {}}


def _gate_payload(results: list[qa.GateResult]) -> list[dict]:
    return [{"name": result.name, "passed": bool(result.passed),
             "skipped": bool(result.skipped), "detail": result.detail}
            for result in results]


def _build_one(corpus_id: str, loaded: dict, table: CompoundTable,
               index) -> dict:
    source_id = paths.CORPUS_TO_SOURCE[corpus_id]
    source = loaded[source_id]
    print(f"Building {corpus_id} from {source_id} "
          f"({source.meta.shape[0]:,} source wells)")
    started = time.time()
    expected = (source.meta.shape[0],
                int(source.meta["n_raw_rows"].sum()))
    adata = build_corpus(corpus_id, source, table, index)
    stats = corpus_stats(adata, corpus_id)
    gates = qa.corpus_structure_gates(
        adata, corpus_id, paths.EXPECTED_CORPUS_COUNTS[corpus_id])
    gates += qa.source_value_check({source_id: source})
    failed = [gate for gate in gates if not gate.passed and not gate.skipped]
    if failed:
        raise BuildError(
            f"{corpus_id}: QA gates failed: "
            + "; ".join(f"{gate.name}: {gate.detail}" for gate in failed))
    digest = content_digest(adata)
    write_info = write_corpus(adata, paths.corpus_h5ad(corpus_id))
    wall = time.time() - started
    peak_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    entry = {
        "dataset_id": corpus_id,
        "source_dataset": source_id,
        "accession": paths.source_accession(source_id),
        "rows": stats["rows"],
        "perturbed": stats["perturbed"],
        "control": stats["control"],
        "blank": stats["blank"],
        "chemicals": stats["chemicals"],
        "axis_length": stats["axis"],
        "distinct_wells": stats["distinct_wells"],
        "source": {
            "wells": expected[0],
            "raw_gene_ids": len(source.gene_ids),
            "resolved_genes": len(source.gene_map.resolved_ids),
            "raw_design_rows": expected[1],
        },
        "n_raw_rows": {
            "sum": int(adata.obs["n_raw_rows"].sum()),
            "collapsed_wells": int((adata.obs["n_raw_rows"] > 1).sum()),
            "max": int(adata.obs["n_raw_rows"].max()),
        },
        "container": "dense int32 X + dense float32 layer",
        "dtypes": {"X": str(adata.X.dtype),
                   "log1p_cp10k": str(adata.layers["log1p_cp10k"].dtype)},
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
        "content_sha256": digest,
        "file": write_info["file"],
        "file_sha256": write_info["file_sha256"],
        "bytes": write_info["bytes"],
        "built_at": _iso_now(),
        "wall_seconds": round(wall, 1),
        "peak_rss_kib": int(peak_kib),
        "peak_rss_gib": round(peak_kib / 1024 ** 2, 2),
        "gates": _gate_payload(gates),
    }
    print(f"  {stats['rows']:,} rows x {stats['axis']:,} genes; "
          f"{stats['perturbed']:,} perturbed + {stats['control']:,} control + "
          f"{stats['blank']:,} blank; {stats['chemicals']:,} chemicals")
    print(f"  content_sha256 {digest}")
    print(f"  file_sha256    {write_info['file_sha256']} "
          f"({write_info['bytes']:,} B); wall {wall:.1f}s, "
          f"peak RSS {peak_kib / 1024 ** 2:.2f} GiB")
    print(f"  {sum(1 for gate in gates if gate.passed)}/{len(gates)} gates passed")
    del adata
    gc.collect()
    loaded[source_id].counts = None
    return entry


def build(args) -> None:
    started = time.time()
    available = _available_gib()
    print(f"available memory: {available:.1f} GiB (need >= {MIN_FREE_GIB:.0f} GiB)")
    if available < MIN_FREE_GIB:
        raise BuildError(
            f"only {available:.1f} GiB available; refusing to start a build "
            f"that needs at least {MIN_FREE_GIB:.0f} GiB")

    print("Loading pinned annotation assets")
    index, annotation = load_assets()
    print(f"  HGNC genes: {len(index.genes):,}; "
          f"compound annotation rows: {len(annotation.table):,}")

    print("Loading sources")
    loaded = _load_sources(index)

    print("Building compound table")
    table = build_compound_table(
        {dataset_id: loaded[dataset_id].meta
         for dataset_id in paths.SOURCE_DATASETS}, annotation)
    print(f"  {table.compounds['drug_key'].nunique():,} drug keys, "
          f"{len(table.compounds):,} key/name rows; tiers {tier_counts(table)}")
    print(f"  exclusions (not in the corpus): {len(table.exclusions)}")

    targets = (paths.CORPUS_IDS if args.corpus == "all" else [args.corpus])
    record = _load_record()
    record["ticket"] = "#15"
    record["adr"] = "0033"
    record["step"] = "stage2/02_dataset_corpus"
    record["updated_at"] = _iso_now()
    record["build"] = {
        "source_commit": paths.SOURCE_COMMIT,
        "writer": "track_times=False, gzip, .tmp -> atomic rename",
        "content_digest_recipe": (
            "sha256(X.tobytes() | layers['log1p_cp10k'].tobytes() | "
            "obs.to_csv(index=False, lineterminator='\\n') | "
            "var.to_csv(lineterminator='\\n'))  # v1 recipe minus uns"),
        "row_rule": ("vehicle controls + technical blanks + every "
                     "SMILES-resolvable small molecule; no tier filter"),
        "well_id_policy": (
            "obs['well'] is the loader string verbatim, so it keeps the "
            "historical id prefix (screen_ctox::, temp_env::); dataset_id and "
            "control_group use the v2 ids (sctox::, tenv::). Ticket #17 "
            "remaps those two for the v1 comparison, never the well id."),
        "pharmacology_tiers": (
            "pharm_atc_moa/pharm_mesh are copied per row; all three tier "
            "booleans live in tables/compounds.csv for a consumer to filter "
            "on; no tier is applied when building a corpus"),
        "environment": _environment(),
    }

    for corpus_id in targets:
        record["corpora"][corpus_id] = _build_one(
            corpus_id, loaded, table, index)

    exclusions = _excluded_compounds(table)
    _write_frame_atomic(exclusions, paths.EXCLUDED_COMPOUNDS_RECORD)
    record["excluded_compounds"] = {
        "file": str(paths.EXCLUDED_COMPOUNDS_RECORD.relative_to(paths.REPO_ROOT)),
        "keys": int((exclusions["entity_type"] != "technical_blank").sum()),
        "wells": int(exclusions.loc[
            exclusions["entity_type"] != "technical_blank", "n_wells"].sum()),
        "blank_wells": int(exclusions.loc[
            exclusions["entity_type"] == "technical_blank", "n_wells"].sum()),
    }
    record["corpora"] = {key: record["corpora"][key]
                         for key in paths.CORPUS_IDS if key in record["corpora"]}
    _write_json_atomic(record, paths.CORPUS_BUILD_RECORD)

    print(f"\nWrote {len(targets)} corpus file(s) in {time.time() - started:.1f}s:")
    for corpus_id in targets:
        entry = record["corpora"][corpus_id]
        print(f"  {entry['file']}  {entry['rows']:,} x {entry['axis_length']:,}  "
              f"{entry['chemicals']:,} chemicals  {entry['bytes']:,} B")
    print(f"Build record: "
          f"{paths.CORPUS_BUILD_RECORD.relative_to(paths.REPO_ROOT)}")


def repeat_check(corpus_id: str) -> int:
    """Gate G5: build one corpus twice, require identical file sha256."""
    index, annotation = load_assets()
    loaded = _load_sources(index)
    table = build_compound_table(
        {dataset_id: loaded[dataset_id].meta
         for dataset_id in paths.SOURCE_DATASETS}, annotation)
    hashes, digests, sizes = [], [], []
    scratch = Path(tempfile.mkdtemp(prefix="g5_", dir=HERE))
    try:
        for run in (1, 2):
            print(f"G5 run {run}/2 for {corpus_id}")
            adata = build_corpus(
                corpus_id, loaded[paths.CORPUS_TO_SOURCE[corpus_id]], table,
                index)
            digests.append(content_digest(adata))
            info = write_corpus(adata, scratch / f"run{run}_{corpus_id}.h5ad")
            hashes.append(info["file_sha256"])
            sizes.append(info["bytes"])
            del adata
            gc.collect()
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    identical = hashes[0] == hashes[1] and digests[0] == digests[1]
    payload = {
        "ticket": "#15",
        "gate": "G5",
        "corpus": corpus_id,
        "checked_at": _iso_now(),
        "runs": 2,
        "content_sha256": digests[0],
        "file_sha256_run1": hashes[0],
        "file_sha256_run2": hashes[1],
        "bytes": sizes[0],
        "identical_content": digests[0] == digests[1],
        "identical_file_sha256": hashes[0] == hashes[1],
        "result": "PASS" if identical else "FAIL",
    }
    _write_json_atomic(payload, paths.REPRODUCIBILITY_RECORD)
    print(f"G5 {payload['result']}: file sha256 run1={hashes[0]}")
    print(f"                               file sha256 run2={hashes[1]}")
    print(f"record: {paths.REPRODUCIBILITY_RECORD.relative_to(paths.REPO_ROOT)}")
    return 0 if identical else 2


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--corpus", default="all",
                        choices=["all", *paths.CORPUS_IDS],
                        help="corpus to build (default: all four, one at a time)")
    parser.add_argument("--check", action="store_true",
                        help="print the four raw-source counts and exit")
    parser.add_argument("--repeat-check", action="store_true",
                        help="G5: build the selected corpus twice and compare "
                             "its file sha256 (does not write output/)")
    parser.add_argument("--verify", action="store_true",
                        help="G1-G6 gates; reserved for ticket #17")
    args = parser.parse_args(argv)

    if args.check:
        return sources.check_sources()
    if args.verify:
        print("TODO(#17): the --verify entry point (G1-G6, including the v1 "
              "reconstruction gate) is implemented by ticket #17, not #15.",
              file=sys.stderr)
        return 3
    if args.repeat_check:
        if args.corpus == "all":
            print("--repeat-check needs one corpus (G5 checks one small "
                  "corpus twice); pass --corpus <id>.", file=sys.stderr)
            return 2
        if _available_gib() < MIN_FREE_GIB:
            print(f"only {_available_gib():.1f} GiB available; aborting.",
                  file=sys.stderr)
            return 2
        try:
            return repeat_check(args.corpus)
        except (CorpusError, CompoundError, AnnotationError,
                sources.SourceError, BuildError) as exc:
            print(f"\nREPEAT CHECK FAILED: {exc}", file=sys.stderr)
            return 2
    try:
        build(args)
    except (CorpusError, CompoundError, AnnotationError,
            sources.SourceError, BuildError) as exc:
        print(f"\nBUILD FAILED: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
