"""Hermetic fixtures for the raw-count and corpus release tests.

Builds a tiny but contract-shaped source tree (two genes, one tenv plate, two
Tahoe shards), the two accepted-style raw payloads through the public builders,
and the frozen axis/crosswalk/drug-metadata inputs. Every value is hand-worked in
the tests; the fixture exists only to drive the real build/replay boundary.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def _write_seal(source, seal):
    lines = []
    for path in sorted(source.rglob('*.parquet')):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {path.relative_to(source)}")
    seal.write_text("\n".join(lines) + "\n")


def write_axis(root):
    axis = Path(root) / "genes.csv"
    pd.DataFrame({
        "symbol": ["A", "B"],
        "hgnc_id": ["HGNC:1", "HGNC:2"],
        "ensembl_id": ["E1", "E2"],
    }).to_csv(axis, index=False)
    crosswalk = Path(root) / "gene_crosswalk.csv"
    pd.DataFrame([
        {"dataset_id": "tahoe", "source_gene_id": "E1", "hgnc_id": "HGNC:1", "match_type": "direct"},
        {"dataset_id": "tahoe", "source_gene_id": "3", "hgnc_id": "HGNC:1", "match_type": "direct"},
        {"dataset_id": "tahoe", "source_gene_id": "E2", "hgnc_id": "HGNC:2", "match_type": "direct"},
        {"dataset_id": "tahoe", "source_gene_id": "5", "hgnc_id": "HGNC:2", "match_type": "direct"},
    ]).to_csv(crosswalk, index=False)
    return axis, crosswalk


def write_sources(root):
    """Write the tiny immutable source trees under ``root/data``."""
    data = Path(root) / "data"

    # --- tenv: one plate, 2 genes, treated + DMSO + Media + Blank wells ---
    tenv = data / "GSE262419/source"
    tenv.mkdir(parents=True)
    pd.DataFrame({
        "Genes": ["A_1", "B_3"],
        "Plate1-A1": [1, 2],   # treated
        "Plate1-A2": [3, 4],   # DMSO vehicle
        "Plate1-A3": [5, 6],   # Media vehicle
        "Plate1-A4": [0, 1],   # technical blank
    }).to_csv(tenv / "GSE262419_Plate1.csv.gz", index=False)
    pd.DataFrame({
        "Plate_ID": [1, 1, 1, 1],
        "Well_ID": ["A01", "A02", "A03", "A04"],
        "Chemical_name": ["ethanol", "VEH", "Media", "Blank"],
        "Chemical_Concentration_uM": [1, 0, 0, 0],
        "Chemical_Index": ["0", "1", "2", "3"],
    }).to_csv(tenv / "GSE262419_hash.csv.gz", index=False)
    annotation = data / "Annotation/source/hgnc_complete_set.txt"
    annotation.parent.mkdir(parents=True)
    pd.DataFrame({
        "hgnc_id": ["HGNC:1", "HGNC:2"],
        "symbol": ["A", "B"],
        "alias_symbol": ["", ""],
        "prev_symbol": ["", ""],
        "entrez_id": ["1", "3"],
        "ensembl_gene_id": ["E1", "E2"],
    }).to_csv(annotation, sep="\t", index=False)
    compounds = data / "Annotation/local/cm_compound_annotation.csv"
    compounds.parent.mkdir(parents=True)
    pd.DataFrame([{
        "lookup_name": "ethanol", "dataset_inchikey": "",
        "parent_inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "smiles": "CCO",
    }]).to_csv(compounds, index=False)
    tenv_manifest = data / "GSE262419/manifest.sha256"
    tenv_manifest.write_text(
        f"{hashlib.sha256((tenv / 'GSE262419_Plate1.csv.gz').read_bytes()).hexdigest()}"
        "  source/GSE262419_Plate1.csv.gz\n"
        f"{hashlib.sha256((tenv / 'GSE262419_hash.csv.gz').read_bytes()).hexdigest()}"
        "  source/GSE262419_hash.csv.gz\n")

    # --- Tahoe: two shards, plate1 has a control, plate2 does not ---
    tahoe = data / "Tahoe-100M/source"
    (tahoe / "metadata").mkdir(parents=True)
    (tahoe / "data").mkdir()
    pd.DataFrame([
        {"sample": "A", "plate": "plate1", "drugname_drugconc": "[('drug', 1, 'uM')]"},
        {"sample": "C", "plate": "plate1", "drugname_drugconc": "[('DMSO_TF', 0, 'uM')]"},
        {"sample": "B", "plate": "plate2", "drugname_drugconc": "[('drug', 1, 'uM')]"},
    ]).to_parquet(tahoe / "metadata/sample_metadata.parquet", index=False)
    pd.DataFrame([{"Cell_ID_Cellosaur": "L1", "cell_name": "line1"}]).to_parquet(
        tahoe / "metadata/cell_line_metadata.parquet", index=False)
    pd.DataFrame([
        {"gene_symbol": "G1", "ensembl_id": "E1", "token_id": 3},
        {"gene_symbol": "G2", "ensembl_id": "E2", "token_id": 5},
    ]).to_parquet(tahoe / "metadata/gene_metadata.parquet", index=False)
    shards = [
        [{"plate": "plate1", "sample": "A", "cell_line_id": "L1",
          "genes": [0, 3, 5], "expressions": [-2, 2, 7]},
         {"plate": "plate2", "sample": "B", "cell_line_id": "L1",
          "genes": [0, 3], "expressions": [-2, 8]}],
        [{"plate": "plate1", "sample": "C", "cell_line_id": "L1",
          "genes": [0, 3], "expressions": [-2, 1]}],
    ]
    for index, shard in enumerate(shards):
        pq.write_table(pa.Table.from_pylist(shard),
                       tahoe / f"data/train-{index:05d}-of-00002.parquet")
    _write_seal(tahoe, data / "Tahoe-100M/manifest.sha256")

    drug_metadata = data / "drug_metadata.parquet"
    pd.DataFrame([{"drug": "drug", "canonical_smiles": "CCO"}]).to_parquet(
        drug_metadata, index=False)
    return {"data_root": data, "drug_metadata": drug_metadata}


def write_payloads(root, axis, crosswalk):
    """Build the two accepted-style raw payloads with the public builders."""
    from bases.native import build_tenv
    from bases.tahoe import build_tahoe

    sources = write_sources(root)
    data = sources["data_root"]
    tenv_payload = Path(root) / "tenv_payload"
    build_tenv(data / "GSE262419/source",
               data / "Annotation/source/hgnc_complete_set.txt",
               data / "Annotation/local/cm_compound_annotation.csv", tenv_payload)
    mapping = pd.DataFrame([
        {"gene_index": 0, "symbol": "A", "hgnc_id": "HGNC:1", "ensembl_id": "E1",
         "native_column": 0, "native_symbol": "G1", "native_ensembl_id": "E1",
         "native_token_id": "3", "resolution": "exact_ensembl"},
        {"gene_index": 1, "symbol": "B", "hgnc_id": "HGNC:2", "ensembl_id": "E2",
         "native_column": 1, "native_symbol": "G2", "native_ensembl_id": "E2",
         "native_token_id": "5", "resolution": "exact_ensembl"},
    ])
    tahoe_payload = Path(root) / "tahoe_payload"
    build_tahoe(data / "Tahoe-100M/source", mapping,
                data / "Tahoe-100M/manifest.sha256", tahoe_payload,
                expected_shards=2, expected_genes=2, batch_size=1)
    return {**sources, "tenv_payload": tenv_payload,
            "tahoe_payload": tahoe_payload}


def write_raw_release(root, axis, crosswalk):
    from bases import raw_counts

    fixture = write_payloads(root, axis, crosswalk)
    release = Path(root) / "raw_counts"
    raw_counts.build(release, mode="source-accepted-replay",
                     data_root=fixture["data_root"], axis_path=axis,
                     crosswalk_path=crosswalk,
                     tenv_payload=fixture["tenv_payload"],
                     tahoe_payload=fixture["tahoe_payload"], fixture=True)
    raw_counts.verify(release)
    return {**fixture, "raw_counts": release}


def write_corpus(root, axis, crosswalk):
    from bases import corpus

    fixture = write_raw_release(root, axis, crosswalk)
    release = Path(root) / "corpus"
    corpus.build(release, raw_counts_dir=fixture["raw_counts"], axis_path=axis,
                 drug_metadata=fixture["drug_metadata"], fixture=True)
    corpus.verify(release, fixture=True)
    return {**fixture, "corpus": release}
