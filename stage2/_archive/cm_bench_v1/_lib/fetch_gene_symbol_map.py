#!/usr/bin/env python3
"""
Regenerate the bundled Ensembl-gene-ID -> HGNC-symbol map (`input/`).

`gse244740` (ScreenSeq) identifies its genes in Ensembl gene-ID space; the other
two CM datasets (`dtoxs_primary`, `dtoxs_ipsc`) use HGNC symbols. To build the
common whole-transcriptome gene universe we need an Ensembl -> symbol map, which
is bundled as a frozen input asset so the benchmark itself never touches the
network and the universe is byte-stable.

The map is derived from the Ensembl GRCh38 release-104 `gene` GTF (gse244740's
series matrix documents Ensembl annotation v104):

    https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz

Only the gene_id + gene_name attributes are read. This is an asset-generation
step, not part of the benchmark run.

Usage:
    python _lib/fetch_gene_symbol_map.py [--gtf /path/Homo_sapiens.GRCh38.104.gtf.gz]
"""

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import CM_DATA_ROOT, INPUT_DIR  # noqa: E402

GTF_URL = ("https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/"
           "Homo_sapiens.GRCh38.104.gtf.gz")
DEFAULT_OUT = INPUT_DIR / "ensembl_to_symbol.tsv"


def gse244740_ensembl_ids():
    df = pd.read_csv(CM_DATA_ROOT / "gse244740"
                     / "GSE244740_processed_data_counts.txt.gz",
                     sep="\t", usecols=["Ensembl"])
    return {str(e).split(".")[0] for e in df["Ensembl"]}


def parse_gtf(gtf_path, keep):
    """gene_id -> gene_name from the gene feature lines of a GTF."""
    out = {}
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = {}
            for part in fields[8].split(";"):
                part = part.strip()
                if not part:
                    continue
                k, _, v = part.partition(" ")
                attrs[k] = v.strip('"')
            gid = attrs.get("gene_id")
            if gid in keep:
                out[gid] = attrs.get("gene_name", "")
                if len(out) == len(keep):
                    break
    return out


def download_gtf(dest):
    import urllib.request

    print(f"Downloading {GTF_URL}")
    urllib.request.urlretrieve(GTF_URL, dest)
    return dest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", default=None,
                    help="cached GRCh38.104 GTF; downloads if omitted")
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()

    ids = gse244740_ensembl_ids()
    gtf = Path(args.gtf) if args.gtf else Path("/tmp/Homo_sapiens.GRCh38.104.gtf.gz")
    if not gtf.exists():
        download_gtf(gtf)
    raw = parse_gtf(gtf, ids)
    rows = [{"ensembl_id": e, "gene_symbol": raw.get(e, "")}
            for e in sorted(ids)]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    n = sum(1 for r in rows if r["gene_symbol"])
    print(f"Wrote {out} ({n}/{len(rows)} of gse244740 genes mapped)")


if __name__ == "__main__":
    main()
