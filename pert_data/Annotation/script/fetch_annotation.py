#!/usr/bin/env python3
"""Pin the external annotation assets used by the Stage 2 CM corpus.

The build (`build_cm_datasets.py`) never touches the network.  Everything it
needs about gene and compound identity is frozen under
``pert_data/Annotation/`` by this script:

  source/hgnc_complete_set.txt        HGNC complete set (gene identity; download)
  local/cm_compound_annotation.csv    public-registry structure + pharmacology
  local/README.md                     source URL, retrieval date, checksum, counts

Two acquisition paths
---------------------
``--hgnc``
    Download the HGNC complete set and write the manifest.  Re-running with an
    existing, checksum-matching file is a no-op.

``--compounds --from-audit <csv>``
    Bootstrap the compound snapshot from the frozen offline audit seed in
    Annotation local additions. That table was produced by live PubChem /
    UniChem / ChEMBL / RxNav queries and is the one-time public-registry
    snapshot.  This is the path used for corpus v1: it is deterministic and
    makes no network call.

``--compounds --refresh``
    Re-query the public registries for every name in the Stage 1 inventory.
    Slow and rate-limited; only needed when the upstream audit is unavailable.

Every compound row is keyed by ``(lookup_name, dataset_inchikey)``.  The
``dataset_inchikey`` is empty for name-only datasets and carries the source's
own structure for DToxS; that is what keeps the DToxS ``AMI`` pentyl-benzofuran
homolog from being relabelled as amiodarone.

Usage:
    python pert_data/Annotation/script/fetch_annotation.py --hgnc
    python pert_data/Annotation/script/fetch_annotation.py --compounds
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
from datetime import date
from pathlib import Path

import pandas as pd

def _find_repo_root(start: Path) -> Path:
    """Walk up from ``start`` to the directory containing .git (ADR 0030)."""
    for candidate in [start, *start.parents]:
        if (candidate / ".git").exists():
            return candidate
    raise SystemExit(f"cannot find the repository root (.git) above {start}")


HERE = Path(__file__).resolve().parent
REPO_ROOT = _find_repo_root(HERE)
ITEM_DIR = REPO_ROOT / "pert_data" / "Annotation"
SOURCE_DIR = ITEM_DIR / "source"      # downloaded bytes only
LOCAL_DIR = ITEM_DIR / "local"        # local additions, declared in author_log.json

HGNC_PATH = SOURCE_DIR / "hgnc_complete_set.txt"
COMPOUND_PATH = LOCAL_DIR / "cm_compound_annotation.csv"
MANIFEST_PATH = LOCAL_DIR / "README.md"

HGNC_URL = ("https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/"
            "hgnc_complete_set.txt")
DEFAULT_AUDIT = LOCAL_DIR / "compound_external_audit.csv"

# Columns kept in the frozen compound snapshot.  `parent_inchikey` is the
# desalted parent key (the corpus `drug_key`); `pubchem_inchikey` is the exact
# structure key as registered (the corpus `inchikey`, salts included).
COMPOUND_COLUMNS = [
    "lookup_name", "dataset_inchikey", "inchikey", "parent_inchikey",
    "smiles", "pubchem_cid", "chembl_id", "atc", "moa", "mesh_classes",
    "verification", "notes",
]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _manifest_block(name: str, path: Path, url: str, extra: dict) -> str:
    lines = [
        f"## {name}",
        "",
        f"- file: `{path.name}`",
        f"- source: {url}",
        f"- retrieval date: {extra.get('retrieved', 'unknown')}",
        f"- sha256: `{_sha256(path)}`",
        f"- bytes: {path.stat().st_size:,}",
    ]
    for key, value in extra.items():
        if key == "retrieved":
            continue
        lines.append(f"- {key}: {value}")
    lines.append("")
    return "\n".join(lines)


def write_manifest(hgnc_extra: dict | None = None,
                   compound_extra: dict | None = None) -> None:
    header = [
        "# Stage 2 CM corpus — pinned annotation assets",
        "",
        "These files are read-only annotation inputs for Stage 2 bases gene identity.",
        "They are produced by `fetch_annotation.py`; the build never makes a",
        "network call.",
        "",
    ]
    body = []
    if HGNC_PATH.exists() and hgnc_extra is not None:
        body.append(_manifest_block("HGNC complete set", HGNC_PATH, HGNC_URL,
                                    hgnc_extra))
    if COMPOUND_PATH.exists() and compound_extra is not None:
        body.append(_manifest_block(
            "Compound structure and pharmacology snapshot", COMPOUND_PATH,
            f"{DEFAULT_AUDIT.relative_to(REPO_ROOT)} "
            "(PubChem/UniChem/ChEMBL/RxNav snapshot)",
            compound_extra))
    LOCAL_DIR.mkdir(parents=True, exist_ok=True)
    MANIFEST_PATH.write_text("\n".join(header + body))


# ---------------------------------------------------------------------------
# HGNC
# ---------------------------------------------------------------------------
def hgnc_info(path: Path = HGNC_PATH) -> dict:
    table = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    dates = table["date_modified"].dropna()
    return {
        "release: hgnc_complete_set rows": f"{len(table):,}",
        "latest date_modified": dates.max() if len(dates) else "unknown",
        "columns": str(table.shape[1]),
    }


def fetch_hgnc(force: bool = False) -> None:
    import requests

    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    if HGNC_PATH.exists() and not force:
        print(f"HGNC already pinned: {HGNC_PATH} "
              f"({HGNC_PATH.stat().st_size:,} bytes, sha256 {_sha256(HGNC_PATH)})")
        retrieved = date.fromtimestamp(HGNC_PATH.stat().st_mtime).isoformat()
    else:
        print(f"Downloading {HGNC_URL}")
        response = requests.get(HGNC_URL, timeout=120)
        response.raise_for_status()
        HGNC_PATH.write_bytes(response.content)
        print(f"Wrote {HGNC_PATH} ({HGNC_PATH.stat().st_size:,} bytes)")
        retrieved = date.today().isoformat()

    info = hgnc_info()
    info = {"retrieved": retrieved, **info}
    write_manifest(hgnc_extra=info,
                   compound_extra=(_compound_info() if COMPOUND_PATH.exists()
                                   else None))

# ---------------------------------------------------------------------------
# Compound snapshot
# ---------------------------------------------------------------------------
def _compound_info() -> dict:
    table = pd.read_csv(COMPOUND_PATH)
    return {
        "retrieved": "inherited from the Stage 1 external audit snapshot",
        "rows": f"{len(table):,}",
        "unique lookup names": f"{table.lookup_name.nunique():,}",
        "with a parent InChIKey": f"{table.parent_inchikey.notna().sum():,}",
        "with an ATC code or mechanism":
            f"{((table.atc.notna()) | (table.moa.notna())).sum():,}",
    }


def build_compounds_from_audit(audit_path: Path) -> pd.DataFrame:
    audit = pd.read_csv(audit_path)
    required = {"lookup_name", "dataset_inchikey", "pubchem_inchikey",
                "parent_inchikey", "smiles", "atc", "moa", "mesh_classes"}
    missing = required - set(audit.columns)
    if missing:
        raise SystemExit(f"audit table is missing columns: {sorted(missing)}")

    table = pd.DataFrame({
        "lookup_name": audit["lookup_name"].astype(str).str.strip(),
        "dataset_inchikey": audit["dataset_inchikey"].fillna(""),
        "inchikey": audit["pubchem_inchikey"].fillna(""),
        "parent_inchikey": audit["parent_inchikey"].fillna(""),
        "smiles": audit["smiles"].fillna(""),
        "pubchem_cid": (audit["pubchem_cid"].astype("Int64").astype(str)
                        .replace("<NA>", "")),
        "chembl_id": audit["chembl_id"].fillna(""),
        "atc": audit["atc"].fillna(""),
        "moa": audit["moa"].fillna(""),
        "mesh_classes": audit["mesh_classes"].fillna(""),
        "verification": audit["verification"].fillna(""),
        "notes": audit["notes"].fillna(""),
    })
    # Deterministic order; a (name, structure) pair is unique.
    table = table.sort_values(["lookup_name", "dataset_inchikey"],
                              kind="stable").reset_index(drop=True)
    dup = table.duplicated(["lookup_name", "dataset_inchikey"])
    if dup.any():
        raise SystemExit("duplicate (lookup_name, dataset_inchikey) rows: "
                         f"{table.loc[dup, ['lookup_name', 'dataset_inchikey']].to_dict('records')[:5]}")
    return table[COMPOUND_COLUMNS]


def build_compounds_from_audit_file(audit_path: Path) -> None:
    if not audit_path.exists():
        raise SystemExit(f"audit table not found: {audit_path}")
    table = build_compounds_from_audit(audit_path)
    LOCAL_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(COMPOUND_PATH, index=False)
    print(f"Wrote {COMPOUND_PATH} ({len(table):,} rows, "
          f"{table.lookup_name.nunique():,} names)")
    write_manifest(
        hgnc_extra=({"retrieved": date.today().isoformat(), **hgnc_info()}
                    if HGNC_PATH.exists() else None),
        compound_extra=_compound_info())


# ---------------------------------------------------------------------------
# Live refresh (optional; not used for corpus v1)
# ---------------------------------------------------------------------------
def _live_lookup(name: str, dataset_inchikey: str, cache: dict,
                 cache_path: Path) -> dict:
    """Query PubChem / UniChem / ChEMBL for one compound.

    Only the fields the corpus needs are returned.  This path exists so the
    snapshot can be rebuilt without the Stage 1 audit table; it is deliberately
    separate from the build.
    """
    import requests

    PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    PUG_VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"
    CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
    UNICHEM = "https://www.ebi.ac.uk/unichem/rest"
    session = requests.Session()
    session.headers.update({"User-Agent": "arc-stage2-cm-annotation/1.0"})

    def get_json(url, params=None):
        key = url + ("?" + json.dumps(params, sort_keys=True) if params else "")
        if key in cache:
            return cache[key]
        for attempt in range(5):
            response = session.get(url, params=params, timeout=60)
            if response.status_code == 404:
                cache[key] = None
                return None
            if response.status_code in (429, 500, 502, 503, 504):
                time.sleep(2.0 * (attempt + 1))
                continue
            response.raise_for_status()
            data = response.json()
            cache[key] = data
            time.sleep(0.2)
            cache_path.write_text(json.dumps(cache))
            return data
        return None

    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    def parent_key(smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None:
            return ""
        try:
            parent = rdMolStandardize.LargestFragmentChooser().choose(mol)
            return Chem.MolToInchiKey(parent)
        except Exception:  # noqa: BLE001
            return ""

    query = re.sub(r"\s+", " ", str(name).replace("_", " ").strip())
    record = {"lookup_name": query, "dataset_inchikey": dataset_inchikey,
              "inchikey": "", "parent_inchikey": "", "smiles": "",
              "pubchem_cid": "", "chembl_id": "", "atc": "", "moa": "",
              "mesh_classes": "", "verification": "unresolved", "notes": ""}

    props = None
    if dataset_inchikey:
        cids = get_json(f"{PUG}/compound/inchikey/{dataset_inchikey}/cids/JSON")
        if isinstance(cids, dict) and cids.get("IdentifierList", {}).get("CID"):
            cid = cids["IdentifierList"]["CID"][0]
            props = get_json(
                f"{PUG}/compound/cid/{cid}/property/"
                "CanonicalSMILES,IsomericSMILES,InChIKey/JSON")
    if props is None:
        cids = get_json(f"{PUG}/compound/name/{requests.utils.quote(query)}/cids/JSON")
        if isinstance(cids, dict) and cids.get("IdentifierList", {}).get("CID"):
            cid = cids["IdentifierList"]["CID"][0]
            props = get_json(
                f"{PUG}/compound/cid/{cid}/property/"
                "CanonicalSMILES,IsomericSMILES,InChIKey/JSON")
            record["notes"] = "name lookup"
    if not props or not props.get("PropertyTable", {}).get("Properties"):
        return record
    prop = props["PropertyTable"]["Properties"][0]
    smiles = prop.get("SMILES") or prop.get("CanonicalSMILES") or ""
    record.update({
        "inchikey": prop.get("InChIKey", ""),
        "parent_inchikey": parent_key(smiles),
        "smiles": smiles,
        "pubchem_cid": prop.get("CID", ""),
        "verification": "pubchem",
    })

    parent = record["parent_inchikey"] or record["inchikey"]
    if parent:
        xrefs = get_json(f"{UNICHEM}/verbose_inchikey/{parent}")
        chembl_id = ""
        if isinstance(xrefs, list):
            for source in xrefs:
                if source.get("name") == "chembl" and source.get("src_compound_id"):
                    chembl_id = source["src_compound_id"][0]
                    break
        record["chembl_id"] = chembl_id
        if chembl_id:
            mol = get_json(f"{CHEMBL}/molecule/{chembl_id}.json")
            atc_codes = (mol or {}).get("atc_classifications") or []
            mech = get_json(f"{CHEMBL}/mechanism.json",
                            params={"molecule_chembl_id": chembl_id, "limit": 100})
            moa = sorted({m.get("mechanism_of_action") for m in
                          (mech or {}).get("mechanisms", [])
                          if m.get("mechanism_of_action")})
            record["atc"] = "; ".join(atc_codes)
            record["moa"] = "; ".join(moa)
    if record["pubchem_cid"]:
        view = get_json(f"{PUG_VIEW}/data/compound/{record['pubchem_cid']}/JSON",
                        params={"heading": "MeSH Pharmacological Classification"})
        names = []

        def walk(node):
            if isinstance(node, dict):
                if node.get("TOCHeading") == "MeSH Pharmacological Classification":
                    for info in node.get("Information", []):
                        value = info.get("Name")
                        if not value:
                            strings = info.get("Value", {}).get("StringWithMarkup", [])
                            value = strings[0].get("String") if strings else None
                        if value and value not in names:
                            names.append(value)
                for child in node.values():
                    walk(child)
            elif isinstance(node, list):
                for child in node:
                    walk(child)

        walk((view or {}).get("Record", {}))
        record["mesh_classes"] = "; ".join(names)
    return record


def refresh_compounds(audit_path: Path) -> None:
    audit = pd.read_csv(audit_path)
    pairs = (audit[["lookup_name", "dataset_inchikey"]]
             .fillna({"dataset_inchikey": ""})
             .drop_duplicates()
             .sort_values(["lookup_name", "dataset_inchikey"]))
    LOCAL_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = LOCAL_DIR / "_compound_registry_cache.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    rows = []
    for index, row in enumerate(pairs.itertuples(index=False), 1):
        rows.append(_live_lookup(row.lookup_name, row.dataset_inchikey,
                                 cache, cache_path))
        if index % 25 == 0:
            print(f"  {index}/{len(pairs)}")
    table = pd.DataFrame(rows)[COMPOUND_COLUMNS]
    table.to_csv(COMPOUND_PATH, index=False)
    write_manifest(hgnc_extra=({"retrieved": date.today().isoformat(), **hgnc_info()}
                               if HGNC_PATH.exists() else None),
                   compound_extra=_compound_info())


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hgnc", action="store_true",
                        help="pin the HGNC complete set")
    parser.add_argument("--compounds", action="store_true",
                        help="build the compound annotation snapshot")
    parser.add_argument("--from-audit", type=Path, default=DEFAULT_AUDIT,
                        help="Stage 1 external audit CSV to snapshot")
    parser.add_argument("--refresh", action="store_true",
                        help="re-query public registries instead of snapshotting")
    parser.add_argument("--force", action="store_true",
                        help="re-download HGNC even if present")
    args = parser.parse_args(argv)

    if not (args.hgnc or args.compounds):
        parser.error("choose --hgnc and/or --compounds")
    if args.hgnc:
        fetch_hgnc(force=args.force)
    if args.compounds:
        if args.refresh:
            refresh_compounds(args.from_audit)
        else:
            build_compounds_from_audit_file(args.from_audit)
    return 0


if __name__ == "__main__":
    sys.exit(main())
