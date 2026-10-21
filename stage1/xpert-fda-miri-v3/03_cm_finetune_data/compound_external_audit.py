#!/usr/bin/env python
"""External SMILES + pharmacology audit of the pert_data CM compounds.

The `cm_drugs.html` report used to count only the structures the datasets
themselves shipped. Most GEO deposits ship names only, so that number says
nothing about whether a usable structure *exists*. This module asks the
authoritative public registries instead and records how each answer was
cross-checked:

  1. PubChem PUG-REST  — name (or dataset InChIKey) -> CID + SMILES + InChIKey.
  2. UniChem           — InChIKey -> ChEMBL / DrugBank / ChEBI / KEGG IDs.
  3. ChEMBL API        — molecule record (canonical SMILES, max phase, ATC),
                         mechanism endpoint (MoA + target).
  4. PubChem PUG-View  — MeSH Pharmacological Classification.
  5. RxNav             — WHO ATC classes by drug name (fallback).

Verification tiers
------------------
  verified     : the same structure (InChIKey connectivity block, first 14
                 chars) is registered in >= 2 independent sources — PubChem
                 plus at least one of ChEMBL / DrugBank / ChEBI / KEGG — or
                 the dataset-provided InChIKey agrees with PubChem.
  single       : exactly one registry yields the structure (still a legitimate
                 SMILES, but not independently cross-checked).
  conflict     : registries disagree on the connectivity block.
  unresolved   : no authoritative name/structure match.

Outputs (under 03_cm_finetune_data/output/)
  compound_inventory.csv        one row per dataset-native compound entry
  compound_external_audit.csv   one row per unique lookup name
  compound_external_cache.json  raw API responses, so reruns are free

Usage:
    python 03_cm_finetune_data/compound_external_audit.py
"""

from __future__ import annotations

import glob
import gzip
import json
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
WORKFLOW_ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(WORKFLOW_ROOT / "_lib"))
from utils import CM_DATA_ROOT  # noqa: E402

OUT_DIR = HERE / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CACHE_PATH = OUT_DIR / "compound_external_cache.json"
INVENTORY_PATH = OUT_DIR / "compound_inventory.csv"
AUDIT_PATH = OUT_DIR / "compound_external_audit.csv"

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUG_VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
UNICHEM = "https://www.ebi.ac.uk/unichem/rest"
RXNAV = "https://rxnav.nlm.nih.gov/REST"

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "arc-cm-compound-audit/1.0"})

# Registries in UniChem that count as an independent structural cross-check.
STRUCTURAL_XREFS = {"chembl", "drugbank", "chebi", "kegg", "pubchem"}


# ---------------------------------------------------------------------------
# HTTP with polite retries + on-disk cache
# ---------------------------------------------------------------------------
def _cache_load():
    if CACHE_PATH.exists():
        try:
            return json.loads(CACHE_PATH.read_text())
        except json.JSONDecodeError:
            return {}
    return {}


CACHE = _cache_load()


def _cache_save():
    CACHE_PATH.write_text(json.dumps(CACHE, ensure_ascii=False, indent=1))


_RATE_LOCK = threading.Lock()
_LAST_CALL = {}


def _rate_limit(host: str, min_interval: float):
    """Simple per-host token spacing; keeps the public APIs happy."""
    while True:
        with _RATE_LOCK:
            now = time.time()
            if now - _LAST_CALL.get(host, 0.0) >= min_interval:
                _LAST_CALL[host] = now
                return
            wait = min_interval - (now - _LAST_CALL.get(host, 0.0))
        time.sleep(min(wait, 0.2))


def _get_json(url, params=None, tries=6, pause=0.2):
    key = url + ("?" + json.dumps(params, sort_keys=True) if params else "")
    if key in CACHE:
        cached = CACHE[key]
        return cached if isinstance(cached, (dict, list)) else None
    host = url.split("/")[2]
    _rate_limit(host, 0.40 if "pubchem" in host else 0.15)
    last = None
    for attempt in range(tries):
        try:
            r = SESSION.get(url, params=params, timeout=60)
            if r.status_code == 404:
                # A real "not found" is cached; a transient failure is not.
                CACHE[key] = None
                return None
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(2.5 * (attempt + 1))
                continue
            r.raise_for_status()
            data = r.json()
            if not isinstance(data, (dict, list)):
                # Some endpoints occasionally return a bare JSON string.
                return None
            CACHE[key] = data
            time.sleep(pause)
            return data
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2.5 * (attempt + 1))
    print(f"  [warn] {url}: {last}", file=sys.stderr)
    # Deliberately do NOT cache a transient failure: a rerun retries it.
    return None


# ---------------------------------------------------------------------------
# Name normalisation / query variants
# ---------------------------------------------------------------------------
SALT_SUFFIXES = [
    "hydrochloride", "hydrobromide", "hydroiodide", "sulfate", "sulphate",
    "phosphate", "maleate", "fumarate", "hydrogensulfate", "tosylate",
    "mesylate", "citrate", "oxalate", "sodium", "potassium", "calcium",
    "dihydrate", "monohydrate", "hydrate", "acetate", "tartrate", "besylate",
    "metilsulfate", "saccharide", "diacetate",
]


def clean_name(name: str) -> str:
    """Whitespace/punctuation cleanup only (no chemistry decisions)."""
    n = str(name).replace("_", " ").strip()
    n = re.sub(r"\s+", " ", n)
    return n


# Dataset spellings that do not resolve under a literal PubChem/ChEMBL name
# query.  Each entry is a purely lexical/alias correction (the dataset token
# and the registry token denote the same substance); no structure is inferred.
QUERY_OVERRIDES = {
    "3 Azido 3 Deoxythymidine": [
        "3'-Azido-3'-deoxythymidine", "3-Azido-3-deoxythymidine", "Zidovudine"],
    "5 Fluorouracil": ["5-Fluorouracil"],
    "Acetylsalicylic acid": ["Aspirin", "Acetylsalicylic acid"],
    "D,1 Sotalol": ["Sotalol", "D,L-Sotalol"],
    "Clorophene": ["Chlorophene", "Clorophene"],
    "Glybenclamide": ["Glibenclamide", "Glyburide"],
    "Allura Red C.I.16035": ["Allura Red AC", "Allura Red"],
    "FD&C Green No. 3": ["Fast Green FCF", "FD&C Green No. 3"],
    "tert-Butylhydro-quinone (TBHQ)": ["tert-Butylhydroquinone", "TBHQ"],
}


# Structural fallbacks for the handful of names whose live PubChem PUG-REST
# query was rate-limited (HTTP 503) during the audit.  Each structure is
# confirmed by two independent registries; the source is written to the
# `notes` column.  These are used ONLY when the automatic lookup is
# unresolved, so they never override a live result.
MANUAL_FALLBACKS = {
    "Tralkoxydim": {
        "pubchem_cid": 135492483,
        "inchikey": "DQFPEYARZIQXRM-LTGZKZEYSA-N",
        "smiles": "CCO/N=C(\\CC)C1=C(O)CC(c2c(C)cc(C)cc2C)CC1=O",
        "chembl_id": "CHEMBL60556",
        "title": "Tralkoxydim",
        "source": ("ChEBI:81948 (SMILES/InChIKey) + UniChem ChEBI<->ChEMBL "
                   "CHEMBL60556<->PubChem CID 135492483; PubChem PUG-REST "
                   "was 503 rate-limited during the audit"),
    },
    "Rhodamine 6G": {
        "pubchem_cid": 13806,
        "inchikey": "VYXSBFYARXAAKO-UHFFFAOYSA-N",
        "smiles": ("CC/N=c1\\cc2oc3cc(NCC)c(C)cc3c(-c3ccccc3C(=O)OCC)"
                   "c-2cc1C"),
        "chembl_id": "CHEMBL1185241",
        "title": "Rhodamine 6G",
        "source": ("PubChem CID 13806 + ChEMBL CHEMBL1185241; PubChem "
                   "PUG-REST was 503 rate-limited during the audit"),
    },
}


def query_variants(name: str) -> list[str]:
    """Names to try, most faithful first.  Purely lexical variants."""
    out: list[str] = []

    def add(x):
        x = re.sub(r"\s+", " ", str(x)).strip()
        if x and x not in out:
            out.append(x)

    n = clean_name(name)
    for q in QUERY_OVERRIDES.get(n, []):
        add(q)
    add(n)
    # D,1 / D,L / D,L- racemic prefix (dataset spelling)
    m = re.match(r"^(?:d\s*[,.]?\s*[1l]|dl)\s*[-–]?\s*(.+)$", n, re.I)
    if m:
        add("D,L-" + m.group(1))
        add(m.group(1))
    # trailing parenthetical / slash annotations
    if "(" in n:
        add(n.split("(")[0])
    if "/" in n:
        add(n.split("/")[0])
    # salt/solvate suffix -> parent
    low = n.lower()
    for suf in SALT_SUFFIXES:
        if low.endswith(" " + suf):
            add(n[: -len(suf) - 1])
            break
    return out


def _norm(s: str) -> str:
    return "".join(ch for ch in str(s).lower() if ch.isalnum())


# ---------------------------------------------------------------------------
# PubChem
# ---------------------------------------------------------------------------
def pubchem_by_name(name: str):
    cids = []
    for q in query_variants(name):
        data = _get_json(f"{PUG}/compound/name/{requests.utils.quote(q)}/cids/JSON")
        if isinstance(data, dict) and data.get("IdentifierList", {}).get("CID"):
            cids = data["IdentifierList"]["CID"]
            break
    if not cids:
        return {"found": False, "query": name}
    cid = int(cids[0])
    props = _get_json(
        f"{PUG}/compound/cid/{cid}/property/"
        "CanonicalSMILES,IsomericSMILES,InChIKey,MolecularFormula,"
        "MolecularWeight,Title/JSON")
    if (not isinstance(props, dict)
            or not props.get("PropertyTable", {}).get("Properties")):
        return {"found": False, "query": name, "cid": cid}
    p = props["PropertyTable"]["Properties"][0]
    smi = (p.get("SMILES") or p.get("CanonicalSMILES")
           or p.get("ConnectivitySMILES"))
    exact = True
    if len(cids) > 1:
        syn = _get_json(f"{PUG}/compound/cid/{cid}/synonyms/JSON")
        names = [x for x in (syn or {}).get("InformationList", {})
                 .get("Information", [{}])[0].get("Synonym", [])]
        exact = any(_norm(x) == _norm(name) for x in names)
    return {
        "found": True,
        "query": name,
        "cid": cid,
        "n_cids": len(cids),
        "name_exact": exact,
        "title": p.get("Title"),
        "smiles": smi,
        "isomeric_smiles": p.get("IsomericSMILES") or smi,
        "inchikey": p.get("InChIKey"),
        "formula": p.get("MolecularFormula"),
        "mw": p.get("MolecularWeight"),
        "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
    }


def pubchem_by_inchikey(key: str):
    data = _get_json(f"{PUG}/compound/inchikey/{key}/cids/JSON")
    if not isinstance(data, dict) or not data.get("IdentifierList", {}).get("CID"):
        return {"found": False, "query": key}
    cid = int(data["IdentifierList"]["CID"][0])
    # No `Title` here: this is the hot path for structure-backed datasets and
    # the base property payload is already cached.  Name-based lookups below
    # do request the title.
    props = _get_json(
        f"{PUG}/compound/cid/{cid}/property/"
        "CanonicalSMILES,IsomericSMILES,InChIKey,MolecularFormula,"
        "MolecularWeight/JSON")
    p = (props["PropertyTable"]["Properties"][0]
         if isinstance(props, dict) else {})
    smi = (p.get("SMILES") or p.get("CanonicalSMILES")
           or p.get("ConnectivitySMILES"))
    return {
        "found": True,
        "query": key,
        "cid": cid,
        "n_cids": 1,
        "name_exact": True,
        "title": p.get("Title"),
        "smiles": smi,
        "isomeric_smiles": p.get("IsomericSMILES") or smi,
        "inchikey": p.get("InChIKey"),
        "formula": p.get("MolecularFormula"),
        "mw": p.get("MolecularWeight"),
        "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
    }


def pubchem_mesh(cid: int):
    data = _get_json(
        f"{PUG_VIEW}/data/compound/{cid}/JSON",
        params={"heading": "MeSH Pharmacological Classification"})
    if not isinstance(data, dict):
        return []
    out = []

    def walk(sec):
        if isinstance(sec, dict):
            if sec.get("TOCHeading") == "MeSH Pharmacological Classification":
                for info in sec.get("Information", []):
                    # The class name is in `Name`; `Value` holds its definition.
                    s = info.get("Name")
                    if not s:
                        for swm in info.get("Value", {}).get(
                                "StringWithMarkup", []):
                            s = swm.get("String")
                            break
                    if s and s not in out:
                        out.append(s)
            for v in sec.values():
                walk(v)
        elif isinstance(sec, list):
            for v in sec:
                walk(v)

    walk(data.get("Record", {}))
    return out


# ---------------------------------------------------------------------------
# UniChem + ChEMBL
# ---------------------------------------------------------------------------
def unichem_xrefs(inchikey: str):
    data = _get_json(f"{UNICHEM}/verbose_inchikey/{inchikey}")
    if not isinstance(data, list):
        return {}
    out = {}
    for src in data:
        name = src.get("name")
        ids = src.get("src_compound_id") or []
        if name and ids:
            out.setdefault(name, []).extend(ids)
    return out


def _chembl_molecule(chembl_id: str):
    return _get_json(f"{CHEMBL}/molecule/{chembl_id}.json")


def _chembl_mechanisms(chembl_id: str):
    data = _get_json(f"{CHEMBL}/mechanism.json",
                     params={"molecule_chembl_id": chembl_id, "limit": 100})
    return data.get("mechanisms", []) if isinstance(data, dict) else []


def _chembl_atc_descriptions(chembl_id: str):
    data = _get_json(f"{CHEMBL}/drug.json",
                     params={"molecule_chembl_id": chembl_id, "limit": 5})
    out = {}
    for d in (data.get("drugs", []) if isinstance(data, dict) else []):
        for a in d.get("atc_classification") or []:
            if a.get("code"):
                out.setdefault(a["code"], a.get("description", ""))
    return out


def chembl_by_name(name: str):
    for q in query_variants(name):
        data = _get_json(
            f"{CHEMBL}/molecule.json",
            params={"molecule_synonyms__molecule_synonym__iexact": q,
                    "limit": 2})
        mols = (data or {}).get("molecules", [])
        if mols:
            return mols[0].get("molecule_chembl_id")
    # fuzzy fallback: accept only an exact normalised pref_name/synonym
    for q in query_variants(name):
        data = _get_json(f"{CHEMBL}/search.json", params={"q": q, "limit": 5})
        for m in (data or {}).get("molecules", []):
            if _norm(m.get("pref_name") or "") in (_norm(q), _norm(name)):
                return m.get("molecule_chembl_id")
    return None


def _short_atc(desc: str) -> str:
    """Keep the most specific ATC class (last hierarchy segment)."""
    parts = [p.strip() for p in str(desc).split(":") if p.strip()]
    return parts[-1] if parts else str(desc)


def chembl_bundle(chembl_id: str):
    mol = _chembl_molecule(chembl_id)
    if not isinstance(mol, dict):
        return None
    mech = _chembl_mechanisms(chembl_id)
    atc_codes = mol.get("atc_classifications") or []
    atc_desc = _chembl_atc_descriptions(chembl_id) if atc_codes else {}
    return {
        "chembl_id": chembl_id,
        "pref_name": mol.get("pref_name"),
        "molecule_type": mol.get("molecule_type"),
        "max_phase": mol.get("max_phase"),
        "smiles": (mol.get("molecule_structures") or {}).get(
            "canonical_smiles"),
        "atc": [{"code": c, "description": _short_atc(atc_desc.get(c, ""))}
                for c in atc_codes],
        "moa": sorted({m.get("mechanism_of_action")
                       for m in mech if m.get("mechanism_of_action")}),
        "targets": sorted({m.get("target_chembl_id") for m in mech
                           if m.get("target_chembl_id")}),
        "url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}/",
    }


def rxnav_atc(name: str):
    data = _get_json(f"{RXNAV}/rxclass/class/byDrugName.json",
                     params={"drugName": name, "relaSource": "ATC"})
    if not isinstance(data, dict):
        return []
    info = data.get("rxclassDrugInfoList")
    if not isinstance(info, dict):
        return []
    out = []
    for d in info.get("rxclassDrugInfo") or []:
        if not isinstance(d, dict):
            continue
        c = d.get("rxclassMinConceptItem", {})
        if c.get("classId") and c.get("classId") not in [x["code"] for x in out]:
            out.append({"code": c["classId"], "description": c.get("className", "")})
    return out


# ---------------------------------------------------------------------------
# RDKit helpers
# ---------------------------------------------------------------------------
def _rdkit_key(smiles: str):
    from rdkit import Chem
    if not smiles:
        return None
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    try:
        return Chem.MolToInchiKey(m)
    except Exception:  # noqa: BLE001
        return None


def _parent_key(smiles: str):
    """InChIKey of the largest organic fragment (salt/solvate -> parent)."""
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    if not smiles:
        return None
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    try:
        parent = rdMolStandardize.LargestFragmentChooser().choose(m)
        return Chem.MolToInchiKey(parent)
    except Exception:  # noqa: BLE001
        return None


def _prefix(key):
    return (key or "")[:14]


# ---------------------------------------------------------------------------
# Compound inventory
# ---------------------------------------------------------------------------
DTOXS_META = (CM_DATA_ROOT / "dtoxs_primary" / "Metadata"
              / "Small_Molecule_Metadata.txt")
DTOXS_ANTIBODY = CM_DATA_ROOT / "dtoxs_primary" / "Metadata" / "Antibody_Metadata.txt"
DTOXS_PROTEIN = CM_DATA_ROOT / "dtoxs_primary" / "Metadata" / "Protein_Metadata.txt"

DTOXS_BIOLOGIC_CODES = {
    # code: (entity type, display name, source)
    "TRS": ("antibody", "Trastuzumab", "Antibody_Metadata"),
    "IGF": ("protein", "IGF-1", "Protein_Metadata"),
    "TNF": ("protein", "TNF-alpha", "Protein_Metadata"),
    "END": ("peptide", "Endothelin-1", "Small_Molecule_Metadata (EDN-1)"),
}


def _dtoxs_primary_design_codes():
    codes = set()
    for f in glob.glob(str(CM_DATA_ROOT / "dtoxs_primary" / "Data" / "Level3"
                          / "**" / "PromoCell-DGE-Experiment-Design-*.tsv"),
                       recursive=True):
        codes |= set(pd.read_csv(f, sep="\t")["State"].astype(str))
    return {c for c in codes if not c.startswith("CTRL")}


def _ipsc_design_codes():
    f = sorted(glob.glob(str(CM_DATA_ROOT / "dtoxs_ipsc"
                            / "GSE217423*Design*.tsv.gz")))[0]
    codes = set(pd.read_csv(f, sep="\t")["State"].astype(str))
    return {c for c in codes if c != "CTRL"}


def _ipsc_names():
    t = gzip.open(CM_DATA_ROOT / "dtoxs_ipsc"
                  / "GSE217423_series_matrix.txt.gz", "rt").read()
    codes = re.findall(r'"([A-Z0-9-]+) MSN[^"]*"', t)
    names = re.findall(r'"drug name: ([^"]*)"', t)
    out = {}
    for c, n in zip(codes, names):
        out.setdefault(c, n)
    return out


def _gse244740_names():
    t = gzip.open(CM_DATA_ROOT / "gse244740"
                  / "GSE244740_series_matrix.txt.gz", "rt").read()
    line = [l for l in t.splitlines() if l.startswith("!Sample_title")][0]
    titles = re.findall(r'"([^"]*)"', line)
    names = set()
    for s in titles:
        m = re.match(r"(.+?), (Con\d+|NoTreatment), (?:Optimization|Screening) plate", s)
        if m:
            names.add(m.group(1))
    return sorted(names - {"DMSO", "NoTreatment"})


def _gse247541_names():
    c = pd.read_csv(CM_DATA_ROOT / "gse247541" / "GSE247541_raw_counts.csv.gz",
                    index_col=0, nrows=0)
    return sorted({col.rsplit("_", 1)[0] for col in c.columns
                   if not col.startswith("DMSO")})


def _gse262419_names():
    h = pd.read_csv(CM_DATA_ROOT / "gse262419" / "GSE262419_hash.csv.gz")
    names = {str(x).strip() for x in h["Chemical_name"]}
    return sorted(names - {"VEH", "Media", "Blank"})


def build_inventory() -> pd.DataFrame:
    meta = pd.read_csv(DTOXS_META, sep="\t", dtype=str)
    by_code = {r["sm_center_batch_id"]: r for _, r in meta.iterrows()}
    rows = []

    # --- dtoxs_primary -----------------------------------------------------
    for code in sorted(_dtoxs_primary_design_codes()):
        if code in DTOXS_BIOLOGIC_CODES:
            typ, name, src = DTOXS_BIOLOGIC_CODES[code]
            row = {"dataset": "dtoxs_primary", "compound_code": code,
                   "compound_name": name, "entity_type": typ,
                   "dataset_smiles": None, "dataset_inchikey": None,
                   "note": f"code {code}; {src}"}
        elif code == "SUN+":
            r = by_code.get("SUN")
            row = {"dataset": "dtoxs_primary", "compound_code": code,
                   "compound_name": "Sunitinib", "entity_type": "small_molecule",
                   "dataset_smiles": r["sm_smiles_parent"] if r is not None else None,
                   "dataset_inchikey": r["sm_inchi_key_parent"] if r is not None else None,
                   "note": "code SUN+; duplicate of SUN (same compound)"}
        elif code in by_code:
            r = by_code[code]
            name = r["sm_name"]
            note = ""
            if name.upper().startswith("ZINC"):
                alt = (r.get("sm_alternative_name") or "").split(";")[0].strip()
                if alt:
                    note = f"metadata sm_name={name}; displayed as {alt}"
                    name = alt
            row = {"dataset": "dtoxs_primary", "compound_code": code,
                   "compound_name": name, "entity_type": "small_molecule",
                   "dataset_smiles": r["sm_smiles_parent"],
                   "dataset_inchikey": r["sm_inchi_key_parent"],
                   "note": note}
        elif code == "DIC":
            row = {"dataset": "dtoxs_primary", "compound_code": code,
                   "compound_name": "Diclofenac", "entity_type": "small_molecule",
                   "dataset_smiles": None, "dataset_inchikey": None,
                   "note": "code DIC; no metadata row (name inferred as diclofenac)"}
        else:
            row = {"dataset": "dtoxs_primary", "compound_code": code,
                   "compound_name": code, "entity_type": "unknown",
                   "dataset_smiles": None, "dataset_inchikey": None,
                   "note": "unmapped code"}
        rows.append(row)

    # --- dtoxs_ipsc --------------------------------------------------------
    ipsc_names = _ipsc_names()
    for code in sorted(_ipsc_design_codes()):
        if code in ("BEV", "CTX", "RTX", "TRS"):
            rows.append({
                "dataset": "dtoxs_ipsc", "compound_code": code,
                "compound_name": ipsc_names.get(code, code),
                "entity_type": "antibody", "dataset_smiles": None,
                "dataset_inchikey": None,
                "note": "monoclonal antibody (series matrix)"})
            continue
        r = by_code.get(code)
        rows.append({
            "dataset": "dtoxs_ipsc", "compound_code": code,
            "compound_name": ipsc_names.get(code, code),
            "entity_type": "small_molecule",
            "dataset_smiles": r["sm_smiles_parent"] if r is not None else None,
            "dataset_inchikey": r["sm_inchi_key_parent"] if r is not None else None,
            "note": "" if r is not None else "no DToxS metadata row"})

    # --- gse244740 ---------------------------------------------------------
    for name in _gse244740_names():
        rows.append({"dataset": "gse244740", "compound_code": "",
                     "compound_name": name, "entity_type": "small_molecule",
                     "dataset_smiles": None, "dataset_inchikey": None,
                     "note": "name-only GEO deposit"})
    # --- gse247541 ---------------------------------------------------------
    for name in _gse247541_names():
        rows.append({"dataset": "gse247541", "compound_code": "",
                     "compound_name": name, "entity_type": "small_molecule",
                     "dataset_smiles": None, "dataset_inchikey": None,
                     "note": "name-only GEO deposit"})
    # --- gse262419 ---------------------------------------------------------
    for name in _gse262419_names():
        if name.upper() == "ET-1":
            rows.append({"dataset": "gse262419", "compound_code": "",
                         "compound_name": "ET-1", "entity_type": "peptide",
                         "dataset_smiles": None, "dataset_inchikey": None,
                         "note": "endothelin-1 peptide; not a small molecule"})
        else:
            rows.append({"dataset": "gse262419", "compound_code": "",
                         "compound_name": name, "entity_type": "small_molecule",
                         "dataset_smiles": None, "dataset_inchikey": None,
                         "note": "name-only TempO-Seq chemical library"})

    inv = pd.DataFrame(rows)
    inv["lookup_name"] = inv["compound_name"].map(clean_name)
    return inv


# ---------------------------------------------------------------------------
# Per-compound external lookup
# ---------------------------------------------------------------------------
def audit_one(rec: dict) -> dict:
    name = rec["lookup_name"]
    ds_smi = rec.get("dataset_smiles")
    ds_key = _rdkit_key(ds_smi) if isinstance(ds_smi, str) and ds_smi.strip() else None
    meta_key = rec.get("dataset_inchikey")
    meta_key = (meta_key if isinstance(meta_key, str) and meta_key.strip()
                else None)

    notes = []
    target_key = ds_key or meta_key
    # The DToxS metadata carries both a SMILES and an InChIKey. If they
    # disagree, the curated InChIKey wins and the discrepancy is recorded.
    if ds_key and meta_key and _prefix(ds_key) != _prefix(meta_key):
        notes.append("metadata InChIKey 与 SMILES 不一致（按 InChIKey）")
        target_key = meta_key
    target_prefix = _prefix(target_key)
    target_parent = _prefix(_parent_key(ds_smi)) if ds_smi else None

    # --- structure resolution ---------------------------------------------
    pc = None
    if target_key:
        pc = pubchem_by_inchikey(target_key)
        if not pc.get("found"):
            pc = None
    if pc is None:
        pcn = pubchem_by_name(name)
        if pcn.get("found"):
            if target_prefix and _prefix(pcn["inchikey"]) != target_prefix:
                notes.append("名称查询命中不同结构，未采用")
                pc = None
            else:
                pc = pcn
                if target_key:
                    notes.append("数据集结构键未在 PubChem 命中，名称命中且键一致")
    pc = pc or {"found": False}
    pub_key = pc.get("inchikey")
    pub_smi = pc.get("smiles")
    pub_parent = _prefix(_parent_key(pub_smi)) if pub_smi else None

    # --- ChEMBL / UniChem --------------------------------------------------
    xrefs = unichem_xrefs(pub_key) if pub_key else {}
    chembl_id = (xrefs.get("chembl") or [None])[0]
    if not chembl_id and pub_parent:
        # salt/solvate form not registered: try the parent InChIKey
        parent_key = _parent_key(pub_smi)
        if parent_key:
            pxrefs = unichem_xrefs(parent_key)
            xrefs = {**pxrefs, **xrefs}
            chembl_id = (pxrefs.get("chembl") or [None])[0]
    alias_mismatch = False
    if not chembl_id:
        by_name = chembl_by_name(name)
        if by_name:
            cb_try = chembl_bundle(by_name)
            try_key = _rdkit_key((cb_try or {}).get("smiles"))
            try_parent = _prefix(_parent_key((cb_try or {}).get("smiles")))
            matches = ((pub_key and _prefix(try_key) == _prefix(pub_key))
                       or (pub_parent and try_parent == pub_parent)
                       or (target_prefix and _prefix(try_key) == target_prefix))
            if matches:
                chembl_id = by_name
            else:
                alias_mismatch = True
                notes.append("名称在 ChEMBL 命中不同结构，已弃用其分类")
    cb = chembl_bundle(chembl_id) if chembl_id else None
    cb_key = _rdkit_key((cb or {}).get("smiles"))
    cb_parent = _prefix(_parent_key((cb or {}).get("smiles")))

    # ChEMBL-only hit -> ask PubChem by the ChEMBL structure.
    if not pc.get("found") and cb_key:
        pc = pubchem_by_inchikey(cb_key)
        pub_key = pc.get("inchikey") if pc.get("found") else None
        pub_smi = pc.get("smiles") if pc.get("found") else None

    # --- verification tier -------------------------------------------------
    full_match = bool(pub_key and cb_key
                      and _prefix(pub_key) == _prefix(cb_key))
    parent_match = bool(pub_parent and cb_parent and pub_parent == cb_parent)
    if full_match or parent_match:
        tier = "verified"
        if parent_match and not full_match:
            notes.append("盐/母体形式匹配")
    elif pub_key and any(x in xrefs for x in STRUCTURAL_XREFS - {"pubchem"}):
        tier = "verified"
    elif pub_key or cb_key:
        tier = "single"
    else:
        tier = "unresolved"
    if (pub_key and cb_key and _prefix(pub_key) != _prefix(cb_key)
            and not parent_match):
        tier = "conflict"
        notes.append("PubChem 与 ChEMBL 结构不一致")

    # --- pharmacology ------------------------------------------------------
    mesh = pubchem_mesh(pc["cid"]) if pc.get("found") else []
    atc = (cb or {}).get("atc") or []
    moa = (cb or {}).get("moa") or []
    if not atc:
        # RxNav only on the resolved structure's registry name, never on an
        # unverified dataset alias (e.g. DToxS "Amiodarone" = ZINC150345845).
        rx_name = (cb or {}).get("pref_name") or pc.get("title")
        if isinstance(rx_name, str) and rx_name.strip():
            atc = rxnav_atc(rx_name)
    has_pharm = bool(atc or moa or mesh)

    smiles = (cb or {}).get("smiles") or pub_smi
    if smiles is None and ds_smi:
        smiles = ds_smi

    result = {
        "lookup_name": name,
        "entity_type": rec["entity_type"],
        "dataset_inchikey": meta_key,
        "dataset_smiles": ds_smi,
        "pubchem_cid": pc.get("cid"),
        "pubchem_title": pc.get("title"),
        "pubchem_smiles": pub_smi,
        "pubchem_inchikey": pub_key,
        "pubchem_url": pc.get("url"),
        "pubchem_name_exact": pc.get("name_exact"),
        "chembl_id": (cb or {}).get("chembl_id"),
        "chembl_pref_name": (cb or {}).get("pref_name"),
        "chembl_type": (cb or {}).get("molecule_type"),
        "chembl_max_phase": (cb or {}).get("max_phase"),
        "chembl_smiles": (cb or {}).get("smiles"),
        "chembl_url": (cb or {}).get("url"),
        "unichem_sources": ",".join(sorted(k for k in xrefs
                                           if k in STRUCTURAL_XREFS)),
        "smiles": smiles,
        "parent_inchikey": _parent_key(smiles),
        "verification": tier,
        "alias_mismatch": alias_mismatch,
        "atc": "; ".join(f"{a['code']} {a.get('description', '')}".strip()
                         for a in atc),
        "moa": "; ".join(moa),
        "moa_targets": "; ".join((cb or {}).get("targets") or []),
        "mesh_classes": "; ".join(mesh),
        "has_pharmacology": has_pharm,
        "notes": "; ".join(notes),
    }

    fb = MANUAL_FALLBACKS.get(name)
    if fb and result["verification"] in ("unresolved", "error"):
        fbcb = chembl_bundle(fb["chembl_id"]) if fb.get("chembl_id") else None
        fatc = (fbcb or {}).get("atc") or []
        fmoa = (fbcb or {}).get("moa") or []
        result.update({
            "pubchem_cid": fb["pubchem_cid"],
            "pubchem_title": fb["title"],
            "pubchem_smiles": fb["smiles"],
            "pubchem_inchikey": fb["inchikey"],
            "pubchem_url": ("https://pubchem.ncbi.nlm.nih.gov/compound/"
                            f"{fb['pubchem_cid']}"),
            "chembl_id": (fbcb or {}).get("chembl_id") or fb.get("chembl_id"),
            "chembl_pref_name": (fbcb or {}).get("pref_name"),
            "chembl_type": (fbcb or {}).get("molecule_type"),
            "chembl_max_phase": (fbcb or {}).get("max_phase"),
            "chembl_smiles": (fbcb or {}).get("smiles"),
            "chembl_url": (fbcb or {}).get("url"),
            "smiles": fb["smiles"],
            "parent_inchikey": _parent_key(fb["smiles"]),
            "verification": "verified",
            "atc": "; ".join(f"{a['code']} {a.get('description', '')}".strip()
                             for a in fatc),
            "moa": "; ".join(fmoa),
            "moa_targets": "; ".join((fbcb or {}).get("targets") or []),
            "has_pharmacology": bool(fatc or fmoa or mesh),
            "notes": "; ".join([n for n in notes if n] + [fb["source"]]),
        })

    return result


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--refresh-none", action="store_true",
                    help="drop cached nulls (404/transient failures) and refetch")
    args = ap.parse_args()

    if args.refresh_none:
        stale = [k for k, v in CACHE.items() if v is None]
        for k in stale:
            del CACHE[k]
        _cache_save()
        print(f"Dropped {len(stale)} cached nulls")

    inv = build_inventory()
    inv.to_csv(INVENTORY_PATH, index=False)

    # One lookup per (name, dataset InChIKey); shared names across datasets
    # reuse the same result.
    recs = (inv[inv["entity_type"] == "small_molecule"]
            .drop_duplicates(["lookup_name", "dataset_inchikey"])
            .to_dict("records"))
    print(f"Inventory: {len(inv)} entries, {len(recs)} unique lookup keys")

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(audit_one, r): r for r in recs}
        for i, fut in enumerate(as_completed(futs), 1):
            rec = futs[fut]
            try:
                res = fut.result()
            except Exception as exc:  # noqa: BLE001
                print(f"  [error] {rec['lookup_name']}: {exc}", file=sys.stderr)
                res = {"lookup_name": rec["lookup_name"],
                       "entity_type": rec["entity_type"],
                       "dataset_inchikey": rec.get("dataset_inchikey"),
                       "verification": "error"}
            results.append(res)
            if i % 25 == 0:
                print(f"  {i}/{len(recs)}")
                _cache_save()

    _cache_save()
    out = pd.DataFrame(results).sort_values("lookup_name")
    out.to_csv(AUDIT_PATH, index=False)
    print(f"Wrote {INVENTORY_PATH}")
    print(f"Wrote {AUDIT_PATH}")
    print(out["verification"].value_counts().to_string())


if __name__ == "__main__":
    main()
