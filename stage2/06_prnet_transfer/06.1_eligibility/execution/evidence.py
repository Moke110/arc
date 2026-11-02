"""Evidence collection for the temp_env eligibility audit.

Everything this module adds to a decision comes from a documented source and is
cached with its URL and retrieval timestamp so a rerun is offline and auditable:

* PubChem PUG-REST (InChIKey -> CID) and PUG-View (full compound record),
* the Broad Institute Repurposing Hub tables in the local Biomni data lake,
* BindingDB 2024-09 ligand-target pairs (local TSV, streamed, never loaded whole),
* the frozen ``compounds.csv`` (identity + annotation under audit, not evidence).

Network access goes through ``curl``: the shared ``arc`` environment's Python
cannot verify the site certificate (self-signed certificate in the chain) and
ChEMBL is unreachable, which is recorded by ``probe_chembl`` rather than papered
over.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import subprocess
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

CACHE_FORMAT = "jsonl-index-v1"
_CACHE_LOCK = threading.Lock()
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUGVIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"
DATA_LAKE = Path(
    os.environ.get(
        "ARC_BIOMNI_DATALAKE",
        "/home/chang/.local/share/biomni/data/biomni_data/data_lake",
    )
)
DEFAULT_BINDINGDB = DATA_LAKE / "BindingDB_All_202409.tsv"
BROAD_PHASE = DATA_LAKE / "broad_repurposing_hub_phase_moa_target_info.parquet"
BROAD_MOLECULE = DATA_LAKE / "broad_repurposing_hub_molecule_with_smiles.parquet"
CHEMBL_PROBE = "https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL25.json"

# Strings PubChem repeats in section banners rather than in a compound record.
BOILERPLATE = (
    "ClinicalTrials.gov is a registry",
    "For more ",
)


class EvidenceError(RuntimeError):
    """A documented source could not be read and the build must stop."""


def find_repo_root(start=None):
    """Return the closest ancestor of *start* that contains ``.git``."""
    path = Path(start).resolve() if start is not None else Path(__file__).resolve()
    if path.is_file():
        path = path.parent
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    raise EvidenceError(f"no .git directory at or above {path}")


ROOT = find_repo_root()
STEP = Path(__file__).resolve().parent.parent
EXECUTION = Path(__file__).resolve().parent
OUTPUT = STEP / "output"


def now_iso():
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_timestamp(path):
    """ISO timestamp of a local evidence file (its mtime), stable across rebuilds.

    Local sources (Broad Repurposing Hub tables, BindingDB TSV) are evidence only
    at the bytes on disk, so their retrieval timestamp is the file mtime rather
    than the build time: an offline rebuild then reproduces the same record.
    """
    try:
        return (
            datetime.fromtimestamp(Path(path).stat().st_mtime, tz=timezone.utc)
            .astimezone()
            .isoformat(timespec="seconds")
        )
    except OSError:
        return ""


class JsonlCache:
    """Append-only response cache: a JSONL index plus one gzip body per URL.

    The index line keeps ``url``, ``fetched_at``, ``http_status`` and the body
    file name, so a rerun is offline and every value in the audit table can be
    traced to a cached response without re-reading multi-megabyte bodies.
    """

    def __init__(self, path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.bodies = self.path.parent / (self.path.name.replace(".jsonl", "") + "_bodies")
        self.bodies.mkdir(parents=True, exist_ok=True)
        self._records = None

    def load(self):
        if self._records is not None:
            return self._records
        records = {}
        if self.path.exists():
            with open(self.path, encoding="utf-8") as handle:
                for line in handle:
                    line = line.strip()
                    if not line:
                        continue
                    record = json.loads(line)
                    records[record["url"]] = record
        self._records = records
        return records

    def store(self, record):
        name = hashlib.sha1(record["url"].encode("utf-8")).hexdigest() + ".json.gz"
        meta = {key: value for key, value in record.items() if key != "body"}
        meta["body_file"] = name
        with _CACHE_LOCK:
            with gzip.open(self.bodies / name, "wt", encoding="utf-8") as handle:
                handle.write(record.get("body") or "")
            with open(self.path, "a", encoding="utf-8") as handle:
                handle.write(json.dumps(meta, sort_keys=True) + "\n")
            self.load()[record["url"]] = meta
        return meta

    def read(self, meta):
        body = ""
        if meta and meta.get("body_file"):
            try:
                with gzip.open(self.bodies / meta["body_file"], "rt", encoding="utf-8") as fh:
                    body = fh.read()
            except OSError:
                body = ""
        return {**(meta or {}), "body": body}

    def get(self, url, allow_network=True, timeout=120, retries=3):
        meta = self.load().get(url)
        record = self.read(meta)
        if record.get("body"):
            return record
        if not allow_network:
            return record if meta else None
        fresh = fetch_json(url, timeout=timeout, retries=retries)
        self.store(fresh)
        return fresh


def fetch_json(url, timeout=120, retries=3):
    """Fetch one URL with curl; return the JSONL cache record."""
    last = {"url": url, "fetched_at": now_iso(), "http_status": 0, "body": ""}
    for attempt in range(retries):
        proc = subprocess.run(
            [
                "curl",
                "-sS",
                "--max-time",
                str(timeout),
                "-w",
                "\n__HTTP__%{http_code}",
                url,
            ],
            capture_output=True,
            text=True,
        )
        body, marker, code = proc.stdout.rpartition("\n__HTTP__")
        status = int(code) if code.strip().isdigit() else 0
        last = {
            "url": url,
            "fetched_at": now_iso(),
            "http_status": status,
            "body": body if status == 200 else "",
            "stderr": proc.stderr.strip()[:500] if status != 200 else "",
        }
        if status == 200 and body.strip():
            try:
                json.loads(body)
                return last
            except json.JSONDecodeError:
                last["stderr"] = "response was not JSON"
        time.sleep(1.5 * (attempt + 1))
    return last


def probe_chembl(cache_dir, allow_network=True):
    """Record the ChEMBL availability attempt (2026-11-02: unreachable)."""
    cache = JsonlCache(Path(cache_dir) / "chembl_probe.jsonl")
    record = cache.get(CHEMBL_PROBE, allow_network=allow_network, retries=2)
    return {
        "url": CHEMBL_PROBE,
        "http_status": record.get("http_status") if record else None,
        "fetched_at": record.get("fetched_at") if record else None,
        "used_as_evidence": False,
    }


# --------------------------------------------------------------------------
# PubChem
# --------------------------------------------------------------------------
def pubchem_cids(inchikey, cache_dir, allow_network=True):
    cache = JsonlCache(Path(cache_dir) / "pubchem_cids.jsonl")
    record = cache.get(
        f"{PUBCHEM}/compound/inchikey/{inchikey}/cids/JSON", allow_network=allow_network
    )
    cids = []
    if record and record.get("http_status") == 200:
        try:
            cids = json.loads(record["body"])["IdentifierList"]["CID"]
        except (KeyError, json.JSONDecodeError):
            cids = []
    return {
        "cids": [int(cid) for cid in cids],
        "http_status": record.get("http_status") if record else None,
        "fetched_at": record.get("fetched_at") if record else None,
        "url": record.get("url") if record else None,
    }


def pubchem_record(cid, cache_dir, allow_network=True):
    cache = JsonlCache(Path(cache_dir) / "pubchem_pugview.jsonl")
    return cache.get(
        f"{PUGVIEW}/data/compound/{cid}/JSON", allow_network=allow_network
    )


def _value_strings(info):
    value = info.get("Value") or {}
    out = []
    for item in value.get("StringWithMarkup") or []:
        text = (item or {}).get("String")
        if text:
            out.append(str(text))
    if not out:
        number = value.get("Number")
        if isinstance(number, list) and number:
            out.append(str(number[0]))
        elif isinstance(number, (int, float)):
            out.append(str(number))
        elif value.get("String"):
            out.append(str(value["String"]))
    return out


def _walk_sections(section, path, leaves):
    heading = section.get("TOCHeading") or ""
    here = f"{path}/{heading}" if path else heading
    for child in section.get("Section") or []:
        _walk_sections(child, here, leaves)
    for info in section.get("Information") or []:
        for text in _value_strings(info):
            leaves.append((here, text))


def _clean(strings, limit=None):
    seen, out = set(), []
    for text in strings:
        text = re.sub(r"\s+", " ", str(text)).strip()
        if not text or any(flag in text for flag in BOILERPLATE):
            continue
        key = text.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(text)
        if limit and len(out) >= limit:
            break
    return out


def extract_pubchem(record_body):
    """Extract the documented evidence fields from one PUG-View record."""
    try:
        record = json.loads(record_body)["Record"]
    except (KeyError, json.JSONDecodeError):
        return {}
    leaves = []
    for section in record.get("Section") or []:
        _walk_sections(section, "", leaves)

    def paths_ending(suffix):
        return [text for path, text in leaves if path.endswith(suffix)]

    def under(prefix):
        return [text for path, text in leaves if path.startswith(prefix)]

    max_phase = _clean(paths_ending("Drug Development Summary/Max Phase"))
    return {
        "pubchem_title": record.get("RecordTitle") or "",
        "pubchem_max_phase": max_phase[0] if max_phase else "",
        "atc_pubchem": _clean(paths_ending("Pharmacology and Biochemistry/ATC Code")),
        "drug_classes": _clean(paths_ending("Drug and Medication Information/Drug Classes")),
        "drug_indications": _clean(
            paths_ending("Drug and Medication Information/Drug Indication")
        ),
        "therapeutic_uses": _clean(
            paths_ending("Drug and Medication Information/Therapeutic Uses"), limit=20
        ),
        "mechanism_texts": _clean(
            paths_ending("Pharmacology and Biochemistry/Mechanism of Action"), limit=20
        ),
        "drug_targets": _clean(
            paths_ending("Drug and Medication Information/Drug Targets"), limit=40
        ),
        "use_classification": _clean(
            paths_ending("Use and Manufacturing/Uses/Use Classification")
        ),
        "industry_uses": _clean(paths_ending("Use and Manufacturing/Uses/Industry Uses")),
        # The bare HSDB/EPA "Uses" entries are the sentences that document what a
        # compound is actually used for ("Used as broad-spectrum fungicide ...").
        # Without them a pesticide with a mechanism text but no Use Classification
        # string looks undocumented; v1.3 reads them as use evidence.
        "uses_texts": _clean(
            paths_ending("Use and Manufacturing/Uses"), limit=25
        ),
        "fda_approved_drugs": _clean(
            paths_ending("Drug and Medication Information/FDA Approved Drugs")
        ),
        "who_essential_medicines": _clean(
            paths_ending("Drug and Medication Information/WHO Essential Medicines")
        ),
        "synonyms": _clean(under("Names and Identifiers/Synonyms"), limit=400),
    }


# --------------------------------------------------------------------------
# Broad Repurposing Hub
# --------------------------------------------------------------------------
def normalize_name(value):
    return re.sub(r"[^a-z0-9]", "", str(value or "").lower())


class BroadHub:
    """Lookup by full InChIKey, PubChem CID, or normalised name."""

    def __init__(self, phase_path=BROAD_PHASE, molecule_path=BROAD_MOLECULE):
        import pandas as pd

        phase = pd.read_parquet(phase_path)
        molecule = pd.read_parquet(molecule_path)
        phase = phase.drop_duplicates("pert_iname").copy()
        phase["_name"] = phase["pert_iname"].map(normalize_name)
        self.by_name = phase.drop_duplicates("_name").set_index("_name")
        self.phase = phase.set_index("pert_iname")
        molecule = molecule.dropna(subset=["InChIKey"])
        molecule = molecule.assign(
            _key=molecule["InChIKey"].str.upper(),
            _name=molecule["pert_iname"].map(normalize_name),
            _cid=pd.to_numeric(molecule["pubchem_cid"], errors="coerce").astype("Int64"),
        )
        self.by_key = molecule.drop_duplicates("_key").set_index("_key")
        self.by_cid = molecule.dropna(subset=["_cid"]).drop_duplicates("_cid").set_index("_cid")

    def lookup(self, inchikey, cid, name):
        info = {}
        pert = None
        match = None
        key = str(inchikey or "").upper()
        if key in self.by_key.index:
            pert, match = self.by_key.at[key, "pert_iname"], "inchikey"
        elif cid is not None and int(cid) in {int(v) for v in self.by_cid.index}:
            pert, match = self.by_cid.loc[int(cid), "pert_iname"], "pubchem_cid"
        else:
            normalized = normalize_name(name)
            if normalized in self.by_name.index:
                pert, match = self.by_name.at[normalized, "pert_iname"], "name"
        if pert is None:
            return {"broad_hub_member": False, "broad_match": None}
        row = self.phase.loc[pert] if pert in self.phase.index else None
        info = {
            "broad_hub_member": True,
            "broad_match": match,
            "broad_pert_iname": pert,
            "broad_phase": str(row["clinical_phase"]) if row is not None else "",
            "broad_moa": str(row["moa"]) if row is not None else "",
            "broad_target": str(row["target"]) if row is not None else "",
            "broad_indication": str(row["indication"]) if row is not None else "",
            "broad_disease_area": str(row["disease_area"]) if row is not None else "",
        }
        for field in ("broad_phase", "broad_moa", "broad_target", "broad_indication",
                      "broad_disease_area"):
            if info.get(field, "").strip().lower() in {"nan", "none"}:
                info[field] = ""
        return info


# --------------------------------------------------------------------------
# BindingDB (6.2 GB TSV: stream with awk, never load whole)
# --------------------------------------------------------------------------
def bindingdb_matches(keys, cache_dir, path=DEFAULT_BINDINGDB, force=False):
    """Filter BindingDB to the study's InChIKeys and return per-key pairs."""
    import pandas as pd

    path = Path(path)
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    out = cache_dir / "bindingdb_matches.tsv"
    keys_file = cache_dir / "bindingdb_keys.txt"
    keys_file.write_text("\n".join(sorted(keys)) + "\n", encoding="utf-8")
    stamp = cache_dir / "bindingdb_matches.stamp"
    stale = (
        force
        or not out.exists()
        or not stamp.exists()
        or float(stamp.read_text() or 0) < path.stat().st_mtime
    )
    if stale:
        if not path.exists():
            return {"available": False, "path": str(path), "records": {}}
        tmp = out.with_suffix(".tsv.tmp")
        with open(tmp, "w", encoding="utf-8") as handle:
            subprocess.run(
                [
                    "awk",
                    "-F",
                    "\t",
                    "NR==FNR{key[$1]=1;next} ($4 in key)",
                    str(keys_file),
                    str(path),
                ],
                stdout=handle,
                check=True,
            )
        tmp.replace(out)
        stamp.write_text(str(path.stat().st_mtime), encoding="utf-8")
    columns = [
        "Ligand InChI Key",
        "Target Name",
        "BindingDB Ligand Name",
        "Ki (nM)",
        "IC50 (nM)",
        "Kd (nM)",
        "EC50 (nM)",
        "Curation/DataSource",
    ]
    frame = pd.read_csv(
        out, sep="\t", header=None, names=columns, dtype=str, on_bad_lines="skip"
    )
    records = {}
    for key, group in frame.groupby("Ligand InChI Key"):
        targets = [str(t).strip() for t in group["Target Name"].dropna()]
        unique = []
        for target in targets:
            if target and target not in unique:
                unique.append(target)
        records[str(key).strip().upper()] = {
            "n_pairs": int(len(group)),
            "targets": unique[:15],
        }
    return {"available": True, "path": str(path), "records": records}


# --------------------------------------------------------------------------
# Identity
# --------------------------------------------------------------------------
def smiles_inchikey(smiles):
    """Recompute the InChIKey of the frozen SMILES (rdkit, stereo included)."""
    try:
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        return Chem.MolToInchiKey(mol)
    except Exception:  # pragma: no cover - rdkit missing is recorded as a gap
        return ""


def pubchem_atc_codes(strings):
    """ATC codes from PubChem ATC entries (``CODE`` or ``CODE - Name``)."""
    codes = []
    for raw in strings:
        token = str(raw).split(" - ")[0].strip().upper()
        if re.fullmatch(r"[A-Z][A-Z0-9]{2,6}", token) and token not in codes:
            codes.append(token)
    return codes


def identity_record(row, pubchem, broad, bindingdb):
    """Check the frozen key/name against the sources; never merge structures."""
    key = row["drug_key"].strip()
    inchikey = str(row["inchikey"] or "").strip().upper()
    smiles_key = smiles_inchikey(row["smiles"]) if row.get("smiles") else ""
    cids = pubchem["cids"]
    resolved = bool(cids) or bool(broad.get("broad_hub_member")) or bool(
        bindingdb.get("n_pairs")
    )
    notes = []
    if not cids:
        notes.append("no pubchem cid")
    if inchikey and key != inchikey:
        notes.append("registered structure is a salt/parent of drug_key")
    if smiles_key and inchikey and smiles_key != inchikey:
        notes.append("frozen smiles InChIKey differs from frozen inchikey")
    if smiles_key and not inchikey:
        notes.append("frozen inchikey missing")
    name_matches = None
    synonyms = pubchem.get("synonyms") or []
    if synonyms:
        wanted = normalize_name(row["compound_name"])
        name_matches = wanted in {normalize_name(item) for item in synonyms}
        if not name_matches:
            notes.append("frozen name not among pubchem synonyms")
    elif cids:
        notes.append("pubchem synonyms unavailable")
    return {
        "identity_resolved": resolved,
        "pubchem_cids": cids,
        "smiles_key": smiles_key,
        "smiles_key_matches_frozen": bool(smiles_key and smiles_key == inchikey),
        "name_in_pubchem_synonyms": name_matches,
        "identity_notes": "; ".join(notes),
    }


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------
def frozen_keys(path):
    """Read the frozen compounds.csv and collapse 467 name records to 460 keys."""
    import pandas as pd

    frame = pd.read_csv(path, dtype=str)
    temp_env = frame[frame["temp_env_all-smiles"].str.lower() == "true"].copy()
    keys = []
    for key, group in temp_env.groupby("drug_key", sort=False):
        first = group.iloc[0]
        keys.append(
            {
                "drug_key": key,
                "all_names": "; ".join(dict.fromkeys(group["compound_name"].astype(str))),
                "n_name_records": int(len(group)),
                "compound_name": first["compound_name"],
                "inchikey": first["inchikey"],
                "smiles": first["smiles"],
                "atc_frozen": "" if pd.isna(first["atc"]) else str(first["atc"]),
                "moa_frozen": "" if pd.isna(first["moa"]) else str(first["moa"]),
                "mesh_frozen": ""
                if pd.isna(first["mesh_classes"])
                else str(first["mesh_classes"]),
                "old_tier_atc_moa": str(first["pharm_atc_moa"]).lower() == "true",
                "old_tier_atc_moa_mesh": str(first["pharm_mesh"]).lower() == "true",
                "old_tier_all_smiles": str(first["pharm_all_smiles"]).lower() == "true",
                "annotation_source": first.get("annotation_source", ""),
            }
        )
    return keys


def prefetch(rows, cache_dir, allow_network=True, jobs=3):
    """Warm the response cache in parallel; the builder then runs offline-fast."""
    if not allow_network or jobs <= 1:
        return 0
    cache_dir = Path(cache_dir)
    fetched = []

    def warm(row):
        pubchem = pubchem_cids(row["inchikey"], cache_dir, allow_network)
        if pubchem["cids"]:
            record = pubchem_record(pubchem["cids"][0], cache_dir, allow_network)
            return row["drug_key"], bool(record and record.get("body"))
        return row["drug_key"], False

    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        for key, ok in pool.map(warm, rows):
            if ok:
                fetched.append(key)
    return len(fetched)


def build_evidence(rows, cache_dir, allow_network=True, bindingdb_path=DEFAULT_BINDINGDB):
    """Return ``{drug_key: evidence dict}`` for every frozen key."""
    cache_dir = Path(cache_dir)
    broad = BroadHub()
    keys = [row["drug_key"] for row in rows]
    binding = bindingdb_matches(keys, cache_dir, path=bindingdb_path)
    binding_records = binding.get("records", {})

    evidence = {}
    for row in rows:
        key = row["drug_key"]
        pubchem = pubchem_cids(row["inchikey"], cache_dir, allow_network)
        cid = pubchem["cids"][0] if pubchem["cids"] else None
        extracted, record_url, record_ts, record_status = {}, None, None, None
        if cid is not None:
            cache = JsonlCache(cache_dir / "pubchem_pugview.jsonl")
            record = cache.get(
                f"{PUGVIEW}/data/compound/{cid}/JSON",
                allow_network=allow_network,
            )
            if record and record.get("body"):
                extracted = extract_pubchem(record["body"])
                record_url, record_ts = record["url"], record["fetched_at"]
                record_status = record.get("http_status")
        broad_info = broad.lookup(row["inchikey"], cid, row["compound_name"])
        binding_info = binding_records.get(str(row["inchikey"]).upper(), {})
        identity = identity_record(row, {**pubchem, **extracted}, broad_info, binding_info)

        # PubChem lists every ATC hierarchy level as its own entry; only
        # substance-level (7-character) human codes are counted, so the numbers
        # are not inflated by level prefixes.
        atc_pubchem_codes = [
            code for code in pubchem_atc_codes(extracted.get("atc_pubchem", []))
            if len(code) == 7 and not code.startswith("Q")
        ]
        atc_human = list(atc_pubchem_codes)
        urls = []
        if pubchem.get("url"):
            urls.append({"source": "pubchem-cids", "url": pubchem["url"],
                         "retrieved_at": pubchem["fetched_at"],
                         "http_status": pubchem.get("http_status")})
        if record_url:
            urls.append({"source": "pubchem-pugview", "url": record_url,
                         "retrieved_at": record_ts, "http_status": record_status})
        if broad_info.get("broad_hub_member"):
            urls.append(
                {
                    "source": "broad-repurposing-hub",
                    "url": "https://repo-hub.broadinstitute.org/repurposing",
                    "retrieved_at": file_timestamp(BROAD_PHASE) or now_iso(),
                    "http_status": "local",
                }
            )
        if binding_info.get("n_pairs"):
            binding_path = Path(binding.get("path") or DEFAULT_BINDINGDB)
            urls.append(
                {
                    "source": "bindingdb",
                    "url": "https://www.bindingdb.org/",
                    "retrieved_at": file_timestamp(binding_path) or now_iso(),
                    "http_status": "local",
                }
            )
        retrieved = max((u["retrieved_at"] or "" for u in urls), default="")
        evidence[key] = {
            **extracted,
            **broad_info,
            "drug_key": key,
            "inchikey": row["inchikey"],
            "pubchem_cids": pubchem["cids"],
            "pubchem_cid": cid,
            "atc_pubchem_codes": atc_pubchem_codes,
            "atc_human_codes": sorted(set(atc_human)),
            "bindingdb_targets": binding_info.get("targets", []),
            "bindingdb_n_pairs": int(binding_info.get("n_pairs", 0)),
            "evidence_urls": urls,
            "evidence_ok_urls": sum(
                1 for u in urls if u.get("http_status") in (200, "local")
            ),
            "retrieved_at": retrieved,
            **identity,
        }
    return evidence, {"bindingdb": binding}
