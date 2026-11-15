"""Drug-cardiotoxicity class sources, the pinned crosswalk and the class DEGs.

ADR 0066 makes a **mechanism class** the unit of drug cardiotoxicity evidence.
This module is the single live implementation of that science for ticket #189:

* **The pinned crosswalk.** Every drug code the three owned sources carry is
  resolved to a drug name (from the source's own pinned metadata) and to a
  mechanism class from ``pert_data/Annotation/local/cm_compound_annotation.csv``
  (``atc``/``moa``/``mesh_classes``). No class is hand-entered. Codes whose
  annotation row is missing or carries no ATC subgroup are *unresolved* and are
  excluded and recorded; a code with both a dataset-specific structure
  (``dataset_inchikey``) and a name-only row is keyed on the dataset-specific
  row, so a source's own analogue is never relabelled as the named drug.
* **The class rule.** A code's class is its ATC level-4 pharmacological
  subgroup (the five-character ATC prefix). A multi-ATC code is assigned only
  when exactly one candidate class is corroborated by the drug's
  ``mesh_classes``; otherwise it is excluded as ``ambiguous_multi_class``. The
  partition is mandatory: one code, one class, so the same well never enters two
  class DEGs.
* **TOP2 admission.** The topoisomerase-II-poison class is the ``L01DB``
  subgroup, and it enters a source only when every member's ``mesh_classes``
  contains ``Topoisomerase II Inhibitors``; a member missing that evidence is
  excluded from the class and recorded.
* **Eligibility.** A class enters a source only with **at least three
  small-molecule drugs** in that source. A small molecule is a code that
  resolved to a structure-backed annotation row; antibodies and peptides cannot
  satisfy that test and are excluded (their codes are recorded as unresolved).
  ``cardiotox_other`` is the DToxS source's non-TOP2 classes that pass the rule.
* **The class DEG.** The class's drug wells are pooled into one Disease group
  and compared with the same-plate vehicle in each source stratum (plate / cell
  line / time point), using the source's own test: the DToxS deposit's pinned
  edgeR pipeline on its read counts, DESeq2 on GSE243674's counts, and the
  declared Welch test on GSE244740's ScreenSeq counts. The pooled per-stratum
  contrasts then collapse per gene by the maximum effect through
  :mod:`bases.deg_rank` (the ticket-01 contract); a class is therefore never
  treated as several independent experiments.

Products (candidate-then-replace, ADR 0062):

* ``output/disease-deg/cardiotox/`` — the pinned crosswalk, class selection,
  membership and excluded-code evidence, plus its manifest.
* ``output/disease-deg/cardiotox_TOP2/`` — one DEG table + within-dataset rank
  per source, ``catalog.json`` and ``manifest.json``.
* ``output/disease-deg/cardiotox_other/`` — one DEG table + within-dataset rank
  per entered DToxS class, ``catalog.json`` and ``manifest.json``.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import re
import subprocess
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

from . import deg_rank, disease_deg, r_env
from .release import (ReleaseError, candidate, check_manifest, output_digests,
                      read_json, replace, sha256, write_json)

SCHEMA = "stage2-cardiotox-class-v1"
CATALOG_SCHEMA = "stage2-disease-deg-catalog-v1"
AXIS_SIZE = 18164
MIN_CLASS_DRUGS = 3
TOP2_CLASS_KEY = "L01DB"
TOP2_MECHANISM = "Topoisomerase II Inhibitors"
CATEGORY_TOP2 = "cardiotox_TOP2"
CATEGORY_OTHER = "cardiotox_other"
EVIDENCE_DIR = "cardiotox"
OUTPUT_ROOT = "stage2/01_bases/output/disease-deg"

AXIS = "stage2/01_bases/output/gene-identity/tables/genes.csv"
CROSSWALK = "stage2/01_bases/output/gene-identity/tables/gene_crosswalk.csv"
HGNC = "pert_data/Annotation/source/hgnc_complete_set.txt"
ANNOTATION = "pert_data/Annotation/local/cm_compound_annotation.csv"

DTOXS_ITEM = "pert_data/DToxS-iPSC"
DTOXS_SERIES = DTOXS_ITEM + "/source/GSE217421_series_matrix.txt.gz"
DTOXS_COUNTS = DTOXS_ITEM + "/source/GSE217421_Conv-RNAseq-Read-Counts.tsv.gz"
DTOXS_DESIGN = DTOXS_ITEM + "/source/GSE217421_Conv-RNAseq-Experiment-Design.tsv.gz"
DTOXS_PARAMS = DTOXS_ITEM + "/source/GSE217421_Conv-RNAseq-Params.tar.gz"
DTOXS_DGE_PARAMS = DTOXS_ITEM + "/source/GSE217423_DGE-RNAseq-Params.tar.gz"
DTOXS_DGE_SERIES = DTOXS_ITEM + "/source/GSE217423_series_matrix.txt.gz"
DTOXS_DEG_TAR = DTOXS_ITEM + "/source/GSE217421_Conv-DEG-Calc-FDR-0.1.tar.gz"
DTOXS_PROGRAMS = DTOXS_ITEM + "/source/GSE217421_Program-Source-Codes.tar.gz"

GSE243674_ITEM = "disease_data/GSE243674"
GSE243674_COUNTS = GSE243674_ITEM + "/source/GSE243674_Counts_RNA_ERMatthews.txt.gz"
GSE243674_DEG_JSON = GSE243674_ITEM + "/local/deg/GSE243674_deg.json"
GSE243674_DEG_CSV = GSE243674_ITEM + "/local/deg/GSE243674_deg.csv"

GSE244740_ITEM = "pert_data/GSE244740"
GSE244740_COUNTS = GSE244740_ITEM + "/source/GSE244740_processed_data_counts.txt.gz"
GSE244740_SERIES = GSE244740_ITEM + "/source/GSE244740_series_matrix.txt.gz"

TOOLS = "stage2/01_bases/execution/cardiotox_tools.py"

# The R build environment is configuration, not science: `bases.r_env` resolves
# `Rscript` from an explicit argument, `ARC_RSCRIPT`, or `PATH` in that order.
# The documented build is the `biomni_e1` conda env (R 4.5.3, edgeR 4.8.2,
# DESeq2 1.50.2); the `arc` env's R has neither package. Manifests record the
# environment name/version only, never a host path.

EDGE_R_SCRIPT = "stage2/01_bases/output/api/bases/cardiotox_edgeR.R"
DESEQ2_SCRIPT = "stage2/01_bases/output/api/bases/cardiotox_deseq2.R"

SOURCES = {
    "gse243674": {
        "item": GSE243674_ITEM,
        "accession": "GSE243674",
        "species": "Homo sapiens",
        "model_system": "hiPSC-CM monolayer (6 healthy female donors, 0.5 uM, 3 h / 24 h)",
        "shared_control_group": "GSE243674:same-time vehicle wells (6 donors per time point, shared by every drug arm)",
        "vehicle_design": "12 vehicle wells (6 donors x 2 time points); each drug contrast is DESeq2 ~ donor + group against the same-time vehicle wells",
    },
    "dtoxs_ipsc": {
        "item": DTOXS_ITEM,
        "accession": "GSE217421",
        "species": "Homo sapiens",
        "model_system": "DToxS hiPSC-CM (GSE217424 SuperSeries; 48 h, 6 cell lines)",
        "shared_control_group": "DToxS-iPSC:plate DMSO wells (the shipped per-drug Calc tables share one same-plate CTRL set per cell line)",
        "vehicle_design": "each shipped per-cell-line Calc table compares a drug's wells with the same cell line's DMSO CTRL wells; the pooled class contrast reuses that same CTRL set",
    },
    "gse244740": {
        "item": GSE244740_ITEM,
        "accession": "GSE244740",
        "species": "Homo sapiens",
        "model_system": "ScreenSeq hiPSC-CM panel (24 h, plate S4881Nr8)",
        "shared_control_group": "GSE244740:same-plate DMSO wells at 24 h (shared by the three anthracycline/anthraquinone arms)",
        "vehicle_design": "same-plate DMSO_TF wells at 24 h on plate S4881Nr8",
    },
}

WHITESPACE = re.compile(r"\s+")


class CardiotoxError(ReleaseError):
    """A cardiotoxicity class product violates its contract."""


# --------------------------------------------------------------------------- #
# small file helpers
# --------------------------------------------------------------------------- #
def _open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open(newline="")


def _series_characteristics(path):
    """Ordered ``!Sample_characteristics_ch1`` blocks of a GEO series matrix."""
    blocks = []
    with _open_text(path) as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            fields = [value.strip().strip('"') for value in line.rstrip("\n").split("\t")]
            if fields[0] == "!Sample_characteristics_ch1":
                blocks.append(fields[1:])
    return blocks


def dtoxs_series_names(path):
    """``treatment condition`` code -> ``drug name`` from the Conv series matrix."""
    conditions, names = None, None
    for block in _series_characteristics(path):
        if not block:
            continue
        key = block[0].split(":", 1)[0].strip()
        if key == "treatment condition":
            conditions = [value.split(":", 1)[1].strip() for value in block]
        elif key == "drug name":
            names = [value.split(":", 1)[1].strip() for value in block]
    if not conditions or not names or len(conditions) != len(names):
        raise CardiotoxError("the DToxS series matrix lacks paired treatment/drug-name rows")
    mapping = {}
    for code, name in zip(conditions, names):
        mapping.setdefault(code, set()).add(name)
    conflicts = {code: sorted(values) for code, values in mapping.items() if len(values) > 1}
    if conflicts:
        raise CardiotoxError(f"a DToxS code maps to several drug names: {conflicts}")
    return {code: next(iter(values)) for code, values in mapping.items() if code != "CTRL"}


def dtoxs_drug_groups(tar_path, member="Conv-RNAseq-Drug-Groups.tsv"):
    """The pinned code roster of one DToxS platform (``Group``/``Drug``/``Type``)."""
    with tarfile.open(tar_path, "r|gz") as archive:
        for entry in archive:
            if Path(entry.name).name != member:
                continue
            frame = pd.read_csv(io.BytesIO(archive.extractfile(entry).read()),
                                sep="\t", dtype=str, keep_default_na=False)
            if list(frame.columns[:3]) != ["Group", "Drug", "Type"]:
                raise CardiotoxError(f"{member} lacks Group/Drug/Type columns")
            return frame
    raise CardiotoxError(f"{member} not found in {Path(tar_path).name}")


def gse244740_panel(series_path):
    """Treated compound names and well metadata from the GSE244740 series matrix."""
    blocks = _series_characteristics(series_path)
    descriptions, treatments, times = None, None, None
    with _open_text(series_path) as handle:
        for line in handle:
            if line.startswith("!Sample_description"):
                descriptions = [value.strip().strip('"')
                                for value in line.rstrip("\n").split("\t")][1:]
    for block in blocks:
        if not block:
            continue
        key = block[0].split(":", 1)[0].strip()
        if key == "treatment":
            treatments = [value.split(":", 1)[1].strip() for value in block]
        elif key == "time":
            times = [value.split(":", 1)[1].strip() for value in block]
    if not descriptions or not treatments or not times:
        raise CardiotoxError("the GSE244740 series matrix lacks well/treatment/time rows")
    if not (len(descriptions) == len(treatments) == len(times)):
        raise CardiotoxError("the GSE244740 series matrix rows are ragged")
    records = []
    for well, treatment, time in zip(descriptions, treatments, times):
        compound, _, dose = treatment.partition(",")
        records.append({"well": well, "compound": compound.strip(),
                        "dose_index": dose.strip(), "time": time,
                        "plate": well.split("_")[0]})
    return pd.DataFrame(records)


def gse243674_contrast_metadata(json_path):
    """The item's pinned contrast metadata: token -> drug name and design."""
    value = read_json(json_path)
    contrasts = value.get("contrasts")
    if not contrasts:
        raise CardiotoxError("GSE243674 deg.json declares no contrasts")
    return contrasts


# --------------------------------------------------------------------------- #
# annotation / mechanism class
# --------------------------------------------------------------------------- #
def read_annotation(path):
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"lookup_name", "dataset_inchikey", "parent_inchikey", "atc", "moa",
                "mesh_classes"}
    if not required <= set(frame.columns):
        raise CardiotoxError("compound annotation lacks required columns")
    frame["_name_key"] = frame.lookup_name.map(_name_key)
    return frame


def _name_key(value):
    return WHITESPACE.sub(" ", str(value).strip()).casefold()


def _series_key(value):
    """Tolerant key for the deposited drug-name strings (spaces vs hyphens)."""
    return re.sub(r"[^a-z0-9]", "", str(value).casefold())


_ATC_LABELS = re.compile(r"(?:^|;)\s*([A-Z]\d{2}[A-Z]{2}(?:\d{2})?)\b\s*([^;]*)")


def atc_entries(field):
    """``[(code, label)]`` for every ATC code in an annotation ``atc`` field.

    A level-4-only code (``C09CA``) is accepted as-is; a seven-character code
    keeps its level-5 suffix and is reduced to its level-4 subgroup.
    """
    entries = []
    for code, label in _ATC_LABELS.findall(str(field)):
        entries.append((code, label.strip(" ,;")))
    return entries


def atc_level4(code):
    """The ATC level-4 pharmacological subgroup of a seven-character ATC code."""
    code = str(code).strip()
    if len(code) < 5 or not re.fullmatch(r"[A-Z]\d{2}[A-Z]{2}\d{0,2}", code):
        raise CardiotoxError(f"not an ATC code: {code}")
    return code[:5]


def mesh_terms(field):
    return [term.strip() for term in str(field).split(";") if term.strip()]


def _mechanism_corroborates(label, terms):
    """Does a class label share a mechanism word with the drug's mesh terms?

    Used only to disambiguate a *multi-ATC* code; a single-candidate code never
    needs it. The comparison is a documented word-stem overlap: a label word of
    at least five characters must appear inside a mesh term (or the reverse).
    """
    words = [word for word in re.findall(r"[a-z]{5,}", label.casefold())
             if word not in {"other", "related", "substances", "agents"}]
    if not words:
        return False
    lowered = [term.casefold() for term in terms]
    return any(word[:6] in term for word in words for term in lowered)


def annotate_code(name, annotation):
    """Resolve one deposited drug name to its pinned annotation row.

    Returns a dict with ``status`` in ``{'resolved', 'unresolved_no_annotation'}``
    and the comparable ``lookup_name``/structure evidence. The
    ``(lookup_name, dataset_inchikey)`` key is honoured: when the annotation
    carries a dataset-specific structure for the name, that row is the one used,
    so a source's analogue is never relabelled as the named drug.
    """
    name = str(name).strip()
    if not name:
        return {"drug_name": name, "status": "unresolved_no_annotation", "reason": "empty name"}
    rows = annotation[annotation._name_key.eq(_name_key(name))]
    if rows.empty:
        # Deposited spellings occasionally differ from the registry name by
        # punctuation only; a repository-maintained comment field is never used
        # to rewrite a compound.
        return {"drug_name": name, "status": "unresolved_no_annotation",
                "reason": "no annotation row for the deposited name"}
    dataset_rows = rows[rows.dataset_inchikey.str.strip().ne("")]
    pool = dataset_rows if not dataset_rows.empty else rows
    structures = set(pool.dataset_inchikey.str.strip()) - {""}
    structural = structures if structures else set(pool.parent_inchikey.str.strip())
    if len(structural) > 1:
        return {"drug_name": name, "status": "unresolved_ambiguous_structure",
                "reason": "the annotation carries several structures for the name",
                "candidates": sorted(structural)}
    # Deterministic pick among exact duplicates.
    ordered = pool.sort_values(["dataset_inchikey", "parent_inchikey", "atc", "moa",
                                "mesh_classes", "notes"], kind="stable")
    row = ordered.iloc[0]
    return {
        "drug_name": name,
        "status": "resolved",
        "lookup_name": str(row.lookup_name),
        "dataset_inchikey": str(row.dataset_inchikey).strip(),
        "parent_inchikey": str(row.parent_inchikey).strip(),
        "atc": str(row.atc),
        "moa": str(row.moa),
        "mesh_classes": str(row.mesh_classes),
        "notes": str(row.notes),
    }


def classify_row(row):
    """Assign the mechanism class of one resolved annotation row.

    Returns a dict with ``status`` in ``{'assigned', 'unresolved_no_class',
    'ambiguous_multi_class'}`` and the ATC evidence. A code cannot be split
    across classes: the assignment is either unique or refused.
    """
    entries = atc_entries(row.get("atc", ""))
    by_class = {}
    for code, label in entries:
        level4 = atc_level4(code)
        by_class.setdefault(level4, label or by_class.get(level4, ""))
    candidates = sorted(by_class)
    if not candidates:
        return {"status": "unresolved_no_class", "class_id": "",
                "class_label": "", "atc_level4": "",
                "reason": "the annotation row carries no ATC subgroup"}
    if len(candidates) == 1:
        level4 = candidates[0]
        return {"status": "assigned", "class_id": f"atc_l4:{level4}",
                "class_label": by_class[level4], "atc_level4": level4,
                "atc_candidates": [f"atc_l4:{level4}"], "reason": ""}
    terms = mesh_terms(row.get("mesh_classes", ""))
    corroborated = [level4 for level4 in candidates
                    if _mechanism_corroborates(by_class[level4], terms)]
    if len(corroborated) == 1:
        level4 = corroborated[0]
        return {"status": "assigned", "class_id": f"atc_l4:{level4}",
                "class_label": by_class[level4], "atc_level4": level4,
                "atc_candidates": [f"atc_l4:{code}" for code in candidates],
                "reason": ""}
    return {"status": "ambiguous_multi_class", "class_id": "",
            "class_label": "", "atc_level4": "",
            "atc_candidates": [f"atc_l4:{code}" for code in candidates],
            "reason": "several ATC subgroups and none uniquely corroborated by mesh_classes"}


# --------------------------------------------------------------------------- #
# the crosswalk
# --------------------------------------------------------------------------- #
CROSSWALK_COLUMNS = [
    "source", "source_id", "drug_name", "lookup_name", "dataset_inchikey",
    "parent_inchikey", "atc", "moa", "mesh_classes", "atc_level4", "class_id",
    "class_label", "mechanism_terms", "small_molecule", "status",
    "exclusion_reason", "notes",
]


def source_records(step):
    """Every ``(source, source_id, drug_name)`` the three owned sources carry."""
    step = Path(step)
    records = []
    names = dtoxs_series_names(step / DTOXS_SERIES)
    groups = dtoxs_drug_groups(step / DTOXS_PARAMS)
    for code in groups.Drug.astype(str):
        if code not in names:
            raise CardiotoxError(f"DToxS code {code} has no series-matrix drug name")
        records.append({"source": "dtoxs_ipsc", "source_id": code,
                        "drug_name": names[code], "kind": "dtox_code"})
    for entry in gse243674_contrast_metadata(step / GSE243674_DEG_JSON):
        records.append({"source": "gse243674", "source_id": str(entry["drug_token"]),
                        "drug_name": str(entry["drug"]), "kind": "dtox_code"})
    panel = gse244740_panel(step / GSE244740_SERIES)
    compounds = sorted(set(panel.compound) - {"DMSO", "NoTreatment"})
    for name in compounds:
        records.append({"source": "gse244740", "source_id": name,
                        "drug_name": name, "kind": "compound_name"})
    seen = set()
    unique = []
    for record in records:
        key = (record["source"], record["source_id"].casefold())
        if key in seen:
            continue
        seen.add(key)
        unique.append(record)
    return unique


def build_crosswalk(step, annotation):
    """Resolve every carried code to one drug name and one mechanism class."""
    rows = []
    for record in source_records(step):
        annotation_row = annotate_code(record["drug_name"], annotation)
        row = {column: "" for column in CROSSWALK_COLUMNS}
        row.update({"source": record["source"], "source_id": record["source_id"],
                    "drug_name": record["drug_name"]})
        if annotation_row["status"] != "resolved":
            row.update({"status": annotation_row["status"],
                        "exclusion_reason": annotation_row.get("reason", ""),
                        "small_molecule": False})
            rows.append(row)
            continue
        row.update({
            "lookup_name": annotation_row["lookup_name"],
            "dataset_inchikey": annotation_row["dataset_inchikey"],
            "parent_inchikey": annotation_row["parent_inchikey"],
            "atc": annotation_row["atc"],
            "moa": annotation_row["moa"],
            "mesh_classes": annotation_row["mesh_classes"],
            "mechanism_terms": "; ".join(mesh_terms(annotation_row["mesh_classes"])),
            "notes": annotation_row["notes"],
        })
        classification = classify_row(annotation_row)
        row.update({"status": classification["status"],
                    "class_id": classification["class_id"],
                    "class_label": classification["class_label"],
                    "atc_level4": classification["atc_level4"],
                    "exclusion_reason": classification["reason"],
                    # a resolved structure-backed annotation row is the
                    # small-molecule evidence; antibodies/peptides have none.
                    "small_molecule": True})
        rows.append(row)
    frame = pd.DataFrame(rows, columns=CROSSWALK_COLUMNS)
    return frame.sort_values(["source", "source_id"], kind="stable").reset_index(drop=True)


def _top2_corroborated(row):
    return TOP2_MECHANISM in mesh_terms(row.mesh_classes)


def select_classes(crosswalk):
    """Apply the class rule: partition, TOP2 admission, >=3 small molecules.

    Returns ``(classes, members, exclusions)``:

    * ``classes`` — one row per ``(source, class)`` with the counts, the decision
      and the roster category it enters (empty when it does not enter);
    * ``members`` — the member drugs of every entered class;
    * ``exclusions`` — every code excluded from its class/entry with the reason.
    """
    members = crosswalk[crosswalk.small_molecule.astype(bool)
                        & crosswalk.class_id.ne("")]
    records, dropped = [], []
    # Every code the crosswalk could not assign to a class is excluded up front,
    # with the resolution/classification reason recorded verbatim.
    for row in crosswalk[crosswalk.status.ne("assigned")].itertuples(index=False):
        dropped.append({"source": row.source, "source_id": row.source_id,
                        "drug_name": row.drug_name, "class_id": row.class_id,
                        "reason": row.exclusion_reason or row.status})
    for (source, class_id), part in members.groupby(["source", "class_id"], sort=True):
        label = str(part.class_label.iloc[0])
        mechanism = sorted({term for terms in part.mechanism_terms
                            for term in mesh_terms(terms)})
        top2 = class_id == f"atc_l4:{TOP2_CLASS_KEY}"
        if top2:
            corroborated = part[part.apply(_top2_corroborated, axis=1)]
            missing = part[~part.index.isin(corroborated.index)]
            for row in missing.itertuples(index=False):
                dropped.append({"source": source, "source_id": row.source_id,
                                "drug_name": row.drug_name, "class_id": class_id,
                                "reason": "top2_mesh_missing"})
            included = corroborated
        else:
            included = part
        n = int(len(included))
        across = int(crosswalk[(crosswalk.small_molecule.astype(bool))
                               & crosswalk.class_id.eq(class_id)]
                     .source_id.nunique())
        if top2:
            category = CATEGORY_TOP2 if n >= MIN_CLASS_DRUGS else ""
            reason = "" if category else "class_below_three"
        elif source == "dtoxs_ipsc":
            category = CATEGORY_OTHER if n >= MIN_CLASS_DRUGS else ""
            reason = "" if category else "class_below_three"
        else:
            category = ""
            reason = "not_a_roster_category"
        records.append({
            "source": source,
            "class_id": class_id,
            "class_label": label,
            "atc_level4": class_id.split(":", 1)[1] if ":" in class_id else "",
            "mechanism_terms": "; ".join(mechanism),
            "n_small_molecule_drugs": n,
            "n_small_molecule_across_sources": across,
            "n_dropped_mesh_missing": int(len(part) - n) if top2 else 0,
            "members": "; ".join(sorted(included.drug_name)),
            "member_source_ids": "; ".join(sorted(included.source_id)),
            "category": category,
            "decision": "entered" if category else "not_entered",
            "reason": reason,
        })
        if not category:
            for row in included.itertuples(index=False):
                dropped.append({"source": source, "source_id": row.source_id,
                                "drug_name": row.drug_name, "class_id": class_id,
                                "reason": reason})
    classes = pd.DataFrame(records).sort_values(["source", "class_id"],
                                                kind="stable").reset_index(drop=True)
    entered = classes[classes.decision.eq("entered")]
    member_rows = []
    for part in members.merge(entered[["source", "class_id", "category"]],
                              on=["source", "class_id"], how="inner").itertuples(index=False):
        if part.class_id == f"atc_l4:{TOP2_CLASS_KEY}" and not _top2_corroborated(part):
            continue
        member_rows.append({"source": part.source, "class_id": part.class_id,
                            "category": part.category, "source_id": part.source_id,
                            "drug_name": part.drug_name,
                            "lookup_name": part.lookup_name,
                            "small_molecule": bool(part.small_molecule)})
    member_frame = pd.DataFrame(member_rows).sort_values(
        ["source", "class_id", "source_id"], kind="stable").reset_index(drop=True)
    exclusion_frame = pd.DataFrame(dropped).drop_duplicates(
        subset=["source", "source_id", "class_id", "reason"]).sort_values(
        ["source", "source_id", "class_id", "reason"], kind="stable").reset_index(drop=True)
    return classes, member_frame, exclusion_frame


# --------------------------------------------------------------------------- #
# contrast construction: pooled class wells vs the same-plate vehicle
# --------------------------------------------------------------------------- #
def dtox_pooled_strata(tar_path, member_codes):
    """Per-stratum (cell line, plate) control wells and pooled member wells.

    The shipped Calc tables record exactly the wells the deposit's own pipeline
    used after its sample filtering, so pooling their member wells inherits the
    source's QC. Every member table in a stratum must list the same CTRL wells;
    a disagreement is an error rather than a silent choice.
    """
    junk = {"logFC", "logCPM", "PValue", "FDR", "Regulation"}
    pattern = re.compile(r"Human-(.+)-CM-Hour-(\d+)-Plate-(\d+)-Calc\.CTRL-([A-Z0-9]+)\.tsv$")
    strata = {}
    with tarfile.open(tar_path, "r|gz") as archive:
        for entry in archive:
            match = pattern.search(Path(entry.name).name)
            if not match:
                continue
            cell, hour, plate, code = match.groups()
            if code not in member_codes:
                continue
            header = archive.extractfile(entry).readline().decode().rstrip("\n").split("\t")
            columns = [value for value in header
                       if not value.endswith(".Norm") and value not in junk]
            controls = tuple(value for value in columns if value.startswith("CTRL."))
            treated = tuple(value for value in columns if not value.startswith("CTRL."))
            stratum = (cell, plate, hour)
            record = strata.setdefault(stratum, {"controls": controls, "treated": {}})
            if record["controls"] != controls:
                raise CardiotoxError(f"stratum {stratum} member tables use different CTRL wells")
            record["treated"][code] = treated
    if not strata:
        raise CardiotoxError("no DToxS stratum carries the class members")
    return {stratum: {"controls": list(record["controls"]),
                      "treated": {code: list(wells)
                                  for code, wells in sorted(record["treated"].items())}}
            for stratum, record in sorted(strata.items())}


def pooled_contrast_rows(frame, contrasts):
    """Collapse per-stratum pooled contrasts; ``frame`` needs contract columns.

    ``contrasts`` is the list of ``{"contrast", "stratum", "n_case", "n_control"}``
    kept for the catalog; the collapse itself is :mod:`bases.deg_rank`.
    """
    return frame, contrasts


def write_contrast_inputs(root, name, frame, samples):
    """Write the counts/samples tables an R stratum script consumes."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    counts_path = root / f"{name}.counts.tsv"
    samples_path = root / f"{name}.samples.tsv"
    frame.to_csv(counts_path, sep="\t", index=False, lineterminator="\n")
    samples.to_csv(samples_path, sep="\t", index=False, lineterminator="\n")
    return counts_path, samples_path


def run_r(script, arguments, record_path=None, rscript=None):
    """Run one pinned R compute step under the resolved R build environment."""
    command = [r_env.rscript(rscript), str(script)] + [str(argument) for argument in arguments]
    environment = {"OMP_NUM_THREADS": "2", "MKL_NUM_THREADS": "2",
                   "OPENBLAS_NUM_THREADS": "2", "R_DEFAULT_PACKAGES": "datasets,utils,grDevices,graphics,stats,methods"}
    timed = ["/usr/bin/time", "-v", "-o", str(record_path)] if record_path else []
    completed = subprocess.run(timed + command, capture_output=True, text=True,
                               env={**environment, "PATH": "/usr/bin:/bin:/usr/local/bin"},
                               check=False)
    if completed.returncode != 0:
        raise CardiotoxError(
            f"R step failed ({Path(script).name}): {completed.stderr.strip()[-2000:]}")
    return completed


def r_environment(rscript=None):
    """R and package versions that produced the class DEGs."""
    executable = r_env.rscript(rscript)
    command = [executable, "--vanilla", "-e",
               ('cat(R.version.string, "\\n"); for (p in c("edgeR", "DESeq2")) '
                'cat(p, as.character(packageVersion(p)), "\\n")')]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise CardiotoxError(f"cannot read the R environment: {completed.stderr[-800:]}")
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    # The manifest may not carry an absolute path: name the conda env that
    # supplies the R build instead of the host path.
    environment_name = Path(executable).resolve().parent.parent.name
    result = {"interpreter": f"conda:{environment_name}",
              "r": lines[0] if lines else ""}
    for line in lines[1:]:
        name, _, version = line.partition(" ")
        result[name] = version
    return result


def _read_long(path):
    frame = pd.read_csv(path, sep="\t")
    required = {"gene_symbol", "log2fc", "padj"}
    if not required <= set(frame.columns):
        raise CardiotoxError(f"{Path(path).name} lacks {sorted(required)}")
    return frame


def _format(value):
    """Deterministic text for a numeric table cell; blank when missing."""
    if value is None:
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return ""
    return f"{number:.10g}"


# --------------------------------------------------------------------------- #
# published tables
# --------------------------------------------------------------------------- #
# The v2 published compact-table contract (#208, ADR 0071): the authoritative
# pass-any-contrast ``deg`` flag, the DEG-only rank and ``rank_norm = rank/n_DEG``.
# The retired measured-gene rank/denominator is gone.
DEG_COLUMNS = ["gene_index", "gene_symbol", "contrast", "log2fc", "stat",
               "p_value", "padj", "base_mean", "direction", "significant",
               "deg", "n_contrasts", "n_contrasts_supporting",
               "n_contrasts_deg", "effect", "rank", "rank_norm"]
RANK_COLUMNS = ["gene_index", "rank", "rank_norm"]


def class_deg(rows, resolver, symbol_column="gene_symbol", axis_size=AXIS_SIZE):
    """Align, collapse and rank one ``(source, class)`` DEG table.

    ``rows`` carries one row per per-stratum contrast with ``gene_symbol``,
    ``log2fc``, ``padj`` and optionally ``stat``/``p_value``/``base_mean``.
    The collapse and rank come from :mod:`bases.deg_rank` (ticket-01 contract);
    the winning contrast's statistic/p-value/expression level are carried onto
    the published row.
    """
    aligned, report = deg_rank.align_table(resolver, rows, symbol_column=symbol_column)
    if aligned.empty:
        raise CardiotoxError("a class DEG has no axis-resolved measurement")
    for column in ("stat", "p_value", "base_mean"):
        if column not in aligned:
            aligned[column] = np.nan
    effect = deg_rank.deg_effect(aligned.log2fc.to_numpy(dtype=np.float64),
                                 aligned.padj.to_numpy(dtype=np.float64))
    ordered = aligned.assign(effect=effect).sort_values(
        ["gene_index", "effect", "contrast"], ascending=[True, False, True],
        kind="stable")
    table, rank, stats = deg_rank.dataset_deg_rank(
        aligned[["gene_index", "log2fc", "padj", "contrast"]], axis_size=axis_size)
    winning = ordered.drop_duplicates(["gene_index", "contrast"], keep="first")
    winning = winning.drop_duplicates("gene_index", keep="first").set_index("gene_index")
    table = table.assign(
        stat=table.gene_index.map(winning["stat"]),
        p_value=table.gene_index.map(winning["p_value"]),
        base_mean=table.gene_index.map(winning["base_mean"]),
    )
    table = table.sort_values(["rank", "gene_index"], kind="stable").reset_index(drop=True)
    stats["alignment"] = report
    return table, rank, stats


def write_deg_table(path, table, axis):
    """Write the compact published DEG table (contract columns + symbol)."""
    symbols = axis.set_index("gene_index").symbol
    frame = pd.DataFrame({
        "gene_index": table.gene_index.astype(np.int64),
        "gene_symbol": table.gene_index.map(symbols),
        "contrast": table.contrast.astype(str),
        "log2fc": table.log2fc.map(_format),
        "stat": table.stat.map(_format),
        "p_value": table.p_value.map(_format),
        "padj": table.padj.map(_format),
        "base_mean": table.base_mean.map(_format),
        "direction": table.direction.astype(np.int64),
        "significant": table.significant.astype(bool),
        "deg": table.deg.astype(bool),
        "n_contrasts": table.n_contrasts.astype(np.int64),
        "n_contrasts_supporting": table.n_contrasts_supporting.astype(np.int64),
        "n_contrasts_deg": table.n_contrasts_deg.astype(np.int64),
        "effect": table.effect.map(_format),
        "rank": table["rank"].astype(np.int64),
        "rank_norm": table.rank_norm.map(_format),
    })
    path = Path(path)
    frame.to_csv(path, index=False, lineterminator="\n")
    return frame


def write_rank_table(path, rank):
    frame = pd.DataFrame({
        "gene_index": rank.gene_index.astype(np.int64),
        "rank": rank["rank"].astype(np.int64),
        "rank_norm": rank.rank_norm.map(_format),
    })
    path = Path(path)
    frame.to_csv(path, index=False, lineterminator="\n")
    return frame


def _canonical_catalog(path):
    """Return the catalog without its digest field and the claimed digest."""
    catalog = read_json(path)
    claimed = catalog.pop("catalog_sha256", None)
    if not isinstance(claimed, str) or len(claimed) != 64:
        raise CardiotoxError("catalog.json lacks a valid catalog_sha256")
    return catalog, claimed


def catalog_digest(path):
    """Verify a category catalog's payload digest; returns the verified digest.

    ``catalog_sha256`` is the sha256 of the canonical JSON of the catalog with
    that field removed (``write_json``'s sorted, two-space-indented form plus a
    trailing newline) — the only self-consistent meaning for a file digest
    stored inside the file it seals.
    """
    catalog, claimed = _canonical_catalog(path)
    payload = (json.dumps(catalog, indent=2, sort_keys=True) + "\n").encode()
    computed = hashlib.sha256(payload).hexdigest()
    if computed != claimed:
        raise CardiotoxError(
            f"catalog payload digest differs for {Path(path).name}")
    return computed


def _rows(path):
    with Path(path).open(newline="") as handle:
        return sum(1 for _ in handle) - 1


def _deg_columns():
    return list(DEG_COLUMNS)


# --------------------------------------------------------------------------- #
# public readers
# --------------------------------------------------------------------------- #
class DiseaseDeg:
    """Read-only view of the published cardiotoxicity class evidence."""

    def __init__(self, root):
        self.root = Path(root)
        self.manifest = read_json(self.root / "manifest.json")
        if self.manifest.get("schema") != SCHEMA:
            raise CardiotoxError("unsupported cardiotox schema")
        self.axis_size = int(self.manifest["axis_size"])

    def crosswalk(self):
        return pd.read_csv(self.root / "crosswalk.csv", keep_default_na=False, dtype=str)

    def classes(self):
        return pd.read_csv(self.root / "classes.csv", keep_default_na=False, dtype=str)

    def members(self):
        return pd.read_csv(self.root / "class_members.csv", keep_default_na=False, dtype=str)

    def excluded(self):
        return pd.read_csv(self.root / "excluded.csv", keep_default_na=False, dtype=str)

    def eligibility(self):
        return pd.read_csv(self.root / "eligibility.csv", keep_default_na=False, dtype=str)


class CategoryDeg:
    """One published category: its catalog and its per-dataset tables."""

    def __init__(self, root, category):
        self.root = Path(root)
        self.category = category
        self.catalog = read_json(self.root / "catalog.json")
        if self.catalog.get("schema") != CATALOG_SCHEMA:
            raise CardiotoxError("unsupported catalog schema")
        if self.catalog.get("category") != category:
            raise CardiotoxError(f"catalog category differs from {category}")
        self.catalog_digest = catalog_digest(self.root / "catalog.json")

    def dataset_ids(self):
        return [entry["dataset_id"] for entry in self.catalog["datasets"]]

    def dataset(self, dataset_id):
        for entry in self.catalog["datasets"]:
            if entry["dataset_id"] == dataset_id:
                return entry
        raise CardiotoxError(f"unknown dataset: {dataset_id}")

    def table(self, dataset_id):
        entry = self.dataset(dataset_id)
        return pd.read_csv(self.root / entry["deg"], keep_default_na=False)

    def rank(self, dataset_id):
        entry = self.dataset(dataset_id)
        return pd.read_csv(self.root / entry["rank"], keep_default_na=False)


# --------------------------------------------------------------------------- #
# pinned inputs and the frozen axis
# --------------------------------------------------------------------------- #
def add_input(inputs, path, repo):
    """Record one repo-relative input path and its digest (never an absolute path)."""
    path = Path(path)
    key = path.resolve().relative_to(Path(repo).resolve()).as_posix()
    if key not in inputs:
        inputs[key] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    return key


def load_axis(step):
    frame = pd.read_csv(Path(step) / AXIS, dtype=str, keep_default_na=False)
    frame = frame.assign(gene_index=np.arange(len(frame), dtype=np.int64))
    if len(frame) != AXIS_SIZE:
        raise CardiotoxError(f"the frozen axis has {len(frame)} rows, expected {AXIS_SIZE}")
    return frame


def load_resolver(step):
    step = Path(step)
    axis = load_axis(step)
    crosswalk_rows = pd.read_csv(step / CROSSWALK, dtype=str,
                                 keep_default_na=False).to_dict("records")
    hgnc = pd.read_csv(step / HGNC, sep="\t", dtype=str, keep_default_na=False,
                       usecols=["hgnc_id", "alias_symbol", "prev_symbol",
                                "ensembl_gene_id"]).to_dict("records")
    return axis, deg_rank.SymbolResolver(axis, crosswalk_rows, hgnc)


def python_environment():
    import sys

    import numpy
    import pandas
    try:
        from rdkit import rdBase
        rdkit_version = rdBase.rdkitVersion
    except Exception:  # pragma: no cover - rdkit is an optional pin here
        rdkit_version = ""
    return {"python": sys.version.split()[0], "numpy": numpy.__version__,
            "pandas": pandas.__version__, "rdkit_version": rdkit_version}


# --------------------------------------------------------------------------- #
# source-specific pooled contrasts
# --------------------------------------------------------------------------- #
def _well_meta(wells):
    """Decode GSE243674 ``<DRUG>.<donor>.<time>h`` column labels."""
    records = []
    for well in wells:
        match = re.fullmatch(r"([A-Z]+)\.([0-9]+)\.([0-9]+)h", str(well))
        if not match:
            raise CardiotoxError(f"cannot parse GSE243674 well label: {well}")
        token, donor, time = match.groups()
        records.append({"well": str(well), "token": token,
                        "donor": int(donor), "time_h": int(time)})
    return pd.DataFrame(records)


def _read_r_output(path):
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    frame = frame.replace("", np.nan)
    for column in ("logFC", "logCPM", "log2FC", "stat", "PValue", "p_value",
                   "FDR", "padj", "baseMean"):
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def _r_peak(path):
    if path is None or not Path(path).is_file():
        return None
    text = Path(path).read_text(errors="replace")
    match = re.search(r"Maximum resident set size \(kbytes\):\s*([0-9]+)", text)
    return int(match.group(1)) if match else None


def gse243674_class_contrasts(step, member_tokens, class_slug, workdir):
    """Pooled class-vs-vehicle DESeq2 contrasts, one per time point."""
    step, workdir = Path(step), Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(step / GSE243674_COUNTS, sep=" ", index_col=0)
    counts.index = counts.index.astype(str)
    meta = _well_meta(counts.columns)
    members = set(member_tokens)
    unknown = members - set(meta.token)
    if unknown:
        raise CardiotoxError(f"unknown GSE243674 member tokens: {sorted(unknown)}")
    frames, contrasts, peaks = [], [], []
    for time_h in sorted(meta.time_h.unique()):
        treated = meta[(meta.time_h == time_h) & meta.token.isin(members)]
        control = meta[(meta.time_h == time_h) & meta.token.eq("VEH")]
        if treated.empty or control.empty:
            raise CardiotoxError(f"GSE243674 time {time_h} lacks an arm")
        order = list(control.well) + list(treated.well)
        samples = pd.DataFrame({
            "well": order,
            "donor": [int(meta.set_index("well").loc[well, "donor"]) for well in order],
            "group": ["vehicle"] * len(control) + ["CLASS"] * len(treated),
            "time_h": time_h,
        })
        matrix = counts[order].copy()
        matrix.insert(0, "ENTREZID", counts.index)
        counts_path = workdir / f"gse243674_{class_slug}_{time_h}h.counts.tsv"
        samples_path = workdir / f"gse243674_{class_slug}_{time_h}h.samples.tsv"
        run_path = workdir / f"gse243674_{class_slug}_{time_h}h.time.txt"
        out_path = workdir / f"gse243674_{class_slug}_{time_h}h.deg.tsv"
        matrix.to_csv(counts_path, sep="\t", index=False, lineterminator="\n")
        samples.to_csv(samples_path, sep="\t", index=False, lineterminator="\n")
        run_r(step / DESEQ2_SCRIPT,
              [counts_path, samples_path, step / HGNC, out_path],
              record_path=run_path)
        output = _read_r_output(out_path)
        contrast = f"{class_slug}|{time_h}h_vs_vehicle"
        frames.append(pd.DataFrame({
            "contrast": contrast,
            "gene_symbol": output.gene_symbol.astype(str),
            "log2fc": output.log2FC,
            "stat": output.stat,
            "p_value": output.p_value,
            "padj": output.padj,
            "base_mean": output.baseMean,
            "n_case": int(len(treated)),
            "n_control": int(len(control)),
        }))
        contrasts.append({"contrast": contrast, "stratum": f"{time_h}h",
                          "n_case": int(len(treated)), "n_control": int(len(control))})
        peaks.append(_r_peak(run_path))
    return pd.concat(frames, ignore_index=True), contrasts, {"r_peak_rss_kib": max(
        value for value in peaks if value is not None) if any(peaks) else None}


def read_dtoxs_counts(step):
    """The GSE217421 Conv read counts (first column = gene symbol)."""
    step = Path(step)
    with gzip.open(step / DTOXS_COUNTS, "rt", encoding="utf-8") as handle:
        wells = handle.readline().rstrip("\n").split("\t")
    frame = pd.read_csv(step / DTOXS_COUNTS, sep="\t", header=None, skiprows=1,
                        dtype={0: str})
    if frame.shape[1] != len(wells) + 1:
        raise CardiotoxError("the DToxS read-count table is ragged")
    frame.columns = ["gene"] + wells
    frame[wells] = frame[wells].astype(np.int32)
    return frame


def dtoxs_class_contrasts(step, counts, member_codes, class_slug, workdir):
    """Pooled class-vs-vehicle edgeR contrasts, one per (cell line, plate)."""
    step, workdir = Path(step), Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    strata = dtox_pooled_strata(step / DTOXS_DEG_TAR, set(member_codes))
    wells, samples = [], []
    for (cell, plate, hour), record in strata.items():
        stratum = f"{cell}|plate{plate}|{hour}h"
        for well in record["controls"]:
            samples.append({"well": well, "stratum": stratum, "group": "CTRL"})
        for code in sorted(record["treated"]):
            for well in record["treated"][code]:
                samples.append({"well": well, "stratum": stratum, "group": "CLASS"})
        wells.extend(record["controls"])
        wells.extend(well for code in sorted(record["treated"])
                     for well in record["treated"][code])
    if len(set(wells)) != len(wells):
        raise CardiotoxError("a DToxS well is pooled twice for one class")
    samples = pd.DataFrame(samples)
    subset = counts[["gene"] + wells].copy()
    counts_path = workdir / f"dtoxs_{class_slug}.counts.tsv"
    samples_path = workdir / f"dtoxs_{class_slug}.samples.tsv"
    run_path = workdir / f"dtoxs_{class_slug}.time.txt"
    out_path = workdir / f"dtoxs_{class_slug}.deg.tsv"
    subset.to_csv(counts_path, sep="\t", index=False, lineterminator="\n")
    samples.to_csv(samples_path, sep="\t", index=False, lineterminator="\n")
    library = workdir / f"dtoxs_{class_slug}.library"
    if not library.exists():
        with tarfile.open(step / DTOXS_PROGRAMS, "r|gz") as archive:
            for entry in archive:
                if not Path(entry.name).name.endswith(".R"):
                    continue
                target = library / Path(entry.name).name
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(archive.extractfile(entry).read())
    run_r(step / EDGE_R_SCRIPT,
          [counts_path, samples_path, library, out_path], record_path=run_path)
    output = _read_r_output(out_path)
    output["contrast"] = class_slug + "|" + output.stratum.astype(str)
    frame = pd.DataFrame({
        "contrast": output.contrast,
        "gene_symbol": output.gene_symbol.astype(str),
        "log2fc": output.logFC,
        "p_value": output.PValue,
        "padj": output.FDR,
        "base_mean": output.logCPM,
        "n_case": pd.to_numeric(output.n_case, errors="coerce").astype(int),
        "n_control": pd.to_numeric(output.n_control, errors="coerce").astype(int),
    })
    contrasts = [{"contrast": str(row.contrast), "stratum": str(row.stratum),
                  "n_case": int(row.n_case), "n_control": int(row.n_control)}
                 for row in output.drop_duplicates("stratum").itertuples(index=False)]
    return frame, contrasts, {"r_peak_rss_kib": _r_peak(run_path)}


def _benjamini_hochberg(pvalues):
    pvalues = np.asarray(pvalues, dtype=np.float64)
    n = len(pvalues)
    order = np.argsort(pvalues, kind="stable")
    ranked = pvalues[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty(n, dtype=np.float64)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def _welch(values, treated_mask, control_mask):
    treated, control = values[:, treated_mask], values[:, control_mask]
    n_t, n_c = treated.shape[1], control.shape[1]
    mean_t, mean_c = treated.mean(axis=1), control.mean(axis=1)
    var_t = treated.var(axis=1, ddof=1) if n_t > 1 else np.zeros_like(mean_t)
    var_c = control.var(axis=1, ddof=1) if n_c > 1 else np.zeros_like(mean_c)
    log2fc = mean_t - mean_c
    se_t, se_c = var_t / n_t, var_c / n_c
    total = se_t + se_c
    with np.errstate(invalid="ignore", divide="ignore"):
        statistic = np.where(total > 0, log2fc / np.sqrt(total), 0.0)
        degrees = np.where(total > 0,
                           total ** 2 / (se_t ** 2 / max(n_t - 1, 1)
                                         + se_c ** 2 / max(n_c - 1, 1)),
                           np.inf)
        pvalue = np.where(total > 0, 2.0 * _t_sf(np.abs(statistic), degrees), 1.0)
    return log2fc, pvalue, statistic


def _t_sf(value, degrees):
    from scipy import stats
    return stats.t.sf(value, degrees)


def gse244740_class_contrasts(step, member_names, class_slug):
    """Pooled class-vs-vehicle Welch contrasts on the ScreenSeq counts."""
    step = Path(step)
    panel = gse244740_panel(step / GSE244740_SERIES)
    members = set(member_names)
    frames, contrasts = [], []
    for plate in sorted(set(panel.plate[panel.compound.isin(members)])):
        treated = panel[(panel.plate == plate) & panel.compound.isin(members)]
        control = panel[(panel.plate == plate) & panel.compound.eq("DMSO")]
        if treated.empty or control.empty:
            continue
        wells = list(control.well) + list(treated.well)
        table = pd.read_csv(step / GSE244740_COUNTS, sep="\t",
                            usecols=["Ensembl"] + wells)
        identifiers = table.pop("Ensembl").astype(str)
        values = table[wells].to_numpy(dtype=np.float64)
        library = values.sum(axis=0)
        if (library <= 0).any():
            raise CardiotoxError("GSE244740 has an empty library")
        expression = np.log2(values / library[None, :] * 1e6 + 1.0)
        treated_mask = np.array([well in set(treated.well) for well in wells])
        control_mask = np.array([well in set(control.well) for well in wells])
        log2fc, pvalue, statistic = _welch(expression, treated_mask, control_mask)
        padj = _benjamini_hochberg(pvalue)
        contrast = f"{class_slug}|plate{plate}|24h_vs_DMSO"
        frames.append(pd.DataFrame({
            "contrast": contrast,
            "gene_symbol": identifiers,
            "log2fc": log2fc,
            "stat": statistic,
            "p_value": pvalue,
            "padj": padj,
            "base_mean": expression.mean(axis=1),
            "n_case": int(treated_mask.sum()),
            "n_control": int(control_mask.sum()),
        }))
        contrasts.append({"contrast": contrast, "stratum": f"plate{plate}",
                          "n_case": int(treated_mask.sum()),
                          "n_control": int(control_mask.sum())})
    if not frames:
        raise CardiotoxError("GSE244740 yields no pooled class contrast")
    return pd.concat(frames, ignore_index=True), contrasts, {"r_peak_rss_kib": None}


# --------------------------------------------------------------------------- #
# evidence records
# --------------------------------------------------------------------------- #
RECIPE_EVIDENCE = (
    "Resolve every drug code carried by GSE243674 (DNR/DOX/EPI/MTX/TRZ), "
    "DToxS-iPSC's GSE217421 Conv roster (54 codes; names from the pinned series "
    "matrix) and GSE244740's treated compound names to one annotation row in "
    "pert_data/Annotation/local/cm_compound_annotation.csv, keyed on "
    "(lookup_name, dataset_inchikey) where present; take the mechanism class "
    "from the ATC level-4 subgroup (single candidate, or the one candidate "
    "corroborated by mesh_classes when several ATC codes exist), require "
    "structure-backed small molecules only, and admit a class per source only "
    "at >=3 small-molecule drugs. TOP2 admission additionally requires every "
    "member's mesh_classes to contain Topoisomerase II Inhibitors. Unresolved, "
    "ambiguous, non-small-molecule and below-three codes are excluded and listed."
)
RECIPE_CATEGORY = (
    "Pool an entered class's drug wells into one Disease group against the "
    "same-plate/same-time vehicle in each source stratum, using the source's own "
    "test (DToxS GSE217421 pinned edgeR pipeline on its read counts; DESeq2 on "
    "GSE243674 counts with the item's donor-paired design; the declared Welch "
    "test on GSE244740 ScreenSeq counts). Align the per-stratum contrasts to the "
    "frozen 18,164 gene_index axis and collapse per gene by the maximum "
    "|log2FC| x -log10(padj) effect (ticket #186's deg_rank contract), then "
    "publish the compact table and the within-dataset rank "
    "rank_norm = rank / n_measured with unmeasured genes at n_measured + 1. No "
    "fusion, no axis and no view is computed."
)


def _slug(class_id):
    return str(class_id).replace(":", "-", 1).replace(":", "_")


def eligibility_table(step, crosswalk, classes, members, exclusions):
    records = []
    for source, spec in SOURCES.items():
        entered = classes[(classes.source == source)
                          & classes.decision.eq("entered")].sort_values("class_id")
        dropped = exclusions[exclusions.source == source]
        notes = []
        if source == "dtoxs_ipsc":
            notes.append("GSE217423 3'-DGE is the same experiment on another platform: "
                         "recorded as not independent and never substituted (see rejected.csv)")
        if source == "gse243674":
            notes.append("item-local contrast metadata (tokens/names) is the tracked "
                         "deg.json; the pooled DESeq2 is this step's analysis on the "
                         "immutable source counts")
        if source == "gse244740":
            notes.append("the three TOP2 compounds share one plate and one 24 h time "
                         "point; all doses are pooled into the class arm")
        notes.append("ADR 0066 declared antibodies (trastuzumab, bevacizumab, "
                     "rituximab) carry no small-molecule annotation row and are "
                     "therefore unresolved/excluded here")
        records.append({
            "source": source,
            "item": spec["item"],
            "accession": spec["accession"],
            "species": spec["species"],
            "model_system": spec["model_system"],
            "evidence_type": "bulk",
            "tier": 2,
            "shared_control_group": spec["shared_control_group"],
            "vehicle_design": spec["vehicle_design"],
            "entered_classes": "; ".join(entered.class_id),
            "n_entered_classes": int(len(entered)),
            "n_excluded_codes": int(len(dropped)),
            "excluded_reasons": "; ".join(sorted(set(dropped.reason))),
            "exploratory": False,
            "notes": " | ".join(notes),
        })
    return pd.DataFrame(records)


def rejected_table():
    return pd.DataFrame([
        {"candidate": "GSE217423 (DToxS 3'-DGE arm)",
         "source_item": DTOXS_ITEM,
         "reason": "same experiment as GSE217421 on a targeted 3'-DGE platform; "
                   "not an independent cohort, so it is not a second DToxS dataset"},
        {"candidate": "GSE242692", "source_item": "(not admitted)",
         "reason": "single-drug doxorubicin series; fails the >=3 small-molecule "
                   "class rule (ADR 0066)"},
        {"candidate": "GSE313746", "source_item": "(not admitted)",
         "reason": "single-drug doxorubicin series; fails the >=3 small-molecule "
                   "class rule (ADR 0066)"},
    ])


# --------------------------------------------------------------------------- #
# product writers
# --------------------------------------------------------------------------- #
def _output_record(path, rows):
    path = Path(path)
    return {"sha256": sha256(path), "bytes": path.stat().st_size, "rows": int(rows)}


def write_evidence_product(root, step, crosswalk, classes, members, exclusions,
                           rejected, inputs, environment, axis_digest):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    eligibility = eligibility_table(step, crosswalk, classes, members, exclusions)
    outputs = {}
    for name, frame in (("crosswalk.csv", crosswalk), ("classes.csv", classes),
                        ("class_members.csv", members), ("excluded.csv", exclusions),
                        ("eligibility.csv", eligibility), ("rejected.csv", rejected)):
        frame.to_csv(root / name, index=False, lineterminator="\n")
        outputs[name] = _output_record(root / name, len(frame))
    manifest = {
        "schema": SCHEMA,
        "product": EVIDENCE_DIR,
        "recipe": RECIPE_EVIDENCE,
        "axis": {"path": AXIS, "sha256": axis_digest, "size": AXIS_SIZE},
        "inputs": inputs,
        "outputs": outputs,
        "counts": {
            "codes": int(len(crosswalk)),
            "resolved": int((crosswalk.status == "assigned").sum()),
            "unresolved": int((crosswalk.status != "assigned").sum()),
            "classes_with_members": int(classes.class_id.nunique()),
            "entered_classes": int(classes.decision.eq("entered").sum()),
            "excluded_codes": int(len(exclusions)),
        },
        "environment": environment,
    }
    write_json(root / "manifest.json", manifest)
    return manifest


def _catalog_entry(record):
    return dict(record["entry"])


def category_tier_ledger(category, records):
    """The v2 tier ledger and level-1 ``DEG2000`` of one cardiotox class (#208).

    Every cardiotox class source is declared at tier 2, so the ladder is one rung;
    the record is still computed by the shared disease contract
    (:func:`bases.disease_deg.tier_ledger` / :func:`bases.disease_deg.fuse_disease`)
    rather than written by hand, exactly as the registry categories do. The
    per-source ``n_deg_in_deg2000`` is written back onto the catalog entries.
    """
    declared = [{"dataset_id": record["entry"]["dataset_id"],
                 "tier": int(record["entry"]["tier"])} for record in records]
    measured, prepared = {}, []
    for record in records:
        dataset_id = record["entry"]["dataset_id"]
        derived = disease_deg.source_deg_rank(record["table"])
        mask = np.asarray(derived["rank"], dtype=np.int64) <= int(derived["n_deg"])
        measured[dataset_id] = {"deg": np.flatnonzero(mask).astype(np.int64),
                                "n_deg": int(derived["n_deg"]),
                                "informative": bool(derived["informative"])}
        prepared.append({"dataset_id": dataset_id,
                         "tier": int(record["entry"]["tier"]),
                         "rank": np.asarray(derived["rank"], dtype=np.int64),
                         "n_deg": int(derived["n_deg"]),
                         "informative": bool(derived["informative"]),
                         "informativeness_reason": derived["reason"]})
    ledger = disease_deg.tier_ledger(category, declared, measured)
    fusion = disease_deg.fuse_disease(prepared, axis_size=AXIS_SIZE)
    if int(fusion["n_deg"]) != int(ledger["n_deg"]):
        raise CardiotoxError(
            f"{category}: the tier ledger union ({ledger['n_deg']}) differs from "
            f"the fusion union ({fusion['n_deg']})")
    by_id = {row["dataset_id"]: row for row in fusion["sources"]}
    for record in records:
        record["entry"]["n_deg_in_deg2000"] = int(
            by_id[record["entry"]["dataset_id"]]["n_deg_in_deg2000"])
    return ledger, fusion


def write_category_product(root, category, records, inputs, environment,
                           axis_digest):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    outputs = {}
    ledger, fusion = category_tier_ledger(category, records)
    for record in records:
        write_deg_table(root / record["deg"], record["table"], record["axis"])
        write_rank_table(root / record["rank_name"], record["rank"])
        outputs[record["deg"]] = _output_record(root / record["deg"],
                                                len(record["table"]))
        outputs[record["rank_name"]] = _output_record(root / record["rank_name"],
                                                      AXIS_SIZE)
    catalog = {
        "schema": CATALOG_SCHEMA,
        "category": category,
        "axis": {"path": AXIS, "sha256": axis_digest, "size": AXIS_SIZE},
        "tier": 2,
        "evidence_type": "bulk",
        "datasets": [_catalog_entry(record) for record in records],
        # ADR 0071 decisions 2/4 (#208): the declared tier ladder, the
        # axis-resolved stop record and the level-1 DEG union/`DEG2000` sizes.
        "n_deg": int(fusion["n_deg"]),
        "deg2000": int(fusion["deg2000"]),
        "tiering": ledger["tiering"],
        "tier_sources": ledger["tier_sources"],
    }
    # A JSON file cannot contain its own file digest. ``catalog_sha256`` is the
    # digest of the canonical catalog bytes *without* that field, so a consumer
    # can strip the field, re-serialise with write_json's rule and compare.
    write_json(root / "catalog.json", catalog)
    catalog_sha256 = sha256(root / "catalog.json")
    catalog["catalog_sha256"] = catalog_sha256
    write_json(root / "catalog.json", catalog)
    outputs["catalog.json"] = _output_record(root / "catalog.json", len(records))
    measured = sum(record["entry"]["n_measured"] for record in records)
    manifest = {
        "schema": SCHEMA,
        "product": category,
        "recipe": RECIPE_CATEGORY,
        "axis": {"path": AXIS, "sha256": axis_digest, "size": AXIS_SIZE},
        "inputs": inputs,
        "outputs": outputs,
        "counts": {"datasets": len(records), "measured_gene_rows": int(measured),
                   "n_deg": int(fusion["n_deg"]),
                   "deg2000": int(fusion["deg2000"])},
        "catalog_sha256": catalog_sha256,
        "environment": environment,
    }
    write_json(root / "manifest.json", manifest)
    return manifest


# --------------------------------------------------------------------------- #
# dataset assembly
# --------------------------------------------------------------------------- #
def _arm_counts(contrasts, key):
    """Numeric arm size for one pooled class dataset and how it was derived.

    A class DEG pools drug wells against the same-plate vehicle in one or more
    source strata. When the per-stratum arm sizes are constant the catalog value
    is that arm size (the per-contrast pooled-well count); when they differ the
    catalog value is the sum over the disjoint strata, so the number stays the
    real count of pooled wells instead of a null. Either way the derivation is
    recorded in ``*_basis``.
    """
    values = [int(item[key]) for item in contrasts]
    if not values:
        raise CardiotoxError("a class dataset declares no contrast")
    if len(set(values)) == 1:
        return values[0], (f"per-contrast pooled wells, constant across "
                           f"{len(values)} contrast(s)")
    return sum(values), (f"sum of pooled wells over {len(values)} disjoint source "
                         f"strata (per-contrast counts in contrasts[])")


def _dataset_entry(source, class_id, class_label, member_drugs, member_source_ids,
                   contrasts, stats, deg_name, rank_name, inputs, test_method):
    spec = SOURCES[source]
    n_case, n_case_basis = _arm_counts(contrasts, "n_case")
    n_control, n_control_basis = _arm_counts(contrasts, "n_control")
    return {
        "dataset_id": source,
        "category": None,
        "class_id": class_id,
        "class_label": class_label,
        "class_slug": _slug(class_id),
        "source_item": spec["item"],
        "accession": spec["accession"],
        "species": spec["species"],
        "model_system": spec["model_system"],
        "evidence_type": "bulk",
        "tier": 2,
        "shared_control_group": spec["shared_control_group"],
        "member_drugs": sorted(member_drugs),
        "n_small_molecule_drugs": len(member_drugs),
        "member_source_ids": sorted(member_source_ids),
        "n_case": n_case,
        "n_control": n_control,
        "n_case_basis": n_case_basis,
        "n_control_basis": n_control_basis,
        "contrasts": contrasts,
        "n_contrasts": int(len(contrasts)),
        "n_measured": int(stats["n_measured"]),
        "n_deg": int(stats["n_deg"]),
        "n_significant": int(stats["n_significant"]),
        # ADR 0068: all-tied effect scores carry no ordering evidence; the record
        # is still published but the gene-axis fusion excludes it.
        "informative": bool(stats["informative"]),
        "informativeness": {
            "rule": stats["informativeness_rule"],
            "n_distinct_effect": int(stats["n_distinct_effect"]),
            "reason": stats["informativeness_reason"],
        },
        "alignment": stats["alignment"],
        "method": test_method,
        "deg": deg_name,
        "rank": rank_name,
        "inputs": inputs,
    }


REJECTED = rejected_table()


def _build_top2(step, resolver, axis, classes, members, workdir, inputs):
    records = []
    class_id = f"atc_l4:{TOP2_CLASS_KEY}"
    label = classes[(classes.class_id == class_id)
                    & classes.decision.eq("entered")].class_label.iloc[0]
    dtoxs_counts = read_dtoxs_counts(step)
    for source in ("gse243674", "dtoxs_ipsc", "gse244740"):
        part = members[(members.source == source) & (members.class_id == class_id)]
        if part.empty:
            continue
        slug = _slug(class_id)
        if source == "gse243674":
            frame, contrasts, peaks = gse243674_class_contrasts(
                step, part.source_id.tolist(), slug, workdir)
            method = {
                "name": "GSE243674 pooled DESeq2",
                "implementation": "cardiotox_deseq2.R",
                "env": "conda:biomni_e1 (R 4.5.3, DESeq2 1.50.2)",
                "deposition": "GSE243674 featureCounts gene counts",
                "model": "DESeq2 ~ donor + group",
                "design": "a class's member wells pooled into one CLASS arm vs the "
                          "same-time vehicle wells; donor is a fixed effect",
                "test": "Wald, unshrunken MLE log2FC, BH padj, DESeq2 independent "
                        "filtering (alpha 0.05)",
                "collapse": "per-gene maximum effect over the time-point contrasts",
                "source_item": GSE243674_ITEM,
            }
        elif source == "dtoxs_ipsc":
            frame, contrasts, peaks = dtoxs_class_contrasts(
                step, dtoxs_counts, part.source_id.tolist(), slug, workdir)
            method = {
                "name": "DToxS pinned edgeR pooled class contrast",
                "implementation": "cardiotox_edgeR.R (Program-Source-Codes pipeline)",
                "env": "conda:biomni_e1 (R 4.5.3, edgeR 4.8.2)",
                "deposition": "GSE217421 Conv read counts + shipped per-drug Calc "
                              "tables (wells/QC)",
                "model": "the deposit's pinned edgeR differential-comparison "
                         "pipeline (Program-Source-Codes)",
                "design": "a class's member wells pooled into one CLASS arm vs the "
                          "same cell line's DMSO CTRL wells per (cell line, plate) "
                          "stratum",
                "test": "edgeR exactTest on TMM-normalized counts with trended "
                        "dispersion, remove_small filter, BH FDR (alpha 0.1)",
                "collapse": "per-gene maximum effect over the per-stratum pooled "
                            "contrasts",
                "source_item": DTOXS_ITEM,
            }
        else:
            frame, contrasts, peaks = gse244740_class_contrasts(
                step, [str(value) for value in part.drug_name], slug)
            method = {
                "name": "GSE244740 pooled Welch class contrast",
                "implementation": "bases.cardiotox.gse244740_class_contrasts",
                "env": "arc python (numpy 2.4.6)",
                "deposition": "GSE244740 ScreenSeq processed counts",
                "model": "two-sample Welch t-test on log2(CPM + 1)",
                "design": "a class's member wells (all doses) pooled into one arm "
                          "vs same-plate DMSO wells at 24 h",
                "test": "Welch t per gene, BH within contrast, library over all "
                        "source genes",
                "collapse": "per-gene maximum effect over the per-plate contrasts",
                "source_item": GSE244740_ITEM,
            }
        table, rank, stats = class_deg(frame, resolver)
        deg_name = f"{source}.{slug}.deg.csv"
        rank_name = f"{source}.{slug}.rank.csv"
        source_inputs = dict(inputs)
        entry = _dataset_entry(source, class_id, label, part.drug_name.tolist(),
                               part.source_id.tolist(), contrasts, stats,
                               deg_name, rank_name, source_inputs, method)
        entry["category"] = CATEGORY_TOP2
        records.append({"deg": deg_name, "rank_name": rank_name, "table": table,
                        "rank": rank, "axis": axis, "entry": entry,
                        "r_peak": peaks.get("r_peak_rss_kib")})
    return records


def _build_other(step, resolver, axis, classes, members, workdir, inputs):
    records = []
    dtoxs_counts = read_dtoxs_counts(step)
    entered = classes[(classes.source == "dtoxs_ipsc")
                      & classes.decision.eq("entered")]
    for row in entered.sort_values("class_id").itertuples(index=False):
        if row.class_id == f"atc_l4:{TOP2_CLASS_KEY}":
            continue
        part = members[(members.source == "dtoxs_ipsc")
                       & (members.class_id == row.class_id)]
        slug = _slug(row.class_id)
        frame, contrasts, peaks = dtoxs_class_contrasts(
            step, dtoxs_counts, part.source_id.tolist(), slug, workdir)
        table, rank, stats = class_deg(frame, resolver)
        dataset_id = f"dtoxs_ipsc__{slug}"
        deg_name = f"dtoxs_ipsc.{slug}.deg.csv"
        rank_name = f"dtoxs_ipsc.{slug}.rank.csv"
        method = {
            "name": "DToxS pinned edgeR pooled class contrast",
            "implementation": "cardiotox_edgeR.R (Program-Source-Codes pipeline)",
            "env": "conda:biomni_e1 (R 4.5.3, edgeR 4.8.2)",
            "deposition": "GSE217421 Conv read counts + shipped per-drug Calc "
                          "tables (wells/QC)",
            "model": "the deposit's pinned edgeR differential-comparison pipeline "
                     "(Program-Source-Codes)",
            "design": "a class's member wells pooled into one CLASS arm vs the "
                      "same cell line's DMSO CTRL wells per (cell line, plate) "
                      "stratum",
            "test": "edgeR exactTest on TMM-normalized counts with trended "
                    "dispersion, remove_small filter, BH FDR (alpha 0.1)",
            "collapse": "per-gene maximum effect over the per-stratum pooled "
                        "contrasts",
            "source_item": DTOXS_ITEM,
        }
        entry = _dataset_entry("dtoxs_ipsc", row.class_id, row.class_label,
                               part.drug_name.tolist(), part.source_id.tolist(),
                               contrasts, stats, deg_name, rank_name,
                               dict(inputs), method)
        entry["category"] = CATEGORY_OTHER
        entry["dataset_id"] = dataset_id
        records.append({"deg": deg_name, "rank_name": rank_name, "table": table,
                        "rank": rank, "axis": axis, "entry": entry,
                        "r_peak": peaks.get("r_peak_rss_kib")})
    return records


# --------------------------------------------------------------------------- #
# build, publish and replay
# --------------------------------------------------------------------------- #
PRODUCTS = (EVIDENCE_DIR, CATEGORY_TOP2, CATEGORY_OTHER)

COMMON_INPUTS = (AXIS, CROSSWALK, HGNC, ANNOTATION)
EVIDENCE_INPUTS = COMMON_INPUTS + (DTOXS_SERIES, DTOXS_PARAMS,
                                   GSE243674_DEG_JSON, GSE244740_SERIES)
DTOXS_INPUTS = (DTOXS_SERIES, DTOXS_PARAMS, DTOXS_COUNTS, DTOXS_DEG_TAR,
                DTOXS_PROGRAMS, EDGE_R_SCRIPT)
TOP2_INPUTS = COMMON_INPUTS + DTOXS_INPUTS + (
    GSE243674_COUNTS, GSE243674_DEG_JSON, DESEQ2_SCRIPT,
    GSE244740_COUNTS, GSE244740_SERIES)
OTHER_INPUTS = COMMON_INPUTS + DTOXS_INPUTS


def _inputs(step, repo, paths):
    inputs = {}
    for path in paths:
        add_input(inputs, Path(step) / path, repo)
    return inputs


def _python_peak():
    import resource
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def build_all(root, repo):
    """Build all three cardiotoxicity products into ``root`` (a candidate dir).

    ``repo`` is the repository root; every pinned input path in this module is
    repo-relative.
    """
    import shutil

    root = Path(root)
    repo = Path(repo)
    step = repo
    workdir = root / "_work"
    workdir.mkdir(parents=True, exist_ok=True)
    try:
        axis, resolver = load_resolver(step)
        axis_digest = sha256(step / AXIS)
        environment = python_environment()
        r_environment_info = r_environment()
        annotation = read_annotation(step / ANNOTATION)
        crosswalk = build_crosswalk(step, annotation)
        classes, members, exclusions = select_classes(crosswalk)
        rejected = rejected_table()
        evidence_inputs = _inputs(step, repo, EVIDENCE_INPUTS)
        evidence = write_evidence_product(
            root / EVIDENCE_DIR, step, crosswalk, classes, members, exclusions,
            rejected, evidence_inputs, environment, axis_digest)

        top2_inputs = _inputs(step, repo, TOP2_INPUTS)
        top2_records = _build_top2(step, resolver, axis, classes, members,
                                   workdir, top2_inputs)
        other_inputs = _inputs(step, repo, OTHER_INPUTS)
        other_records = _build_other(step, resolver, axis, classes, members,
                                     workdir, other_inputs)
        category_environment = dict(environment, r=r_environment_info)
        python_peak = _python_peak()
        r_peaks = [record["r_peak"] for record in top2_records + other_records
                   if record.get("r_peak")]
        resources = {"python_peak_rss_kib": python_peak,
                     "r_peak_rss_kib": max(r_peaks) if r_peaks else None}
        top2 = write_category_product(root / CATEGORY_TOP2, CATEGORY_TOP2,
                                      top2_records, top2_inputs,
                                      category_environment, axis_digest)
        other = write_category_product(root / CATEGORY_OTHER, CATEGORY_OTHER,
                                       other_records, other_inputs,
                                       category_environment, axis_digest)
        return ({EVIDENCE_DIR: evidence, CATEGORY_TOP2: top2,
                 CATEGORY_OTHER: other}, resources)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def _new_candidate(staging):
    import tempfile
    staging = Path(staging)
    staging.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix=".disease-deg-cardiotox-", dir=staging))


def publish_release(repo, staging):
    """Build, verify and replace the three accepted products (ADR 0062)."""
    import shutil
    repo = Path(repo)
    step = repo
    candidate_root = _new_candidate(staging)
    try:
        manifests, resources = build_all(candidate_root, repo)
        for product in PRODUCTS:
            check_manifest(candidate_root / product,
                           read_json(candidate_root / product / "manifest.json"))
        output_root = repo / OUTPUT_ROOT
        for product in PRODUCTS:
            replace(candidate_root / product, output_root / product)
        return manifests, resources
    finally:
        shutil.rmtree(candidate_root, ignore_errors=True)


def replay_release(repo, staging):
    """Rebuild from the pinned inputs and compare every accepted digest."""
    import shutil
    repo = Path(repo)
    step = repo
    candidate_root = _new_candidate(staging)
    try:
        _, resources = build_all(candidate_root, repo)
        output_root = repo / OUTPUT_ROOT
        summary = {}
        for product in PRODUCTS:
            accepted_manifest = output_root / product / "manifest.json"
            accepted = read_json(accepted_manifest)
            check_manifest(output_root / product, accepted)
            replayed = read_json(candidate_root / product / "manifest.json")
            if replayed != accepted:
                raise CardiotoxError(
                    f"{product}: replayed manifest differs from the accepted release")
            check_manifest(candidate_root / product, accepted)
            summary[product] = {
                "manifest_sha256": sha256(accepted_manifest),
                "outputs": sorted(output_digests(accepted)),
            }
        return {"result": "PASS", **summary}, resources
    finally:
        shutil.rmtree(candidate_root, ignore_errors=True)
