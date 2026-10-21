"""Load the pinned annotation assets (HGNC + compound registry snapshot).

Both files live under ``pert_data/annotation/`` and are produced by
``fetch_annotation.py``.  The build only reads them; it never calls the
network.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from . import paths
from .text import normalize_name


class AnnotationError(RuntimeError):
    """The pinned annotation assets are missing or malformed."""


# ---------------------------------------------------------------------------
# HGNC
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class GeneRef:
    hgnc_id: str
    symbol: str
    entrez_id: str
    ensembl_id: str


class HgncIndex:
    """Current-symbol, alias/previous-symbol, Entrez and Ensembl indices."""

    def __init__(self, table: pd.DataFrame):
        self.genes = table
        self._symbol: dict[str, str] = {}
        self._alias: dict[str, list[str]] = {}
        self._previous: dict[str, list[str]] = {}
        self._entrez: dict[str, list[str]] = {}
        self._ensembl: dict[str, list[str]] = {}
        self._by_id: dict[str, GeneRef] = {}

        for row in table.itertuples(index=False):
            hgnc_id = row.hgnc_id
            self._symbol[row.symbol] = hgnc_id
            self._by_id[hgnc_id] = GeneRef(
                hgnc_id=hgnc_id, symbol=row.symbol,
                entrez_id=str(row.entrez_id) if pd.notna(row.entrez_id) else "",
                ensembl_id=(str(row.ensembl_gene_id)
                            if pd.notna(row.ensembl_gene_id) else ""))
            if row.alias_symbol:
                for alias in str(row.alias_symbol).split("|"):
                    if alias:
                        self._alias.setdefault(alias, []).append(hgnc_id)
            if row.prev_symbol:
                for previous in str(row.prev_symbol).split("|"):
                    if previous:
                        self._previous.setdefault(previous, []).append(hgnc_id)
            if row.entrez_id:
                self._entrez.setdefault(str(row.entrez_id), []).append(hgnc_id)
            if row.ensembl_gene_id:
                self._ensembl.setdefault(str(row.ensembl_gene_id), []).append(hgnc_id)

    @classmethod
    def load(cls, path=paths.HGNC_PATH) -> "HgncIndex":
        if not path.exists():
            raise AnnotationError(
                f"HGNC asset missing: {path}. Run fetch_annotation.py --hgnc")
        table = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        required = {"hgnc_id", "symbol", "alias_symbol", "prev_symbol",
                    "entrez_id", "ensembl_gene_id"}
        missing = required - set(table.columns)
        if missing:
            raise AnnotationError(f"HGNC asset missing columns: {sorted(missing)}")
        return cls(table)

    def ref(self, hgnc_id: str) -> GeneRef:
        return self._by_id[hgnc_id]

    def _unique(self, mapping: dict[str, list[str]], key: str) -> tuple[str | None, str]:
        hits = mapping.get(key, [])
        if len(hits) == 1:
            return hits[0], "ok"
        if len(hits) > 1:
            return None, "ambiguous"
        return None, "unmapped"

    def resolve_symbol(self, symbol: str) -> tuple[str | None, str]:
        if not symbol:
            return None, "unmapped"
        if symbol in self._symbol:
            return self._symbol[symbol], "direct"
        hgnc_id, status = self._unique(self._previous, symbol)
        if hgnc_id:
            return hgnc_id, "previous"
        hgnc_id, status = self._unique(self._alias, symbol)
        if hgnc_id:
            return hgnc_id, "alias"
        return None, status

    def resolve_entrez(self, entrez: str) -> tuple[str | None, str]:
        if not entrez:
            return None, "unmapped"
        return self._unique(self._entrez, str(entrez))

    def resolve_ensembl(self, ensembl: str) -> tuple[str | None, str]:
        if not ensembl:
            return None, "unmapped"
        return self._unique(self._ensembl, str(ensembl).split(".")[0])


# ---------------------------------------------------------------------------
# Compounds
# ---------------------------------------------------------------------------
# The only manual pharmacology supplement in corpus v1.  It applies only to
# the two curated panels (DToxS / screen_ctox); temp_env (GSE262419) stays
# automatic.
MANUAL_SUPPLEMENT = {
    # SB-431542: a real ALK5/TGF-betaR1 inhibitor missing from ChEMBL's
    # mechanism table at snapshot time.
    "FHYUGAJXYORMHI-UHFFFAOYSA-N": {
        "moa": "TGF-beta receptor type I (ALK5) inhibitor",
        "citation": ("Inman GJ, Nicolas FJ, Callahan JF, et al. "
                     "SB-431542 is a potent and specific inhibitor of "
                     "transforming growth factor-beta superfamily type I "
                     "activin receptor-like kinase (ALK) receptors ALK4, "
                     "ALK5, and ALK7. Mol Pharmacol 2002;62(1):65-74. "
                     "PMID 12065753"),
    },
}
MANUAL_SOURCE_PREFIX = "manual:"


class CompoundAnnotation:
    """Frozen public-registry structure + pharmacology table."""

    def __init__(self, table: pd.DataFrame):
        table = table.copy()
        table["lookup_name"] = table["lookup_name"].map(normalize_name)
        table["dataset_inchikey"] = table["dataset_inchikey"].fillna("")
        self.table = table
        self._index: dict[tuple[str, str], dict] = {}
        for row in table.to_dict("records"):
            key = (row["lookup_name"], row["dataset_inchikey"])
            self._index[key] = row

    @classmethod
    def load(cls, path=paths.COMPOUND_ANNOTATION_PATH) -> "CompoundAnnotation":
        if not path.exists():
            raise AnnotationError(
                f"compound annotation asset missing: {path}. "
                "Run fetch_annotation.py --compounds")
        return cls(pd.read_csv(path, dtype=str, keep_default_na=False))

    def get(self, name: str, dataset_inchikey: str = "") -> dict | None:
        """Look up one compound.

        A structure-bearing row wins; a name-only row is used only when no
        structure row exists.  This is what keeps DToxS ``AMI`` (a
        pentyl-benzofuran homolog) distinct from name-based amiodarone.
        """
        name = normalize_name(name)
        if dataset_inchikey:
            row = self._index.get((name, dataset_inchikey))
            if row is not None:
                return row
        return self._index.get((name, ""))


def load_assets() -> tuple[HgncIndex, CompoundAnnotation]:
    return HgncIndex.load(), CompoundAnnotation.load()
