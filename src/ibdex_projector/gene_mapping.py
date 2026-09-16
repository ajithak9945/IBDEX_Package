from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

from .artifacts import artifact_path


ENSEMBL_RE = re.compile(r"^ENSG\d+(?:\.\d+)?$", re.IGNORECASE)


@dataclass
class GeneMappingReport:
    total_input_genes: int
    mapped_genes: int
    unmapped_genes: int
    duplicated_symbols_after_mapping: int


class GeneResolver:
    """Resolve Ensembl, Entrez and symbol-like identifiers to HGNC symbols."""

    def __init__(self, annotation_path: str | Path | None = None):
        self.annotation_path = Path(annotation_path) if annotation_path else artifact_path(
            "Human.GRCh38.p13.annot.tsv.gz"
        )
        self._mapping = self._load_mapping()

    def _load_mapping(self) -> dict[str, str]:
        cols = ["GeneID", "Symbol", "Synonyms", "EnsemblGeneID"]
        ann = pd.read_csv(self.annotation_path, sep="\t", usecols=cols, dtype=str)
        mapping: dict[str, str] = {}
        for row in ann.itertuples(index=False):
            symbol = _clean_symbol(row.Symbol)
            if not symbol:
                continue
            mapping[symbol.upper()] = symbol
            if isinstance(row.GeneID, str) and row.GeneID:
                mapping[row.GeneID.strip()] = symbol
            if isinstance(row.EnsemblGeneID, str) and row.EnsemblGeneID:
                for ens in re.split(r"[|;, ]+", row.EnsemblGeneID):
                    ens = _strip_ensembl_version(ens)
                    if ens:
                        mapping[ens.upper()] = symbol
            if isinstance(row.Synonyms, str) and row.Synonyms:
                for alias in re.split(r"[|;,]", row.Synonyms):
                    alias = _clean_symbol(alias)
                    if alias:
                        mapping.setdefault(alias.upper(), symbol)
        return mapping

    def resolve_many(self, gene_ids: Iterable[object]) -> tuple[list[str | None], GeneMappingReport]:
        resolved: list[str | None] = []
        input_ids = [str(g).strip() for g in gene_ids]
        for gid in input_ids:
            resolved.append(self.resolve(gid))
        mapped = [x for x in resolved if x]
        duplicated = int(pd.Series(mapped).duplicated().sum()) if mapped else 0
        report = GeneMappingReport(
            total_input_genes=len(input_ids),
            mapped_genes=len(mapped),
            unmapped_genes=len(input_ids) - len(mapped),
            duplicated_symbols_after_mapping=duplicated,
        )
        return resolved, report

    def resolve(self, gene_id: str) -> str | None:
        raw = str(gene_id).strip()
        if not raw:
            return None
        candidates = [raw, raw.upper()]
        if ENSG := _strip_ensembl_version(raw):
            candidates.append(ENSG.upper())
        if "|" in raw:
            candidates.extend(part.strip() for part in raw.split("|"))
        for cand in candidates:
            if cand in self._mapping:
                return self._mapping[cand]
        if _looks_like_symbol(raw):
            return _clean_symbol(raw)
        return None


def _strip_ensembl_version(value: str) -> str | None:
    value = str(value).strip()
    if not value:
        return None
    if "|" in value:
        for part in value.split("|"):
            stripped = _strip_ensembl_version(part)
            if stripped:
                return stripped
    match = ENSEMBL_RE.match(value)
    if match:
        return value.split(".")[0]
    return None


def _clean_symbol(value: object) -> str | None:
    if not isinstance(value, str):
        return None
    value = value.strip()
    if not value or value.lower() == "nan":
        return None
    return value


def _looks_like_symbol(value: str) -> bool:
    if ENSEMBL_RE.match(value):
        return False
    if value.isdigit():
        return False
    return bool(re.match(r"^[A-Za-z0-9_.-]+$", value))
