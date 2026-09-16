from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from .artifacts import artifact_path
from .gene_mapping import GeneResolver


def load_counts(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def counts_to_samples_by_symbol(
    counts: pd.DataFrame,
    resolver: GeneResolver | None = None,
    gene_col: str | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Return samples x HGNC-symbol raw counts."""

    resolver = resolver or GeneResolver()
    matrix, genes_are_rows = _orient_counts(counts, gene_col=gene_col)

    if genes_are_rows:
        genes = matrix.index
        values = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    else:
        genes = matrix.columns
        values = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0).T

    symbols, mapping_report = resolver.resolve_many(genes)
    values = values.copy()
    values["__symbol__"] = symbols
    values = values.dropna(subset=["__symbol__"])
    collapsed = values.groupby("__symbol__", sort=False).sum(numeric_only=True)
    samples_by_gene = collapsed.T
    report = asdict(mapping_report)
    report["unique_symbols_after_collapse"] = int(samples_by_gene.shape[1])
    report["orientation"] = "genes_x_samples" if genes_are_rows else "samples_x_genes"
    return samples_by_gene, report


def preprocess_counts(
    counts: pd.DataFrame,
    gene_col: str | None = None,
    annotation_path: str | Path | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Convert raw counts to normalized IBDEX model input."""

    resolver = GeneResolver(annotation_path)
    raw, gene_report = counts_to_samples_by_symbol(counts, resolver=resolver, gene_col=gene_col)
    raw = raw.clip(lower=0)
    library_sizes = raw.sum(axis=1).replace(0, np.nan)
    cpm = raw.div(library_sizes, axis=0) * 1_000_000.0
    log_cpm = np.log2(cpm.fillna(0.0) + 1.0)

    stats = pd.read_csv(artifact_path("normalization_stats_train_only.csv"))
    panel = stats["gene"].astype(str).tolist()
    mean = stats.set_index("gene")["train_mean_pre_zscore"].astype(float)
    std = stats.set_index("gene")["std_used"].astype(float).replace(0, 1.0)

    aligned = pd.DataFrame(index=log_cpm.index, columns=panel, dtype=float)
    common = [g for g in panel if g in log_cpm.columns]
    aligned.loc[:, common] = log_cpm.loc[:, common]
    missing = [g for g in panel if g not in log_cpm.columns]
    for gene in missing:
        aligned[gene] = mean.loc[gene]

    normalized = (aligned[panel] - mean.loc[panel]) / std.loc[panel]
    normalized = normalized.clip(lower=-5.0, upper=5.0).astype("float32")
    report = {
        **gene_report,
        "n_samples": int(normalized.shape[0]),
        "n_panel_genes": len(panel),
        "n_panel_genes_present": len(common),
        "n_panel_genes_missing": len(missing),
        "missing_panel_genes": missing,
    }
    return normalized, report


def _orient_counts(counts: pd.DataFrame, gene_col: str | None = None) -> tuple[pd.DataFrame, bool]:
    if gene_col:
        if gene_col not in counts.columns:
            raise ValueError(f"gene column not found: {gene_col}")
        matrix = counts.set_index(gene_col)
        return matrix, True

    first = counts.columns[0]
    first_lc = str(first).strip().lower()
    numeric_rest = counts.drop(columns=[first]).apply(pd.to_numeric, errors="coerce")
    if first_lc in {"sample", "sample_id", "sampleid", "id", "barcode"} and numeric_rest.notna().mean().mean() > 0.8:
        return counts.set_index(first), False

    first_col = counts[first].astype(str)
    first_looks_gene = first_col.str.match(r"^(ENSG\d+|[A-Za-z][A-Za-z0-9_.-]*|\d+)(\.\d+)?$").mean()
    if first_looks_gene > 0.8 and numeric_rest.notna().mean().mean() > 0.8:
        return counts.set_index(first), True

    numeric_all = counts.apply(pd.to_numeric, errors="coerce")
    if numeric_all.notna().mean().mean() > 0.8:
        return counts, False
    raise ValueError("Could not infer count-matrix orientation. Pass --gene-col for genes x samples input.")
