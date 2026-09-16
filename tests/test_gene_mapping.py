import pandas as pd

from ibdex_projector.preprocess import counts_to_samples_by_symbol


def test_ensembl_versions_are_collapsed_to_symbols():
    counts = pd.DataFrame({
        "gene_id": ["ENSG00000169245.6", "ENSG00000169245", "CXCL10"],
        "sample_001": [1, 2, 3],
    })
    out, report = counts_to_samples_by_symbol(counts, gene_col="gene_id")
    assert "CXCL10" in out.columns
    assert out.loc["sample_001", "CXCL10"] == 6
    assert report["duplicated_symbols_after_mapping"] >= 1
