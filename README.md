# IBDEX projector

Project raw intestinal bulk RNA-seq count matrices into the frozen IBDEX
(IBD Explorer) latent reference.

The command-line tool performs the full public projection path:

1. Convert gene identifiers to HGNC symbols.
2. Collapse duplicate symbols by summing raw counts.
3. Convert raw counts to log2(CPM + 1).
4. Restrict to the frozen 3,000-gene IBDEX panel.
5. Fill missing panel genes at the training mean.
6. Apply training z-score normalization and clipping.
7. Encode samples with the frozen CVAE using dataset condition zero.
8. Compute calibrated latent pathway scores.
9. Write an explorer-ready JSON file.

## Install locally

```bash
pip install -e .
```

## Input format

Counts may be genes x samples or samples x genes. If genes are rows, the first
column should contain gene identifiers. Identifiers may be HGNC symbols, Ensembl
gene IDs with or without version suffixes, NCBI Entrez GeneIDs, or a mixture.

Metadata must contain one row per sample and a tissue column. Accepted tissue
values include `colon`, `ileum`, `left_colon`, `right_colon`, `rectum`,
`sigmoid`, `transverse_colon` and `colon_unspecified`.

## Use

```bash
ibdex project \
  --counts examples/counts_smoke.csv \
  --metadata examples/metadata.csv \
  --sample-col sample_id \
  --tissue-col tissue \
  --out ibdex_projection.json
```

The JSON can be loaded by the IBDEX HTML explorer.

## Important interpretation

External samples are projected into the frozen reference without retraining.
All dataset condition indicators are set to zero. This matches the public
projection setting: the incoming cohort is not one of the four training datasets.

Pathway scores are latent decoder-derived pathway coordinates, standardized
against the internal IBDEX reference distribution. They are not ssGSEA scores.
