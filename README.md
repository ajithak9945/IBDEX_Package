# IBDEX projector

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-3.9%2B-blue)
[![Explorer](https://img.shields.io/badge/explorer-ibdexmd.org-0f766e)](https://ibdexmd.org)

Project raw intestinal bulk RNA-seq counts into IBDEX, a frozen latent
reference built from 3,168 harmonised colon and ileum samples across four
IBD cohorts. A tissue-conditioned variational autoencoder places each new
sample in that reference space without retraining, and a decoder-native
scoring step reads off 16 calibrated pathway scores directly from the
latent geometry.

## Links

- Explorer at [ibdexmd.org](https://ibdexmd.org)
- Source on [github.com/ajithak9945/IBDEX_Package](https://github.com/ajithak9945/IBDEX_Package)
- Contact at [ajithak9945@gmail.com](mailto:ajithak9945@gmail.com)

## What it does

1. Convert gene identifiers to HGNC symbols.
2. Collapse duplicate symbols by summing raw counts.
3. Convert raw counts to log2(CPM + 1).
4. Restrict to the frozen 3,000-gene IBDEX panel.
5. Fill missing panel genes at the training mean.
6. Apply training z-score normalization and clipping.
7. Encode samples with the frozen CVAE using dataset condition zero.
8. Compute calibrated latent pathway scores.
9. Place each sample onto the reference UMAP layout and assign it to its
   nearest discovered endotype.
10. Write an explorer-ready JSON file.

## Model

| | |
|---|---|
| Architecture | tissue-conditioned variational autoencoder, gradient-reversal adversaries on dataset and tissue-site identity |
| Latent dimensions | 16 |
| Training cohort | 3,168 samples, 4 harmonised IBD cohorts (GSE193677, IBDome, HMP2, GSE165512) |
| Gene panel | 3,000 genes, multi-evidence weighted selection (STRING, Reactome, KEGG, GO-BP, Open Targets, GWAS) |
| Discovered endotypes | 4 colon (C1-C4), 4 ileum (I1-I4), consensus clustering over four base algorithms |
| Pathway axes | 16, decoder-native and tissue-specific |

## Prerequisites

- Python 3.9 or newer
- pip
- git, to clone the repository

The package pulls in its own dependencies on install, including numpy >= 1.23,
pandas >= 1.5 and torch >= 2.0. A GPU is not required; the frozen model runs
on CPU (`--device cpu`, the default).

## Install

```bash
git clone https://github.com/ajithak9945/IBDEX_Package.git
cd IBDEX_Package
pip install -e .
```

This installs the `ibdex` command and the `ibdex_projector` Python package,
along with the frozen model weights and reference data bundled under
`src/ibdex_projector/artifacts/`.

## Prepare your input files

Two files are needed, a raw count matrix and a sample metadata table.

**Counts file.** A CSV or TSV of raw (unnormalised) counts, either genes as
rows and samples as columns, or samples as rows and genes as columns.
Orientation is detected automatically from the first column, whether it
looks like a sample identifier, a gene identifier, or the sheet has no ID
column at all. Pass `--gene-col <column name>` yourself only if that
detection fails or you want to force genes-as-rows explicitly. Gene
identifiers can be HGNC symbols, Ensembl gene IDs (with or without a
version suffix, for example `ENSG00000135100.10`), NCBI Entrez GeneIDs, or
a mixture of these; the tool resolves them all to HGNC symbols
automatically and sums counts for any duplicates that collapse onto the
same symbol.

**Metadata file.** A CSV with one row per sample. Two columns matter most.

| Flag | Default column name | What it holds |
|---|---|---|
| `--sample-col` | `sample_id` | must match the sample identifiers used in the counts file |
| `--tissue-col` | `tissue` | one of `colon`, `ileum`, `left_colon`, `right_colon`, `rectum`, `sigmoid`, `transverse_colon`, `colon_unspecified` |

If your files already use these column names, you do not need to pass the
flags at all. If they use different names, for example `patient_id` or
`biopsy_site`, point the tool at them with `--sample-col patient_id
--tissue-col biopsy_site`; no renaming of your files is needed.

Two further metadata columns are optional and passed straight through into
the output JSON, where the explorer shows them on the patient card.

- `disease` (also accepted as `Disease` or `diagnosis`)
- `cohort` (also accepted as `Cohort` or `dataset`)

## Run it

### From the command line

```bash
ibdex project \
  --counts your_counts.csv \
  --metadata your_metadata.csv \
  --sample-col sample_id \
  --tissue-col tissue \
  --out ibdex_projection.json
```

`examples/counts_smoke.csv` and `examples/metadata.csv` are included as a
working example; run the command above against those two files first to
confirm the install works before pointing it at your own data.

### From Python

```python
from ibdex_projector.project import project_counts_to_json

project_counts_to_json(
    "your_counts.csv",
    "your_metadata.csv",
    sample_col="sample_id",
    tissue_col="tissue",
    out="ibdex_projection.json",
)
```

`project_counts_to_json` writes the file and returns its path.
`project_counts` (same arguments, without `out`) returns the same result as
a plain Python dict instead, if you want to inspect or post-process it
before writing anything to disk.

## Getting the JSON into the explorer

Whichever way you ran it, you now have a JSON file (`ibdex_projection.json`
in the examples above). Open the [explorer](https://ibdexmd.org), click
"Upload external JSON", and pick that file. Each sample in it carries the
`sc` (calibrated pathway score) vector, `u1`/`u2` plot coordinates, and a
`cluster` / `cluster_color` / `cluster_confidence` endotype assignment, all
in the schema the explorer expects, so it renders immediately alongside the
reference cohort. Coordinates are produced by locating each new sample's
nearest neighbours in pathway-score space among the 3,168-sample reference
cohort and taking a distance-weighted average of their reference UMAP
positions, since the original UMAP model used to lay out the reference
cohort itself was not retained. `latent` (raw 16-D CVAE coordinates) and the
full `pathway_scores` dictionary are also included in the JSON for
downstream analysis outside the explorer.

## Repository structure

```
IBDEX_Package/
├── src/ibdex_projector/    the installable package
│   ├── artifacts/          frozen model weights, gene panel, reference embedding
│   ├── project.py          end-to-end projection pipeline
│   ├── embedding.py        placement onto the reference UMAP layout
│   ├── model.py            the CVAE architecture and loader
│   ├── pathways.py         decoder-native pathway scoring
│   ├── gene_mapping.py      gene identifier resolution to HGNC symbols
│   └── cli.py              the `ibdex` command line entry point
├── explorer/               the interactive HTML explorer
├── examples/               a small smoke-test dataset
└── tests/
```

## Important interpretation

External samples are projected into the frozen reference without retraining.
All dataset condition indicators are set to zero. This matches the public
projection setting, where the incoming cohort is not one of the four
training datasets.

Pathway scores are latent decoder-derived pathway coordinates, standardized
against the internal IBDEX reference distribution. They are not ssGSEA scores.

## Citation

A manuscript describing IBDEX is in preparation. This section will be
updated with the full citation once it is available.

## License

MIT. See [LICENSE](LICENSE).
