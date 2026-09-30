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
- Worked example (raw counts &rarr; JSON, real GSE66207 run) at [ibdexmd.org/usage.html](https://ibdexmd.org/usage.html)
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
pandas >= 1.5, torch >= 2.0, scikit-learn >= 1.3, scipy >= 1.10, and
umap-learn >= 0.5 (the last three power endotype/UMAP placement, see
"Important interpretation" below). A GPU is not required; the frozen model
runs on CPU (`--device cpu`, the default).

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

## Step-by-step walkthrough

`project_counts_to_json` above is a convenience wrapper. Internally it calls
six functions in sequence, and each one can be called on its own, which is
useful for debugging or for building a custom pipeline. The output below is
real, captured by running each step against `examples/counts_smoke.csv`.

**1. Load the raw counts file.**

```python
from ibdex_projector.preprocess import load_counts

counts = load_counts("examples/counts_smoke.csv")
print(counts.shape)      # (3000, 3)
print(counts.columns.tolist())   # ['gene_id', 'sample_001', 'sample_002']
```

**2. Resolve gene identifiers to HGNC symbols, normalize, and restrict to the 3,000-gene panel.**

```python
from ibdex_projector.preprocess import preprocess_counts

x, qc = preprocess_counts(counts)
print(x.shape)   # (2, 3000), samples x panel genes
print(qc)
# {'total_input_genes': 3000, 'mapped_genes': 3000, 'unmapped_genes': 0,
#  'duplicated_symbols_after_mapping': 0, 'unique_symbols_after_collapse': 3000,
#  'orientation': 'genes_x_samples', 'n_samples': 2, 'n_panel_genes': 3000,
#  'n_panel_genes_present': 3000, 'n_panel_genes_missing': 0, 'missing_panel_genes': []}
```

**3. Load metadata and build the tissue-conditioning vectors the encoder expects.**

```python
import pandas as pd
from ibdex_projector.conditions import build_conditions

metadata = pd.read_csv("examples/metadata.csv").set_index("sample_id")
enc_cond, dec_cond, tissues = build_conditions(metadata["tissue"].tolist())
print(metadata["tissue"].tolist())   # ['colon', 'ileum']
print(tissues)                       # ['colon_unspecified', 'ileum']
print(enc_cond.shape)                # (2, 6)
```

`colon` resolves to `colon_unspecified` here. That is expected, not an
error, the model conditions on specific colon sub-sites (ascending,
descending, sigmoid, and so on) and falls back to an unspecified colon
category when the metadata does not name one.

**4. Load the frozen model and encode into the 16-dimensional latent space.**

```python
import torch
from ibdex_projector.model import load_model

model, config = load_model(device="cpu")
print(config["best_params"]["latent_dim"])   # 16

with torch.no_grad():
    tx = torch.tensor(x.to_numpy(dtype="float32"))
    tc = torch.tensor(enc_cond)
    mu, logvar = model.encode(tx, tc)
z = mu.numpy()
print(z.shape)          # (2, 16)
print(z[0][:4])         # [1.3847817, -2.0976980, 1.4880019, 0.2173322]
```

**5. Score pathways from the latent codes.**

```python
from ibdex_projector.pathways import score_pathways

pathway_scores = score_pathways(z, tissues)
print(pathway_scores.shape)   # (2, 65), samples x scored pathways
print(pathway_scores.iloc[0][["HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                               "HALLMARK_INFLAMMATORY_RESPONSE",
                               "HALLMARK_IL6_JAK_STAT3_SIGNALING"]].to_dict())
# {'HALLMARK_TNFA_SIGNALING_VIA_NFKB': 2.185941, 'HALLMARK_INFLAMMATORY_RESPONSE': 0.478269,
#  'HALLMARK_IL6_JAK_STAT3_SIGNALING': 0.720717}
```

Each pathway's raw score is a dot product of the sample's latent code `z`
against that pathway's direction vector (`pathway_direction_vectors.csv`),
then standardized `(raw - mean) / std` against the internal reference
cohort's own distribution of that same raw score for that tissue
(`pathway_score_calibration.csv`, `n=2198` colon / `n=970` ileum). Sixty-five
pathways get scored in total, but only the sixteen the explorer displays get
carried into the output JSON's `sc` vector.

**6. Build the ordered `sc` vector and place the sample onto the reference UMAP layout.**

```python
from ibdex_projector.embedding import sc_vector, place_sample

row0 = pathway_scores.iloc[0].to_dict()
sc = sc_vector(row0)
placement = place_sample(z[0], tissues[0])
print(sc)
# [0.478, 2.186, 0.721, 0.019, 1.27, 1.376, 2.568, 0.964,
#  1.902, 1.765, 2.548, 2.449, 0.339, 0.051, 2.673, 1.266]
print(placement)
# {'u1': 11.040817, 'u2': 1.026376, 'cluster': 'C3',
#  'cluster_color': '#6ACC65', 'cluster_confidence': 0.031}
```

`place_sample` reproduces the original scientific pipeline's placement
directly, from the sample's raw 16-D latent code `z` (not from `sc`).
Cluster assignment is nearest-centroid in standardized latent space: a
`StandardScaler` fit on the tissue-specific reference cohort's own latent
codes (`reference_latent.npz`), centroids computed as the mean standardized
position of each pre-discovered endotype (C1-C4 colon, I1-I4 ileum), and the
new sample assigned to whichever centroid is closest by Euclidean distance.
`u1`/`u2` come from a UMAP reducer (`n_neighbors=15, min_dist=0.3,
metric="euclidean", random_state=42`) fit on the full 3,168-sample reference
cohort's standardized latent codes, with the new sample placed via that
fitted reducer's own `.transform()`. `cluster_confidence` is the relative
margin between the nearest and second-nearest centroid (closer to 1.0 means
clearly inside one cluster; closer to 0 means near a cluster boundary) --
note the low value above is expected here too, since `counts_smoke.csv` is
synthetic placeholder data with no real biological cluster structure, not a
real sample.

That `sc` vector, together with `placement`, is exactly what
`project_counts_to_json` assembles into each sample's entry in the output
JSON.

## Getting the JSON into the explorer

Whichever way you ran it, you now have a JSON file (`ibdex_projection.json`
in the examples above). Open the [explorer](https://ibdexmd.org), click
"Upload external JSON", and pick that file. Each sample in it carries the
`sc` (calibrated pathway score) vector, `u1`/`u2` plot coordinates, and a
`cluster` / `cluster_color` / `cluster_confidence` endotype assignment, all
in the schema the explorer expects, so it renders immediately alongside the
reference cohort. Coordinates and endotype assignment reproduce the original scientific
pipeline directly: cluster assignment is nearest-centroid in standardized
16-D CVAE latent space (against the tissue-specific reference cohort's
discovered C1-C4 / I1-I4 centroids), and `u1`/`u2` come from a UMAP reducer
fit on the full 3,168-sample reference cohort's latent codes, transformed
for each new sample. `latent` (raw 16-D CVAE coordinates) and the full
`pathway_scores` dictionary are also included in the JSON for downstream
analysis outside the explorer.

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
├── examples/               example inputs and outputs
│   ├── counts_smoke.csv, metadata.csv   a small 2-sample synthetic smoke test
│   └── gse66207/           a real 33-sample public cohort run end to end,
│                           see the worked example linked above
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
