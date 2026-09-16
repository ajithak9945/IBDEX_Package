from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from .conditions import build_conditions
from .model import load_model
from .pathways import score_pathways
from .preprocess import load_counts, preprocess_counts


def project_counts(
    counts_path: str | Path,
    metadata_path: str | Path,
    sample_col: str = "sample_id",
    tissue_col: str = "tissue",
    gene_col: str | None = None,
    device: str = "cpu",
) -> dict:
    counts = load_counts(counts_path)
    metadata = pd.read_csv(metadata_path)
    if sample_col not in metadata.columns:
        raise ValueError(f"metadata sample column not found: {sample_col}")
    if tissue_col not in metadata.columns:
        raise ValueError(f"metadata tissue column not found: {tissue_col}")

    x, qc = preprocess_counts(counts, gene_col=gene_col)
    metadata = metadata.set_index(sample_col)
    missing_meta = [s for s in x.index.astype(str) if s not in metadata.index.astype(str)]
    metadata.index = metadata.index.astype(str)
    x.index = x.index.astype(str)
    if missing_meta:
        raise ValueError(f"metadata is missing {len(missing_meta)} samples, including {missing_meta[:5]}")
    metadata = metadata.loc[x.index]

    enc_cond, _dec_cond, tissues = build_conditions(metadata[tissue_col].tolist())
    model, config = load_model(device=device)
    with torch.no_grad():
        tx = torch.tensor(x.to_numpy(dtype="float32"), device=device)
        tc = torch.tensor(enc_cond, device=device)
        mu, logvar = model.encode(tx, tc)
    z = mu.cpu().numpy()
    pathway_scores = score_pathways(z, tissues)

    samples = []
    for i, sample_id in enumerate(x.index.tolist()):
        samples.append({
            "id": sample_id,
            "tissue": tissues[i],
            "latent": {f"z{j}": float(z[i, j]) for j in range(z.shape[1])},
            "pathway_scores": {
                k: (None if np.isnan(v) else float(v))
                for k, v in pathway_scores.iloc[i].to_dict().items()
            },
        })

    return {
        "ibdex_version": "0.1.0",
        "model": {
            "name": "IBDEX frozen tissue-conditioned CVAE",
            "latent_dim": int(config["best_params"]["latent_dim"]),
            "dataset_condition": "all_zero_external_projection",
        },
        "input_qc": qc,
        "samples": samples,
    }


def project_counts_to_json(*args, out: str | Path, **kwargs) -> Path:
    result = project_counts(*args, **kwargs)
    out = Path(out)
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return out
