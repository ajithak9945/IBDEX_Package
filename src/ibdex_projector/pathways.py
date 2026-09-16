from __future__ import annotations

import numpy as np
import pandas as pd

from .artifacts import artifact_path


def score_pathways(z: np.ndarray, tissues: list[str], variant: str = "all_nes") -> pd.DataFrame:
    vectors = pd.read_csv(artifact_path("pathway_direction_vectors.csv"))
    cal = pd.read_csv(artifact_path("pathway_score_calibration.csv"))
    dim_cols = [c for c in vectors.columns if c.startswith("dim_")]
    out = pd.DataFrame(index=range(z.shape[0]))
    for tissue in sorted(set(_broad_tissue(t) for t in tissues)):
        sample_idx = [i for i, t in enumerate(tissues) if _broad_tissue(t) == tissue]
        sub = vectors[(vectors["tissue"].astype(str) == tissue) & (vectors["weighting_variant"] == variant)]
        if sub.empty:
            continue
        weights = sub[dim_cols].to_numpy(dtype="float32")
        raw_scores = z[sample_idx] @ weights.T
        for j, pathway in enumerate(sub["pathway"].astype(str).tolist()):
            c = cal[(cal["tissue"].astype(str) == tissue) & (cal["pathway"].astype(str) == pathway)]
            mean = float(c["mean"].iloc[0]) if not c.empty else 0.0
            std = float(c["std"].iloc[0]) if not c.empty and float(c["std"].iloc[0]) != 0 else 1.0
            out.loc[sample_idx, pathway] = (raw_scores[:, j] - mean) / std
    return out.astype(float)


def _broad_tissue(tissue: str) -> str:
    return "ileum" if str(tissue).lower() == "ileum" else "colon"
