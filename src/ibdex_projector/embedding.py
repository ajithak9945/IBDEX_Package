"""Place newly-projected samples onto the reference UMAP layout used by the
IBDEX explorer, and assign each one to its nearest discovered endotype.

The explorer's 2D layout (u1, u2 per reference sample) was produced once,
offline, by fitting UMAP on the full training cohort's latent space. That
fitted UMAP reducer was not saved, so new external samples cannot be passed
through its `.transform()`. Instead we reuse the reference cohort's own
calibrated pathway-score vectors (which we do have, bundled as
`reference_embedding.json.gz`) as anchors, and place each new sample at the
distance-weighted average position of its k nearest reference samples in
that 15-dimensional pathway-score space. This keeps new samples visually
close to reference samples with a similar molecular profile, which is the
property the explorer's scatter plot relies on, without needing the
original UMAP model object.
"""
from __future__ import annotations

import gzip
import json
from functools import lru_cache

import numpy as np

from .artifacts import artifact_path

# Pathway ids, in the exact order the explorer's radar/scatter axes (`AX`)
# expect the `sc` array to be built in. Keep this in sync with the `AX`
# array in explorer/IBDEX_refined.html if that file's axis list ever changes.
AX_ORDER = [
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
    "REACTOME_INTERLEUKIN_23_SIGNALING",
    "REACTOME_NOD1_NOD2_SIGNALING_PATHWAY",
    "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
    "KEGG_BILE_SECRETION",
    "HALLMARK_BILE_ACID_METABOLISM",
    "REACTOME_TERMINATION_OF_O_GLYCAN_BIOSYNTHESIS",
]


@lru_cache(maxsize=1)
def _load_reference() -> dict:
    with gzip.open(artifact_path("reference_embedding.json.gz"), "rt", encoding="utf-8") as fh:
        ref = json.load(fh)
    out = {}
    for tissue, rows in ref.items():
        sc = np.array([r["sc"] for r in rows], dtype="float64")
        u1 = np.array([r["u1"] for r in rows], dtype="float64")
        u2 = np.array([r["u2"] for r in rows], dtype="float64")
        et = [r["et"] for r in rows]
        ec = {r["et"]: r["ec"] for r in rows}
        out[tissue] = {"sc": sc, "u1": u1, "u2": u2, "et": et, "ec": ec}
    return out


def broad_tissue(tissue: str) -> str:
    return "ileum" if str(tissue).lower() == "ileum" else "colon"


def sc_vector(pathway_scores: dict) -> list[float]:
    """Build the `sc` array in explorer AX order from a pathway_scores dict."""
    return [float(pathway_scores.get(p, 0.0) or 0.0) for p in AX_ORDER]


def place_sample(sc: list[float], tissue: str, k: int = 15) -> dict:
    """Return {u1, u2, cluster, cluster_color, cluster_confidence} for one sample."""
    ref = _load_reference()[broad_tissue(tissue)]
    diffs = ref["sc"] - np.asarray(sc, dtype="float64")
    dists = np.sqrt((diffs ** 2).sum(axis=1))
    k = min(k, len(dists))
    nn_idx = np.argsort(dists)[:k]
    nn_dists = dists[nn_idx]
    weights = 1.0 / (nn_dists + 1e-6)
    weights = weights / weights.sum()

    u1 = float((weights * ref["u1"][nn_idx]).sum())
    u2 = float((weights * ref["u2"][nn_idx]).sum())

    vote: dict[str, float] = {}
    for w, i in zip(weights, nn_idx):
        cluster = ref["et"][i]
        vote[cluster] = vote.get(cluster, 0.0) + float(w)
    cluster = max(vote, key=vote.get)
    confidence = round(vote[cluster], 4)

    return {
        "u1": u1,
        "u2": u2,
        "cluster": cluster,
        "cluster_color": ref["ec"].get(cluster, "#FFD700"),
        "cluster_confidence": confidence,
    }
