"""Place newly-projected samples onto the reference UMAP layout used by the
IBDEX explorer, and assign each one to its nearest discovered endotype.

This reproduces the methodology used in the original scientific pipeline
(ExternalCVAEValidation.ipynb) exactly, rather than approximating it:

  Cluster assignment: nearest-centroid in *standardized* 16-D CVAE latent
  space. A StandardScaler is fit on the tissue-specific reference cohort's
  raw latent codes (Z_colon / Z_ileum), each pre-discovered endotype's
  centroid is the mean of its standardized reference samples, and a new
  sample is assigned to whichever centroid is closest by Euclidean
  distance after being scaled with that same fitted scaler.

  UMAP placement: a StandardScaler is fit on the full reference cohort's
  raw latent codes (Z_all, colon + ileum together), UMAP is fit on the
  scaled result with the same hyperparameters and random seed used when
  the reference layout was built (n_neighbors=15, min_dist=0.3,
  metric="euclidean", random_state=42), and a new sample is placed with
  that fitted reducer's own `.transform()` after being scaled the same
  way. UMAP's exact pixel coordinates can vary slightly across library
  versions/environments even with a fixed seed, but the reducer is fit
  fresh, once, from the same reference data and hyperparameters the
  original notebook used, so relative placement is faithful.

Earlier versions of this module placed samples using a k-nearest-neighbor
average in *pathway-score* space, as a workaround for not having the raw
reference latent codes bundled. That approximation is no longer needed or
used now that Z_all/Z_colon/Z_ileum and their real endotype labels are
bundled as `reference_latent.npz`.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.spatial.distance import cdist
from sklearn.preprocessing import StandardScaler

from .artifacts import artifact_path

UMAP_N_NEIGHBORS = 15
UMAP_MIN_DIST = 0.3
UMAP_METRIC = "euclidean"
UMAP_SEED = 42

# Pathway ids, in the exact order the explorer's radar axes (`AX`) expect
# the `sc` array to be built in. Kept for the `sc` field in the output
# JSON (informational / radar chart use) -- no longer used for placement.
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


def broad_tissue(tissue: str) -> str:
    return "ileum" if str(tissue).lower() == "ileum" else "colon"


def sc_vector(pathway_scores: dict) -> list[float]:
    """Build the `sc` array in explorer AX order from a pathway_scores dict."""
    return [float(pathway_scores.get(p, 0.0) or 0.0) for p in AX_ORDER]


@lru_cache(maxsize=1)
def _load_latent_reference() -> dict:
    with np.load(artifact_path("reference_latent.npz"), allow_pickle=True) as npz:
        Z = np.asarray(npz["Z"], dtype="float64")
        broad = np.asarray(npz["broad_tissue"])
        cluster_label = np.asarray(npz["cluster_label"])

    out = {}
    for tissue in ("colon", "ileum"):
        mask = broad == tissue
        Z_t = Z[mask]
        labels_t = cluster_label[mask]

        cluster_scaler = StandardScaler().fit(Z_t)
        Z_t_scaled = cluster_scaler.transform(Z_t)
        uniq = sorted(set(labels_t.tolist()))
        centroids = np.stack([Z_t_scaled[labels_t == k].mean(axis=0) for k in uniq])

        out[tissue] = {
            "cluster_scaler": cluster_scaler,
            "centroids": centroids,
            "centroid_labels": np.array(uniq),
        }

    # UMAP is fit once, on the full reference cohort (colon + ileum together),
    # matching the original notebook's Z_all-based fit.
    umap_scaler = StandardScaler().fit(Z)
    Z_umap_input = umap_scaler.transform(Z).astype("float32")

    import umap as umap_lib

    reducer = umap_lib.UMAP(
        n_neighbors=UMAP_N_NEIGHBORS,
        min_dist=UMAP_MIN_DIST,
        metric=UMAP_METRIC,
        random_state=UMAP_SEED,
    )
    reducer.fit(Z_umap_input)

    out["umap_scaler"] = umap_scaler
    out["umap_reducer"] = reducer

    # Cluster colors, matching the explorer's fixed C1-C4 / I1-I4 palette.
    out["cluster_colors"] = {
        "C1": "#4878D0", "C2": "#EE854A", "C3": "#6ACC65", "C4": "#D65F5F",
        "I1": "#4878D0", "I2": "#EE854A", "I3": "#6ACC65", "I4": "#D65F5F",
    }
    return out


def place_sample(z: list[float] | np.ndarray, tissue: str) -> dict:
    """Return {u1, u2, cluster, cluster_color, cluster_confidence} for one sample.

    `z` is the raw 16-D CVAE latent code for the sample (model.encode output),
    not the calibrated pathway-score vector.
    """
    ref = _load_latent_reference()
    t = broad_tissue(tissue)
    z = np.asarray(z, dtype="float64").reshape(1, -1)

    # Cluster: nearest centroid in standardized tissue-specific latent space.
    tissue_ref = ref[t]
    z_cluster = tissue_ref["cluster_scaler"].transform(z)
    dists = cdist(z_cluster, tissue_ref["centroids"], metric="euclidean")[0]
    nearest_idx = int(dists.argmin())
    cluster = str(tissue_ref["centroid_labels"][nearest_idx])

    # Confidence: relative margin between the nearest and second-nearest
    # centroid (1.0 = right on the centroid with a clear runner-up gap,
    # lower = closer to a cluster boundary). Not part of the original
    # notebook's output (which left this field null); provided here as a
    # useful, well-defined addition for the explorer.
    sorted_dists = np.sort(dists)
    d1 = float(sorted_dists[0])
    d2 = float(sorted_dists[1]) if len(sorted_dists) > 1 else d1
    confidence = round(1.0 - d1 / d2, 4) if d2 > 0 else 1.0

    # UMAP: transform through the reducer fit on the full reference cohort.
    z_umap_in = ref["umap_scaler"].transform(z).astype("float32")
    u = ref["umap_reducer"].transform(z_umap_in)[0]

    return {
        "u1": float(u[0]),
        "u2": float(u[1]),
        "cluster": cluster,
        "cluster_color": ref["cluster_colors"].get(cluster, "#FFD700"),
        "cluster_confidence": confidence,
    }
