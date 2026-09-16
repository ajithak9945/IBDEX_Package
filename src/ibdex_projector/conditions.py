from __future__ import annotations

import numpy as np
import pandas as pd

from .artifacts import artifact_path


TISSUE_ALIASES = {
    "colon": "colon_unspecified",
    "large_bowel": "colon_unspecified",
    "colonic": "colon_unspecified",
    "ileum": "ileum",
    "ileal": "ileum",
    "terminal_ileum": "ileum",
    "left": "left_colon",
    "left_colon": "left_colon",
    "right": "right_colon",
    "right_colon": "right_colon",
    "rectum": "rectum",
    "rectal": "rectum",
    "sigmoid": "sigmoid",
    "sigmoid_colon": "sigmoid",
    "transverse": "transverse_colon",
    "transverse_colon": "transverse_colon",
    "colon_unspecified": "colon_unspecified",
}


def load_condition_columns() -> tuple[list[str], list[str]]:
    enc = pd.read_csv(artifact_path("encoder_condition_columns.csv")).iloc[:, 0].astype(str).tolist()
    dec = pd.read_csv(artifact_path("decoder_condition_columns.csv")).iloc[:, 0].astype(str).tolist()
    return enc, dec


def normalize_tissue(value: object) -> str:
    key = str(value).strip().lower().replace(" ", "_").replace("-", "_")
    if key not in TISSUE_ALIASES:
        allowed = ", ".join(sorted(TISSUE_ALIASES))
        raise ValueError(f"unknown tissue value {value!r}; accepted values include {allowed}")
    return TISSUE_ALIASES[key]


def build_conditions(tissues: list[object]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    enc_cols, dec_cols = load_condition_columns()
    enc = np.zeros((len(tissues), len(enc_cols)), dtype="float32")
    dec = np.zeros((len(tissues), len(dec_cols)), dtype="float32")
    normalized = [normalize_tissue(t) for t in tissues]

    for i, tissue in enumerate(normalized):
        broad = "ileum" if tissue == "ileum" else "colon"
        enc_name = f"broad_tissue_model_{broad}"
        if enc_name in enc_cols:
            enc[i, enc_cols.index(enc_name)] = 1.0
        dec_name = f"tissue_raw_model_{tissue}"
        if dec_name in dec_cols:
            dec[i, dec_cols.index(dec_name)] = 1.0
    return enc, dec, normalized
