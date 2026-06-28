"""
RANO outcome helpers for CFB-GBM v3 supplementary data.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from src.config import RANO_TSV


def is_rano_pd(label: object) -> bool:
    """Return True if RANO label indicates progressive disease."""
    if label is None or (isinstance(label, float) and np.isnan(label)):
        return False
    return "Progressive Disease" in str(label)


def rano_controlled(label: object) -> float:
    """
    Binary tumor-control proxy from RANO label.

    1 = non-PD (SD/MR/PR/CR), 0 = PD, NaN = missing.
    """
    if label is None or (isinstance(label, float) and np.isnan(label)):
        return np.nan
    text = str(label).strip()
    if not text:
        return np.nan
    return 0.0 if is_rano_pd(text) else 1.0


def load_rano(path=RANO_TSV) -> pd.DataFrame:
    """
    Load RANO criteria TSV (v3).

    Returns
    -------
    pd.DataFrame
        One row per patient with snake_case RANO columns.
    """
    df = pd.read_csv(path, sep="\t")
    rename = {
        "id_patient": "patient_id",
        "size_t0 (cm3)": "size_t0_cm3",
        "size_t1 (cm3)": "size_t1_cm3",
        "size_t2 (cm3)": "size_t2_cm3",
        "reduction_rate_t0_to_t1": "reduction_rate_t0_t1",
        "rano_t0_to_t1": "rano_t0_t1",
        "reduction_rate_t0_to_t2": "reduction_rate_t0_t2",
        "rano_t0_to_t2": "rano_t0_t2",
        "reduction_rate_t1_to_t2": "reduction_rate_t1_t2",
        "rano_t1_to_t2": "rano_t1_t2",
    }
    df = df.rename(columns=rename)
    for col in (
        "size_t0_cm3",
        "size_t1_cm3",
        "size_t2_cm3",
        "reduction_rate_t0_t1",
        "reduction_rate_t0_t2",
        "reduction_rate_t1_t2",
    ):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    for rano_col, out_col in (
        ("rano_t0_t1", "rano_controlled_t1"),
        ("rano_t0_t2", "rano_controlled_t2"),
        ("rano_t1_t2", "rano_controlled_t1_t2"),
    ):
        if rano_col in df.columns:
            df[out_col] = df[rano_col].map(rano_controlled)

    return df


def rano_summary(frame: pd.DataFrame, rano_col: str = "rano_t0_t1") -> pd.DataFrame:
    """Count RANO categories among non-null patients."""
    if rano_col not in frame.columns:
        return pd.DataFrame()
    counts = frame[rano_col].value_counts(dropna=True)
    return counts.rename_axis("category").reset_index(name="n")
