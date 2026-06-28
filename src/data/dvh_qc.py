"""
Post-hoc DVH quality control: exclude patients with non-physical dose in GTV.

Usage
-----
    python -m src.data.dvh_qc
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.config import COHORT_CSV, FEATURES_CSV

# Minimum mean GTV dose (Gy) for an RTDOSE/GTV pair to be modeling-eligible.
MIN_DMEAN_GY: float = 1.0

DVH_QC_EXCLUSION_REASON = "invalid RTDOSE (zero dose in GTV)"


def apply_dvh_qc_to_cohort(
    cohort: pd.DataFrame,
    features: pd.DataFrame,
    min_dmean_gy: float = MIN_DMEAN_GY,
) -> tuple[pd.DataFrame, list[str]]:
    """
    Mark included patients as excluded when DVH metrics indicate invalid RTDOSE.

    Returns updated cohort and list of newly excluded patient IDs.
    """
    cohort = cohort.copy()
    cohort["patient_id"] = cohort["patient_id"].astype(str)
    features = features.copy()
    features["patient_id"] = features["patient_id"].astype(str)

    dvh = features[["patient_id", "Dmean_gy"]]
    merged = cohort.merge(dvh, on="patient_id", how="left")

    fail_mask = merged["included"].fillna(False) & (merged["Dmean_gy"].fillna(0.0) < min_dmean_gy)
    failed_ids = merged.loc[fail_mask, "patient_id"].tolist()

    for patient_id in failed_ids:
        idx = cohort.index[cohort["patient_id"] == patient_id]
        if idx.empty:
            continue
        row_idx = idx[0]
        cohort.at[row_idx, "included"] = False
        existing = cohort.at[row_idx, "exclusion_reason"]
        if pd.isna(existing) or str(existing).strip() == "":
            cohort.at[row_idx, "exclusion_reason"] = DVH_QC_EXCLUSION_REASON
        elif DVH_QC_EXCLUSION_REASON not in str(existing):
            cohort.at[row_idx, "exclusion_reason"] = f"{existing}; {DVH_QC_EXCLUSION_REASON}"

    return cohort, failed_ids


def run_dvh_qc(
    cohort_path: Path = COHORT_CSV,
    features_path: Path = FEATURES_CSV,
    min_dmean_gy: float = MIN_DMEAN_GY,
    write_cohort: bool = True,
) -> list[str]:
    """Apply DVH QC and optionally write updated cohort.csv."""
    cohort = pd.read_csv(cohort_path)
    features = pd.read_csv(features_path)
    cohort, failed_ids = apply_dvh_qc_to_cohort(cohort, features, min_dmean_gy=min_dmean_gy)

    if write_cohort:
        cohort.to_csv(cohort_path, index=False)

    if failed_ids:
        print(f"DVH QC excluded {len(failed_ids)} patient(s): {', '.join(failed_ids)}")
    else:
        print("DVH QC: no additional exclusions.")

    return failed_ids


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Apply DVH QC exclusions to cohort.csv")
    parser.add_argument("--cohort", type=Path, default=COHORT_CSV)
    parser.add_argument("--features", type=Path, default=FEATURES_CSV)
    parser.add_argument("--min-dmean", type=float, default=MIN_DMEAN_GY, metavar="GY")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    run_dvh_qc(
        cohort_path=args.cohort,
        features_path=args.features,
        min_dmean_gy=args.min_dmean,
    )
