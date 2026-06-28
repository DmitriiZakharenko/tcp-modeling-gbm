"""
Build a compact modeling dataset for team members (no raw NIfTI required).

Merges cohort.csv + features.csv for included patients and writes
data/processed/modeling_table.csv.

Usage
-----
    python -m src.data.export_modeling_dataset
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.config import COHORT_CSV, DATA_PROCESSED, FEATURES_CSV
from src.data.dvh_qc import apply_dvh_qc_to_cohort

MODELING_TABLE_CSV = DATA_PROCESSED / "modeling_table.csv"


def export_modeling_dataset(
    cohort_path: Path = COHORT_CSV,
    features_path: Path = FEATURES_CSV,
    output_path: Path = MODELING_TABLE_CSV,
    apply_qc: bool = True,
) -> pd.DataFrame:
    """Merge included cohort rows with DVH features for TCP modeling."""
    if not cohort_path.exists():
        raise FileNotFoundError(f"Cohort not found: {cohort_path}")
    if not features_path.exists():
        raise FileNotFoundError(
            f"Features not found: {features_path}\n"
            "Run: python -m src.data.feature_builder --workers 4"
        )

    cohort = pd.read_csv(cohort_path)
    features = pd.read_csv(features_path)
    cohort["patient_id"] = cohort["patient_id"].astype(str)
    features["patient_id"] = features["patient_id"].astype(str)

    if apply_qc:
        cohort, excluded = apply_dvh_qc_to_cohort(cohort, features)
        if excluded:
            cohort.to_csv(cohort_path, index=False)
            print(f"Updated {cohort_path.name}: excluded {len(excluded)} patient(s) after DVH QC")

    included = cohort[cohort["included"] == True].copy()  # noqa: E712
    table = included.merge(features, on="patient_id", how="inner", validate="one_to_one")

    missing_features = set(included["patient_id"]) - set(table["patient_id"])
    if missing_features:
        raise RuntimeError(
            f"{len(missing_features)} included patient(s) missing from features.csv: "
            f"{sorted(missing_features)[:10]}"
        )

    table = table.sort_values("patient_id").reset_index(drop=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)

    print(f"Modeling table: {output_path} ({len(table)} patients, {len(table.columns)} columns)")
    return table


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export merged cohort + DVH features for modeling")
    parser.add_argument("--cohort", type=Path, default=COHORT_CSV)
    parser.add_argument("--features", type=Path, default=FEATURES_CSV)
    parser.add_argument("--output", type=Path, default=MODELING_TABLE_CSV)
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    export_modeling_dataset(
        cohort_path=args.cohort,
        features_path=args.features,
        output_path=args.output,
    )
