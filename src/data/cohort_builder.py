"""
Cohort builder: merge CFB-GBM clinical and treatment TSV files, apply
inclusion/exclusion criteria, and export a clean cohort table.

CFB-GBM v3 (2026-06) adds RANO criteria, MRI/CT availability, and extended
clinical / imaging metadata. See ``src/config.py`` for file paths.

Inclusion criteria
------------------
1. RTDOSE available at t0 (rtdose == 1 in imaging availability TSV, temporality == t0).
2. GTV segmentation available at t0 (gtv == 1, temporality == t0).
3. Radiotherapy dose is known (dose not NaN).
4. Number of fractions is known (fractions_number not NaN).

Output
------
data/processed/cohort.csv
"""

from __future__ import annotations

import pandas as pd
import numpy as np

from src.config import (
    CLINICAL_TSV,
    TREATMENT_TSV,
    IMAGING_AVAILABILITY_TSV,
    MRI_AVAILABILITY_TSV,
    CT_AVAILABILITY_TSV,
    RANO_TSV,
    DATA_PROCESSED,
    COHORT_CSV,
    ALPHA_BETA_GBM,
    CLINICAL_FILES,
)
from src.data.rano_loader import load_rano


def require_clinical_files() -> None:
    """
    Raise FileNotFoundError if any required clinical TSV is missing.

    The error message lists missing files and their download URLs.
    """
    missing = [f for f in CLINICAL_FILES if f.required and not f.path.exists()]
    if not missing:
        return

    lines = ["Missing clinical TSV file(s). Download with:"]
    lines.append("  python -m src.data.download_clinical_data")
    lines.append("")
    for clinical_file in missing:
        lines.append(f"  {clinical_file.path.name}")
        lines.append(f"    {clinical_file.url}")
    raise FileNotFoundError("\n".join(lines))


def _build_exclusion_reasons(cohort: pd.DataFrame) -> pd.Series:
    """Vectorized exclusion reason strings for each patient."""
    reasons = pd.Series("", index=cohort.index, dtype="object")

    if "has_rtdose" in cohort.columns:
        mask = ~cohort["has_rtdose"].fillna(False)
        reasons = reasons.mask(mask, reasons + "missing RTDOSE; ")

    if "has_gtv" in cohort.columns:
        mask = ~cohort["has_gtv"].fillna(False)
        reasons = reasons.mask(mask, reasons + "missing GTV; ")

    if "rt_dose_gy" in cohort.columns:
        mask = cohort["rt_dose_gy"].isna()
        reasons = reasons.mask(mask, reasons + "unknown RT dose; ")

    if "n_fractions" in cohort.columns:
        mask = cohort["n_fractions"].isna()
        reasons = reasons.mask(mask, reasons + "unknown n_fractions; ")

    return reasons.str.rstrip("; ").fillna("")


def compute_eqd2(total_dose_gy: float, n_fractions: float, alpha_beta: float = ALPHA_BETA_GBM) -> float:
    """
    Compute EQD2 (Equivalent Dose in 2 Gy fractions) using the LQ model.

    EQD2 = D_total * (d_fraction + alpha_beta) / (2 + alpha_beta)
    """
    if pd.isna(total_dose_gy) or pd.isna(n_fractions) or n_fractions <= 0:
        return np.nan
    d_fraction = total_dose_gy / n_fractions
    return total_dose_gy * (d_fraction + alpha_beta) / (2.0 + alpha_beta)


def load_clinical(path=CLINICAL_TSV) -> pd.DataFrame:
    """Load clinical TSV (v3: adds height, weight, who_guideline)."""
    df = pd.read_csv(path, sep="\t")
    df = df.rename(
        columns={
            "id_patient": "patient_id",
            "survival (weeks)": "survival_weeks",
            "age_at_t0 (years)": "age",
            "who_performance_status": "who_status",
            "gender": "sex",
            "height (cm)": "height_cm",
            "weight (kg)": "weight_kg",
        }
    )
    cols = [
        c
        for c in [
            "patient_id",
            "survival_weeks",
            "age",
            "sex",
            "who_status",
            "who_guideline",
            "height_cm",
            "weight_kg",
        ]
        if c in df.columns
    ]
    df = df[cols].copy()
    df["survival_weeks"] = pd.to_numeric(df["survival_weeks"], errors="coerce")
    df["age"] = pd.to_numeric(df["age"], errors="coerce")
    df["who_status"] = pd.to_numeric(df["who_status"], errors="coerce")
    for col in ("height_cm", "weight_kg"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if "height_cm" in df.columns and "weight_kg" in df.columns:
        df["bmi"] = df["weight_kg"] / (df["height_cm"] / 100.0) ** 2
    return df


def load_treatment(path=TREATMENT_TSV) -> pd.DataFrame:
    """Load treatment TSV including RT delay."""
    df = pd.read_csv(path, sep="\t")
    df = df.rename(
        columns={
            "id_patient": "patient_id",
            "dose (Gy)": "rt_dose_gy",
            "fractions_number": "n_fractions",
            "delay_t0_to_radiotherapy (weeks)": "rt_delay_wk",
        }
    )
    cols = [c for c in ["patient_id", "rt_dose_gy", "n_fractions", "rt_delay_wk"] if c in df.columns]
    df = df[cols].copy()
    df["rt_dose_gy"] = pd.to_numeric(df["rt_dose_gy"], errors="coerce")
    df["n_fractions"] = pd.to_numeric(df["n_fractions"], errors="coerce")
    if "rt_delay_wk" in df.columns:
        df["rt_delay_wk"] = pd.to_numeric(df["rt_delay_wk"], errors="coerce")
    return df


def load_imaging_availability(path=IMAGING_AVAILABILITY_TSV) -> pd.DataFrame:
    """Load t0 RTDOSE/GTV flags plus machine/TPS metadata (v3)."""
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()

    if "temporality" in df.columns:
        df = df[df["temporality"].str.strip().str.lower() == "t0"].copy()

    df = df.rename(
        columns={
            "id_patient": "patient_id",
            "rtdose": "has_rtdose",
            "gtv": "has_gtv",
            "treatment_machine": "rt_machine",
            "tps": "rt_tps",
            "gtv_type": "gtv_segmentation_type",
        }
    )
    cols = [
        c
        for c in [
            "patient_id",
            "has_rtdose",
            "has_gtv",
            "rt_machine",
            "rt_tps",
            "gtv_segmentation_type",
        ]
        if c in df.columns
    ]
    df = df[cols].copy()

    for col in ("has_rtdose", "has_gtv"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(bool)

    return df


def load_mri_followup(path=MRI_AVAILABILITY_TSV) -> pd.DataFrame:
    """
    Summarise follow-up MRI availability per patient.

    Adds ``mri_t1_weeks``, ``has_mri_t1``, ``has_mri_t2`` from MRI availability TSV.
    """
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()
    df = df.rename(columns={"id_patient": "patient_id", "time_diff_t0 (weeks)": "time_diff_t0_wk"})
    df["temporality"] = df["temporality"].str.strip().str.lower()
    df["time_diff_t0_wk"] = pd.to_numeric(df["time_diff_t0_wk"], errors="coerce")

    rows = []
    for patient_id, sub in df.groupby("patient_id"):
        row = {"patient_id": patient_id}
        for tp in ("t1", "t2"):
            tp_rows = sub[sub["temporality"] == tp]
            row[f"has_mri_{tp}"] = not tp_rows.empty
            if tp == "t1" and not tp_rows.empty:
                row["mri_t1_weeks"] = float(tp_rows["time_diff_t0_wk"].iloc[0])
        rows.append(row)
    return pd.DataFrame(rows)


def load_ct_availability(path=CT_AVAILABILITY_TSV) -> pd.DataFrame:
    """Extract t0 CT availability and scanner type."""
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()
    if "temporality" in df.columns:
        df = df[df["temporality"].str.strip().str.lower() == "t0"].copy()
    df = df.rename(columns={"id_patient": "patient_id", "ct": "has_ct", "ct_machine": "ct_machine"})
    cols = [c for c in ["patient_id", "has_ct", "ct_machine"] if c in df.columns]
    df = df[cols].copy()
    if "has_ct" in df.columns:
        df["has_ct"] = pd.to_numeric(df["has_ct"], errors="coerce").fillna(0).astype(bool)
    return df


def build_cohort(
    clinical_path=CLINICAL_TSV,
    treatment_path=TREATMENT_TSV,
    imaging_path=IMAGING_AVAILABILITY_TSV,
    mri_path=MRI_AVAILABILITY_TSV,
    ct_path=CT_AVAILABILITY_TSV,
    rano_path=RANO_TSV,
    output_path=None,
) -> pd.DataFrame:
    """
    Merge clinical, treatment, imaging, RANO; apply inclusion criteria; export cohort.
    """
    if output_path is None:
        output_path = COHORT_CSV

    require_clinical_files()

    clinical = load_clinical(clinical_path)
    treatment = load_treatment(treatment_path)
    imaging = load_imaging_availability(imaging_path)
    mri_followup = load_mri_followup(mri_path) if mri_path.exists() else pd.DataFrame()
    ct_avail = load_ct_availability(ct_path) if ct_path.exists() else pd.DataFrame()
    rano = load_rano(rano_path) if rano_path.exists() else pd.DataFrame()

    cohort = clinical.merge(treatment, on="patient_id", how="outer")
    cohort = cohort.merge(imaging, on="patient_id", how="outer")
    if not mri_followup.empty:
        cohort = cohort.merge(mri_followup, on="patient_id", how="left")
    if not ct_avail.empty:
        cohort = cohort.merge(ct_avail, on="patient_id", how="left")
    if not rano.empty:
        cohort = cohort.merge(rano, on="patient_id", how="left")

    for col in ("has_rtdose", "has_gtv"):
        if col not in cohort.columns:
            cohort[col] = False
        else:
            cohort[col] = cohort[col].fillna(False)

    cohort["exclusion_reason"] = _build_exclusion_reasons(cohort)
    cohort["included"] = cohort["exclusion_reason"] == ""

    d_fraction = cohort["rt_dose_gy"] / cohort["n_fractions"]
    cohort["eqd2_gy"] = np.where(
        cohort["included"],
        cohort["rt_dose_gy"] * (d_fraction + ALPHA_BETA_GBM) / (2.0 + ALPHA_BETA_GBM),
        np.nan,
    )

    col_order = [
        "patient_id",
        "rt_dose_gy",
        "n_fractions",
        "eqd2_gy",
        "rt_delay_wk",
        "survival_weeks",
        "age",
        "sex",
        "who_status",
        "who_guideline",
        "height_cm",
        "weight_kg",
        "bmi",
        "has_rtdose",
        "has_gtv",
        "gtv_segmentation_type",
        "rt_machine",
        "rt_tps",
        "has_ct",
        "ct_machine",
        "has_mri_t1",
        "mri_t1_weeks",
        "has_mri_t2",
        "size_t0_cm3",
        "size_t1_cm3",
        "size_t2_cm3",
        "reduction_rate_t0_t1",
        "rano_t0_t1",
        "rano_controlled_t1",
        "reduction_rate_t0_t2",
        "rano_t0_t2",
        "rano_controlled_t2",
        "reduction_rate_t1_t2",
        "rano_t1_t2",
        "rano_controlled_t1_t2",
        "included",
        "exclusion_reason",
    ]
    col_order = [c for c in col_order if c in cohort.columns]
    cohort = cohort[col_order].sort_values("patient_id").reset_index(drop=True)

    cohort.to_csv(output_path, index=False)

    n_rano = int(cohort["rano_t0_t1"].notna().sum()) if "rano_t0_t1" in cohort.columns else 0
    print(f"Cohort saved to: {output_path}")
    print(f"  Total patients : {len(cohort)}")
    print(f"  Included       : {cohort['included'].sum()}")
    print(f"  Excluded       : {(~cohort['included']).sum()}")
    if n_rano:
        print(f"  With RANO t0→t1: {n_rano}")
    print("\nExclusion reasons:")
    for reason, count in cohort[~cohort["included"]]["exclusion_reason"].value_counts().items():
        print(f"  {count:3d}  {reason}")

    return cohort


if __name__ == "__main__":
    build_cohort()
