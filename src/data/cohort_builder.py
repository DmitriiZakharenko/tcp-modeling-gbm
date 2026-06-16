"""
Cohort builder: merge CFB-GBM clinical and treatment TSV files, apply
inclusion/exclusion criteria, and export a clean cohort table.

Actual TSV column names (v02, 2026-01-29)
------------------------------------------
clinical_data    : id_patient, survival (weeks), age_at_t0 (years),
                   who_performance_status, gender
treatment_data   : id_patient, delay_t0_to_radiotherapy (weeks),
                   dose (Gy), fractions_number
imaging_avail    : id_patient, temporality, gtv, rtdose, treatment_machine, tps

Inclusion criteria
------------------
1. RTDOSE available at t0 (rtdose == 1 in imaging availability TSV, temporality == t0).
2. GTV segmentation available at t0 (gtv == 1, temporality == t0).
3. Radiotherapy dose is known (dose not NaN).
4. Number of fractions is known (fractions_number not NaN).

Output
------
data/processed/cohort.csv with columns:
    patient_id, rt_dose_gy, n_fractions, eqd2_gy, survival_weeks,
    age, sex, who_status, has_rtdose, has_gtv, included, exclusion_reason
"""

import pandas as pd
import numpy as np

from src.config import (
    CLINICAL_TSV,
    TREATMENT_TSV,
    IMAGING_AVAILABILITY_TSV,
    DATA_PROCESSED,
    ALPHA_BETA_GBM,
)


def compute_eqd2(total_dose_gy: float, n_fractions: float, alpha_beta: float = ALPHA_BETA_GBM) -> float:
    """
    Compute EQD2 (Equivalent Dose in 2 Gy fractions) using the LQ model.

    EQD2 = D_total * (d_fraction + alpha_beta) / (2 + alpha_beta)

    Parameters
    ----------
    total_dose_gy : float
        Total prescribed dose in Gy.
    n_fractions : float
        Number of fractions.
    alpha_beta : float
        Alpha/beta ratio in Gy (default 10.0 for GBM).

    Returns
    -------
    float
        EQD2 in Gy, or NaN if inputs are invalid.
    """
    if pd.isna(total_dose_gy) or pd.isna(n_fractions) or n_fractions <= 0:
        return np.nan
    d_fraction = total_dose_gy / n_fractions
    return total_dose_gy * (d_fraction + alpha_beta) / (2.0 + alpha_beta)


def load_clinical(path=CLINICAL_TSV) -> pd.DataFrame:
    """
    Load and clean the clinical data TSV.

    Parameters
    ----------
    path : str or Path
        Path to CFB-GBM_clinical_data TSV.

    Returns
    -------
    pd.DataFrame
        Columns: patient_id, survival_weeks, age, sex, who_status.
    """
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={
        "id_patient": "patient_id",
        "survival (weeks)": "survival_weeks",
        "age_at_t0 (years)": "age",
        "who_performance_status": "who_status",
        "gender": "sex",
    })
    cols = [c for c in ["patient_id", "survival_weeks", "age", "sex", "who_status"] if c in df.columns]
    df = df[cols].copy()
    df["survival_weeks"] = pd.to_numeric(df["survival_weeks"], errors="coerce")
    df["age"] = pd.to_numeric(df["age"], errors="coerce")
    df["who_status"] = pd.to_numeric(df["who_status"], errors="coerce")
    return df


def load_treatment(path=TREATMENT_TSV) -> pd.DataFrame:
    """
    Load and clean the treatment data TSV.

    Parameters
    ----------
    path : str or Path
        Path to CFB-GBM_treatment_data TSV.

    Returns
    -------
    pd.DataFrame
        Columns: patient_id, rt_dose_gy, n_fractions.
    """
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={
        "id_patient": "patient_id",
        "dose (Gy)": "rt_dose_gy",
        "fractions_number": "n_fractions",
    })
    cols = [c for c in ["patient_id", "rt_dose_gy", "n_fractions"] if c in df.columns]
    df = df[cols].copy()
    df["rt_dose_gy"] = pd.to_numeric(df["rt_dose_gy"], errors="coerce")
    df["n_fractions"] = pd.to_numeric(df["n_fractions"], errors="coerce")
    return df


def load_imaging_availability(path=IMAGING_AVAILABILITY_TSV) -> pd.DataFrame:
    """
    Load imaging availability TSV and extract RTDOSE and GTV flags at t0.

    Filters to temporality == 't0' before extracting flags.

    Parameters
    ----------
    path : str or Path
        Path to CFB-GBM_treatment_imaging_availability TSV.

    Returns
    -------
    pd.DataFrame
        Columns: patient_id, has_rtdose, has_gtv. One row per patient.
    """
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()

    # Keep only t0 rows
    if "temporality" in df.columns:
        df = df[df["temporality"].str.strip().str.lower() == "t0"].copy()

    df = df.rename(columns={"id_patient": "patient_id", "rtdose": "has_rtdose", "gtv": "has_gtv"})
    cols = [c for c in ["patient_id", "has_rtdose", "has_gtv"] if c in df.columns]
    df = df[cols].copy()

    for col in ("has_rtdose", "has_gtv"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(bool)

    return df


def build_cohort(
    clinical_path=CLINICAL_TSV,
    treatment_path=TREATMENT_TSV,
    imaging_path=IMAGING_AVAILABILITY_TSV,
    output_path=None,
) -> pd.DataFrame:
    """
    Merge clinical, treatment, and imaging data; apply inclusion criteria;
    add EQD2; export cohort table.

    Parameters
    ----------
    clinical_path, treatment_path, imaging_path : str or Path
        Paths to the respective TSV files.
    output_path : str or Path, optional
        Export path for cohort CSV. Defaults to data/processed/cohort.csv.

    Returns
    -------
    pd.DataFrame
        Full cohort including excluded patients.
        `included` column (bool) flags eligible patients.
        `exclusion_reason` is empty string for included patients.
    """
    if output_path is None:
        output_path = DATA_PROCESSED / "cohort.csv"

    clinical = load_clinical(clinical_path)
    treatment = load_treatment(treatment_path)
    imaging = load_imaging_availability(imaging_path)

    cohort = clinical.merge(treatment, on="patient_id", how="outer")
    cohort = cohort.merge(imaging, on="patient_id", how="outer")

    for col in ("has_rtdose", "has_gtv"):
        if col not in cohort.columns:
            cohort[col] = False
        else:
            cohort[col] = cohort[col].fillna(False)

    reasons = []
    for _, row in cohort.iterrows():
        r = []
        if not row.get("has_rtdose", False):
            r.append("missing RTDOSE")
        if not row.get("has_gtv", False):
            r.append("missing GTV")
        if pd.isna(row.get("rt_dose_gy")):
            r.append("unknown RT dose")
        if pd.isna(row.get("n_fractions")):
            r.append("unknown n_fractions")
        reasons.append("; ".join(r))

    cohort["exclusion_reason"] = reasons
    cohort["included"] = cohort["exclusion_reason"] == ""

    cohort["eqd2_gy"] = cohort.apply(
        lambda row: compute_eqd2(row["rt_dose_gy"], row["n_fractions"]) if row["included"] else np.nan,
        axis=1,
    )

    col_order = [
        "patient_id", "rt_dose_gy", "n_fractions", "eqd2_gy",
        "survival_weeks", "age", "sex", "who_status",
        "has_rtdose", "has_gtv", "included", "exclusion_reason",
    ]
    col_order = [c for c in col_order if c in cohort.columns]
    cohort = cohort[col_order].sort_values("patient_id").reset_index(drop=True)

    cohort.to_csv(output_path, index=False)

    print(f"Cohort saved to: {output_path}")
    print(f"  Total patients : {len(cohort)}")
    print(f"  Included       : {cohort['included'].sum()}")
    print(f"  Excluded       : {(~cohort['included']).sum()}")
    print("\nExclusion reasons:")
    for reason, count in cohort[~cohort["included"]]["exclusion_reason"].value_counts().items():
        print(f"  {count:3d}  {reason}")

    return cohort


if __name__ == "__main__":
    build_cohort()
