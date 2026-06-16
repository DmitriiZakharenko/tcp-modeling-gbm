"""
Cohort builder: merge CFB-GBM clinical and treatment TSV files, apply
inclusion/exclusion criteria, and export a clean cohort table.

Inclusion criteria
------------------
1. RTDOSE file available (has_rtdose == True in imaging availability TSV).
2. GTV segmentation available (has_gtv == True in imaging availability TSV).
3. Radiotherapy dose is known (rt_dose_gy is not NaN).
4. Number of fractions is known (n_fractions is not NaN).

Exclusion criteria (documented in cohort table as exclusion_reason)
-------------------------------------------------------------------
- Missing RTDOSE
- Missing GTV segmentation
- Unknown radiotherapy dose
- Unknown number of fractions

Output
------
data/processed/cohort.csv with columns:
    patient_id, rt_dose_gy, n_fractions, eqd2_gy, survival_weeks,
    age, sex, who_status, included, exclusion_reason
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


def compute_eqd2(total_dose_gy: float, n_fractions: int, alpha_beta: float = ALPHA_BETA_GBM) -> float:
    """
    Compute EQD2 (Equivalent Dose in 2 Gy fractions) using the LQ model.

    EQD2 = D_total * (d_fraction + alpha_beta) / (2 + alpha_beta)

    Parameters
    ----------
    total_dose_gy : float
        Total prescribed dose in Gy.
    n_fractions : int
        Number of fractions.
    alpha_beta : float
        Alpha/beta ratio in Gy (default 10.0 for GBM).

    Returns
    -------
    float
        EQD2 in Gy, or NaN if inputs are invalid.
    """
    if n_fractions <= 0 or np.isnan(total_dose_gy) or np.isnan(n_fractions):
        return np.nan
    d_fraction = total_dose_gy / n_fractions
    return total_dose_gy * (d_fraction + alpha_beta) / (2.0 + alpha_beta)


def load_clinical(path: str = CLINICAL_TSV) -> pd.DataFrame:
    """
    Load and minimally clean the clinical data TSV.

    Parameters
    ----------
    path : str or Path
        Path to CFB-GBM_clinical_data TSV file.

    Returns
    -------
    pd.DataFrame
        Columns: patient_id, survival_weeks, age, sex, who_status.

    Raises
    ------
    FileNotFoundError
        If the TSV file does not exist at the given path.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = df.columns.str.strip().str.lower().str.replace(" ", "_")

    rename = {
        "patient_id": "patient_id",
        "overall_survival": "survival_weeks",
        "age": "age",
        "sex_at_birth": "sex",
        "who_status_performance": "who_status",
    }
    # Keep only columns that exist; rename them
    available = {k: v for k, v in rename.items() if k in df.columns}
    df = df[list(available.keys())].rename(columns=available)

    df["survival_weeks"] = pd.to_numeric(df["survival_weeks"], errors="coerce")
    df["age"] = pd.to_numeric(df["age"], errors="coerce")
    df["who_status"] = pd.to_numeric(df["who_status"], errors="coerce")

    return df


def load_treatment(path: str = TREATMENT_TSV) -> pd.DataFrame:
    """
    Load and minimally clean the treatment data TSV.

    Parameters
    ----------
    path : str or Path
        Path to CFB-GBM_treatment_data TSV file.

    Returns
    -------
    pd.DataFrame
        Columns: patient_id, rt_dose_gy, n_fractions.

    Raises
    ------
    FileNotFoundError
        If the TSV file does not exist at the given path.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = df.columns.str.strip().str.lower().str.replace(" ", "_")

    rename = {
        "patient_id": "patient_id",
        "radiotherapy_dose": "rt_dose_gy",
        "fractions_number": "n_fractions",
    }
    available = {k: v for k, v in rename.items() if k in df.columns}
    df = df[list(available.keys())].rename(columns=available)

    df["rt_dose_gy"] = pd.to_numeric(df["rt_dose_gy"], errors="coerce")
    df["n_fractions"] = pd.to_numeric(df["n_fractions"], errors="coerce")

    return df


def load_imaging_availability(path: str = IMAGING_AVAILABILITY_TSV) -> pd.DataFrame:
    """
    Load the treatment imaging availability TSV and extract RTDOSE and GTV flags.

    Parameters
    ----------
    path : str or Path
        Path to CFB-GBM_treatment_imaging_availability TSV file.

    Returns
    -------
    pd.DataFrame
        Columns: patient_id, has_rtdose, has_gtv.

    Raises
    ------
    FileNotFoundError
        If the TSV file does not exist at the given path.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = df.columns.str.strip().str.lower().str.replace(" ", "_")

    rename = {
        "patient_id": "patient_id",
        "rtdose": "has_rtdose",
        "gtv": "has_gtv",
    }
    available = {k: v for k, v in rename.items() if k in df.columns}
    df = df[list(available.keys())].rename(columns=available)

    for col in ("has_rtdose", "has_gtv"):
        if col in df.columns:
            df[col] = df[col].str.strip().str.lower().map({"true": True, "1": True, "false": False, "0": False})

    return df


def build_cohort(
    clinical_path: str = CLINICAL_TSV,
    treatment_path: str = TREATMENT_TSV,
    imaging_path: str = IMAGING_AVAILABILITY_TSV,
    output_path: str = None,
) -> pd.DataFrame:
    """
    Merge clinical, treatment, and imaging availability data; apply inclusion
    criteria; add EQD2 column; export cohort table.

    Parameters
    ----------
    clinical_path : str or Path
        Path to clinical data TSV.
    treatment_path : str or Path
        Path to treatment data TSV.
    imaging_path : str or Path
        Path to imaging availability TSV.
    output_path : str or Path, optional
        If provided, export cohort table to this CSV path.
        Defaults to data/processed/cohort.csv.

    Returns
    -------
    pd.DataFrame
        Full cohort table including excluded patients.
        Column `included` (bool) flags eligible patients.
        Column `exclusion_reason` is empty string for included patients.
    """
    if output_path is None:
        output_path = DATA_PROCESSED / "cohort.csv"

    clinical = load_clinical(clinical_path)
    treatment = load_treatment(treatment_path)
    imaging = load_imaging_availability(imaging_path)

    # Merge on patient_id; outer join to preserve all patients for transparency
    cohort = clinical.merge(treatment, on="patient_id", how="outer")
    cohort = cohort.merge(imaging, on="patient_id", how="outer")

    # Fill missing availability flags as False
    for col in ("has_rtdose", "has_gtv"):
        if col not in cohort.columns:
            cohort[col] = False
        else:
            cohort[col] = cohort[col].fillna(False)

    # Determine inclusion/exclusion
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

    # Add EQD2 for included patients
    cohort["eqd2_gy"] = cohort.apply(
        lambda row: compute_eqd2(row["rt_dose_gy"], row["n_fractions"])
        if row["included"]
        else np.nan,
        axis=1,
    )

    # Reorder columns
    cols = [
        "patient_id", "rt_dose_gy", "n_fractions", "eqd2_gy",
        "survival_weeks", "age", "sex", "who_status",
        "has_rtdose", "has_gtv", "included", "exclusion_reason",
    ]
    cols = [c for c in cols if c in cohort.columns]
    cohort = cohort[cols].sort_values("patient_id").reset_index(drop=True)

    cohort.to_csv(output_path, index=False)
    print(f"Cohort saved to {output_path}")
    print(f"  Total patients : {len(cohort)}")
    print(f"  Included       : {cohort['included'].sum()}")
    print(f"  Excluded       : {(~cohort['included']).sum()}")

    excl = cohort[~cohort["included"]]["exclusion_reason"].value_counts()
    print("\nExclusion reasons:")
    for reason, count in excl.items():
        print(f"  {reason}: {count}")

    return cohort


if __name__ == "__main__":
    build_cohort()
