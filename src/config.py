"""
Central configuration: project-wide paths and constants.

All paths are resolved relative to the repository root using pathlib.Path.
Override DATA_RAW by setting the environment variable TCP_DATA_RAW.
"""

import os
from pathlib import Path
from typing import NamedTuple

# Repository root (two levels up from this file: src/config.py -> src/ -> root)
ROOT_DIR: Path = Path(__file__).resolve().parent.parent

# Data directories
DATA_RAW: Path = Path(os.environ.get("TCP_DATA_RAW", ROOT_DIR / "data" / "raw"))
DATA_PROCESSED: Path = ROOT_DIR / "data" / "processed"
DVH_CACHE_DIR: Path = DATA_PROCESSED / "dvh_cache"
DVH_CURVES_DIR: Path = DATA_PROCESSED / "dvh_curves"
DOSE_SLICES_DIR: Path = DATA_PROCESSED / "dose_slices"
DVH_CURVES_NPZ: Path = DATA_PROCESSED / "dvh_curves_all.npz"

# Output directories
FIGURES_DIR: Path = ROOT_DIR / "figures"
REPORTS_DIR: Path = ROOT_DIR / "reports"
NOTEBOOKS_DIR: Path = ROOT_DIR / "notebooks"

# Ensure output directories exist at import time
for _dir in (DATA_PROCESSED, FIGURES_DIR, REPORTS_DIR, DVH_CACHE_DIR, DVH_CURVES_DIR, DOSE_SLICES_DIR):
    _dir.mkdir(parents=True, exist_ok=True)

TCIA_BASE_URL = "https://www.cancerimagingarchive.net/wp-content/uploads"

# Clinical / metadata TSV filenames (stored in data/processed/ after download)
CLINICAL_TSV: Path = DATA_PROCESSED / "CFB-GBM_clinical_data_v03_20260619.tsv"
TREATMENT_TSV: Path = DATA_PROCESSED / "CFB-GBM_treatment_data_v02_20260129.tsv"
IMAGING_AVAILABILITY_TSV: Path = DATA_PROCESSED / "CFB-GBM_treatment_imaging_availability_v03_20260619.tsv"
MRI_AVAILABILITY_TSV: Path = DATA_PROCESSED / "CFB-GBM_mri_availability_v02_20260129.tsv"
CT_AVAILABILITY_TSV: Path = DATA_PROCESSED / "CFB-GBM_ct_availability_v02_20260129.tsv"
RANO_TSV: Path = DATA_PROCESSED / "CFB-GBM_rano_criteria_v03_20260619.tsv"
COLUMNS_DESCRIPTION_TSV: Path = DATA_PROCESSED / "CFB-GBM_column_descriptions_v03_20260619.tsv"
PYRADIOMICS_TSV: Path = DATA_PROCESSED / "CFB_GBM_features_extraction_pyradiomics_v03_20260619.tsv"
COHORT_CSV: Path = DATA_PROCESSED / "cohort.csv"
FEATURES_CSV: Path = DATA_PROCESSED / "features.csv"
RAW_MANIFEST_CSV: Path = DATA_PROCESSED / "raw_data_manifest.csv"


class ClinicalFile(NamedTuple):
    """Local path and remote URL for a CFB-GBM clinical/metadata TSV."""

    path: Path
    url: str
    required: bool = True


CLINICAL_FILES: tuple[ClinicalFile, ...] = (
    ClinicalFile(CLINICAL_TSV, f"{TCIA_BASE_URL}/CFB-GBM_clinical_data_v03_20260619.tsv"),
    ClinicalFile(TREATMENT_TSV, f"{TCIA_BASE_URL}/CFB-GBM_treatment_data_v02_20260129.tsv"),
    ClinicalFile(
        IMAGING_AVAILABILITY_TSV,
        f"{TCIA_BASE_URL}/CFB-GBM_treatment_imaging_availability_v03_20260619.tsv",
    ),
    ClinicalFile(MRI_AVAILABILITY_TSV, f"{TCIA_BASE_URL}/CFB-GBM_mri_availability_v02_20260129.tsv"),
    ClinicalFile(CT_AVAILABILITY_TSV, f"{TCIA_BASE_URL}/CFB-GBM_ct_availability_v02_20260129.tsv"),
    ClinicalFile(RANO_TSV, f"{TCIA_BASE_URL}/CFB-GBM_rano_criteria_v03_20260619.tsv"),
    ClinicalFile(
        COLUMNS_DESCRIPTION_TSV,
        f"{TCIA_BASE_URL}/CFB-GBM_column_descriptions_v03_20260619.tsv",
    ),
    ClinicalFile(
        PYRADIOMICS_TSV,
        f"{TCIA_BASE_URL}/CFB_GBM_features_extraction_pyradiomics_v03_20260619.tsv",
        required=False,
    ),
)

# Aspera / TCIA Faspex package (public credentials embedded in TCIA download URL)
FASPEX_HOST = "faspex.cancerimagingarchive.net"
FASPEX_API_BASE = f"https://{FASPEX_HOST}/aspera/faspex"
FASPEX_PACKAGE_ID = "1196"
FASPEX_PACKAGE_SLUG = "CFB-GBM"
FASPEX_PACKAGE_PATH = f"/packages/{FASPEX_PACKAGE_ID}/{FASPEX_PACKAGE_SLUG}"
FASPEX_PASSCODE = "77a6f7a7258b1a4ef9bf2d13fcd3f1ec87b9c313"
FASPEX_USER = "anon"
FASPEX_CLIENT_ID = "ff9aa63a-72e1-436f-82ef-5677eb1f7aee"
FASPEX_REDIRECT_URI = f"{FASPEX_API_BASE}/token"
FASPEX_PORT = 33001          # TCIA Faspex SSH port (direct ascp only)
FASPEX_UDP_PORT = 33001      # FASP data transfer port (direct ascp only)
FASPEX_FALLBACK_PORTS = (33001, 22)
ASPERA_MAX_RATE_MBPS = 200

RTDOSE_GLOB = "**/*_t0_rtdose.nii.gz"
GTV_GLOB = "**/*_t0_gtv.nii.gz"
GTV_T1_GLOB = "**/*_t1_gtv.nii.gz"

# TCIA Faspex v3 package (June 2026) — follow-up GTV / full imaging
FASPEX_PACKAGE_ID_V3 = "1303"
FASPEX_PACKAGE_SLUG_V3 = "CFB-GBM"
FASPEX_PASSCODE_V3 = "757109fecdc9c7f35da8e99badb3b595bd120522"
FASPEX_PACKAGE_PATH_V3 = f"/packages/{FASPEX_PACKAGE_ID_V3}/{FASPEX_PACKAGE_SLUG_V3}"

# Radiobiological constants
ALPHA_BETA_GBM: float = 10.0  # Gy; standard value for GBM (high alpha/beta tumor)

# Random seed for reproducibility
RANDOM_SEED: int = 42
