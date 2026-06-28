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

# Clinical TSV filenames (stored in data/processed/ after download)
CLINICAL_TSV: Path = DATA_PROCESSED / "CFB-GBM_clinical_data_v02_20260129.tsv"
TREATMENT_TSV: Path = DATA_PROCESSED / "CFB-GBM_treatment_data_v02_20260129.tsv"
IMAGING_AVAILABILITY_TSV: Path = DATA_PROCESSED / "CFB-GBM_treatment_imaging_availability_v02_20260129.tsv"
COLUMNS_DESCRIPTION_TSV: Path = DATA_PROCESSED / "CFB-GBM_columns_description_new_v02_20260129.tsv"
COHORT_CSV: Path = DATA_PROCESSED / "cohort.csv"
FEATURES_CSV: Path = DATA_PROCESSED / "features.csv"
RAW_MANIFEST_CSV: Path = DATA_PROCESSED / "raw_data_manifest.csv"

TCIA_BASE_URL = "https://www.cancerimagingarchive.net/wp-content/uploads"


class ClinicalFile(NamedTuple):
    """Local path and remote URL for a CFB-GBM clinical/metadata TSV."""

    path: Path
    url: str


CLINICAL_FILES: tuple[ClinicalFile, ...] = (
    ClinicalFile(CLINICAL_TSV, f"{TCIA_BASE_URL}/CFB-GBM_clinical_data_v02_20260129.tsv"),
    ClinicalFile(TREATMENT_TSV, f"{TCIA_BASE_URL}/CFB-GBM_treatment_data_v02_20260129.tsv"),
    ClinicalFile(
        IMAGING_AVAILABILITY_TSV,
        f"{TCIA_BASE_URL}/CFB-GBM_treatment_imaging_availability_v02_20260129.tsv",
    ),
    ClinicalFile(
        COLUMNS_DESCRIPTION_TSV,
        f"{TCIA_BASE_URL}/CFB-GBM_columns_description_new_v02_20260129.tsv",
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

# Radiobiological constants
ALPHA_BETA_GBM: float = 10.0  # Gy; standard value for GBM (high alpha/beta tumor)

# Random seed for reproducibility
RANDOM_SEED: int = 42
