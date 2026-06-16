"""
Central configuration: project-wide paths and constants.

All paths are resolved relative to the repository root using pathlib.Path.
Override DATA_RAW by setting the environment variable TCP_DATA_RAW.
"""

import os
from pathlib import Path

# Repository root (two levels up from this file: src/config.py -> src/ -> root)
ROOT_DIR: Path = Path(__file__).resolve().parent.parent

# Data directories
DATA_RAW: Path = Path(os.environ.get("TCP_DATA_RAW", ROOT_DIR / "data" / "raw"))
DATA_PROCESSED: Path = ROOT_DIR / "data" / "processed"

# Output directories
FIGURES_DIR: Path = ROOT_DIR / "figures"
REPORTS_DIR: Path = ROOT_DIR / "reports"
NOTEBOOKS_DIR: Path = ROOT_DIR / "notebooks"

# Ensure output directories exist at import time
for _dir in (DATA_PROCESSED, FIGURES_DIR, REPORTS_DIR):
    _dir.mkdir(parents=True, exist_ok=True)

# Clinical TSV filenames (stored in data/processed/ after download)
CLINICAL_TSV: Path = DATA_PROCESSED / "CFB-GBM_clinical_data_v02_20260129.tsv"
TREATMENT_TSV: Path = DATA_PROCESSED / "CFB-GBM_treatment_data_v02_20260129.tsv"
IMAGING_AVAILABILITY_TSV: Path = DATA_PROCESSED / "CFB-GBM_treatment_imaging_availability_v02_20260129.tsv"

# Radiobiological constants
ALPHA_BETA_GBM: float = 10.0  # Gy; standard value for GBM (high alpha/beta tumor)

# Random seed for reproducibility
RANDOM_SEED: int = 42
