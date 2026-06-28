"""Expected on-disk paths for CFB-GBM RT NIfTI files (no nibabel dependency)."""

from pathlib import Path
from typing import Dict

from src.config import DATA_RAW


def _patient_dir(patient_id: str, data_dir: Path) -> Path:
    return data_dir / str(patient_id) / "t0"


def rtdose_path(patient_id: str, data_dir: Path = DATA_RAW) -> Path:
    return _patient_dir(patient_id, data_dir) / f"{patient_id}_t0_rtdose.nii.gz"


def gtv_path(patient_id: str, data_dir: Path = DATA_RAW) -> Path:
    return _patient_dir(patient_id, data_dir) / f"{patient_id}_t0_gtv.nii.gz"


def gtv_t1_path(patient_id: str, data_dir: Path = DATA_RAW) -> Path:
    return data_dir / str(patient_id) / "t1" / f"{patient_id}_t1_gtv.nii.gz"


def expected_nifti_paths(patient_id: str, data_dir: Path = DATA_RAW) -> Dict[str, Path]:
    return {"rtdose": rtdose_path(patient_id, data_dir), "gtv": gtv_path(patient_id, data_dir)}
