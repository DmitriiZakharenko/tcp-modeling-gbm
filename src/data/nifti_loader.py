"""
NIfTI loader: load RTDOSE and GTV segmentation files for a single patient.

Expected file naming convention under data/raw/<patient_id>/:
    <patient_id>_t0_rtdose.nii.gz
    <patient_id>_t0_gtv.nii.gz
"""

from pathlib import Path
from typing import Tuple

import nibabel as nib
import numpy as np

from src.config import DATA_RAW


def load_rtdose(patient_id: str, data_dir: Path = DATA_RAW) -> Tuple[np.ndarray, np.ndarray, tuple]:
    """
    Load RTDOSE NIfTI file for a patient.

    Parameters
    ----------
    patient_id : str
        Patient identifier matching the directory and filename convention.
    data_dir : Path
        Root directory containing per-patient subdirectories.

    Returns
    -------
    dose_array : np.ndarray
        3D dose array in Gy.
    affine : np.ndarray
        4x4 affine matrix (voxel-to-world transformation).
    voxel_spacing_mm : tuple of float
        Voxel dimensions (dx, dy, dz) in mm.

    Raises
    ------
    FileNotFoundError
        If the expected NIfTI file does not exist.
    """
    raise NotImplementedError("Task EXTRACT-01: implement nifti_loader.load_rtdose")


def load_gtv_mask(patient_id: str, data_dir: Path = DATA_RAW) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load GTV segmentation mask NIfTI file for a patient.

    Parameters
    ----------
    patient_id : str
        Patient identifier.
    data_dir : Path
        Root directory containing per-patient subdirectories.

    Returns
    -------
    mask_array : np.ndarray
        3D binary mask (1 = GTV, 0 = background).
    affine : np.ndarray
        4x4 affine matrix.

    Raises
    ------
    FileNotFoundError
        If the expected NIfTI file does not exist.
    ValueError
        If the mask contains values other than 0 and 1.
    """
    raise NotImplementedError("Task EXTRACT-01: implement nifti_loader.load_gtv_mask")
