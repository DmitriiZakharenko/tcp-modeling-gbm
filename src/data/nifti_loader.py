"""
NIfTI loader: load RTDOSE and GTV segmentation files for a single patient.

File naming convention (CFB-GBM dataset):
    data/raw/<patient_id>/t0/<patient_id>_t0_rtdose.nii.gz
    data/raw/<patient_id>/t0/<patient_id>_t0_gtv.nii.gz

Example
-------
    from src.data.nifti_loader import load_rtdose, load_gtv_mask

    dose, affine, spacing = load_rtdose("1")
    mask, _ = load_gtv_mask("1")
"""

from pathlib import Path
from typing import Tuple

import nibabel as nib
import numpy as np

from src.config import DATA_RAW


def _patient_dir(patient_id: str, data_dir: Path) -> Path:
    """Return the t0 subdirectory for a given patient."""
    return data_dir / str(patient_id) / "t0"


def _voxel_spacing(nii_img: nib.Nifti1Image) -> Tuple[float, float, float]:
    """
    Extract voxel spacing in mm from a NIfTI image header.

    Parameters
    ----------
    nii_img : nib.Nifti1Image
        Loaded NIfTI image.

    Returns
    -------
    tuple of float
        (dx, dy, dz) voxel dimensions in mm.
    """
    zooms = nii_img.header.get_zooms()
    return float(zooms[0]), float(zooms[1]), float(zooms[2])


def load_rtdose(
    patient_id: str,
    data_dir: Path = DATA_RAW,
) -> Tuple[np.ndarray, np.ndarray, Tuple[float, float, float]]:
    """
    Load RTDOSE NIfTI file for a patient.

    Parameters
    ----------
    patient_id : str
        Patient identifier (e.g. "1", "42").
    data_dir : Path
        Root raw data directory containing per-patient subdirectories.

    Returns
    -------
    dose_array : np.ndarray, shape (X, Y, Z)
        3D dose array in Gy (float32).
    affine : np.ndarray, shape (4, 4)
        Voxel-to-world affine transformation matrix.
    voxel_spacing_mm : tuple of float
        Voxel dimensions (dx, dy, dz) in mm.

    Raises
    ------
    FileNotFoundError
        If the expected NIfTI file does not exist.
    """
    path = _patient_dir(patient_id, data_dir) / f"{patient_id}_t0_rtdose.nii.gz"
    if not path.exists():
        raise FileNotFoundError(f"RTDOSE not found for patient {patient_id}: {path}")

    nii = nib.load(str(path))
    dose = np.asarray(nii.dataobj, dtype=np.float32)
    return dose, nii.affine, _voxel_spacing(nii)


def load_gtv_mask(
    patient_id: str,
    data_dir: Path = DATA_RAW,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load GTV segmentation mask NIfTI file for a patient.

    Parameters
    ----------
    patient_id : str
        Patient identifier.
    data_dir : Path
        Root raw data directory.

    Returns
    -------
    mask_array : np.ndarray, shape (X, Y, Z)
        3D binary mask (True = GTV voxel, False = background).
    affine : np.ndarray, shape (4, 4)
        Voxel-to-world affine transformation matrix.

    Raises
    ------
    FileNotFoundError
        If the expected NIfTI file does not exist.
    ValueError
        If the mask contains values other than 0 and 1.
    """
    path = _patient_dir(patient_id, data_dir) / f"{patient_id}_t0_gtv.nii.gz"
    if not path.exists():
        raise FileNotFoundError(f"GTV mask not found for patient {patient_id}: {path}")

    nii = nib.load(str(path))
    raw = np.asarray(nii.dataobj, dtype=np.float32)

    unique = np.unique(raw)
    non_binary = unique[~np.isin(unique, [0.0, 1.0])]
    if len(non_binary) > 0:
        raise ValueError(
            f"GTV mask for patient {patient_id} contains unexpected values: {non_binary}. "
            "Expected binary (0/1) mask."
        )

    return raw.astype(bool), nii.affine


def check_shape_match(dose: np.ndarray, mask: np.ndarray, patient_id: str) -> None:
    """
    Verify that dose array and GTV mask have the same shape.

    Parameters
    ----------
    dose : np.ndarray
        Dose array loaded via load_rtdose.
    mask : np.ndarray
        GTV mask loaded via load_gtv_mask.
    patient_id : str
        Used in the error message.

    Raises
    ------
    ValueError
        If shapes do not match.
    """
    if dose.shape != mask.shape:
        raise ValueError(
            f"Shape mismatch for patient {patient_id}: "
            f"dose {dose.shape} vs mask {mask.shape}. "
            "Dose and GTV must be co-registered and in the same voxel space."
        )
