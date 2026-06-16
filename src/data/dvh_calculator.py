"""
DVH calculator: compute cumulative dose-volume histogram and extract
standard DVH metrics from a 3D dose array and a binary GTV mask.

Example
-------
    from src.data.nifti_loader import load_rtdose, load_gtv_mask
    from src.data.dvh_calculator import compute_dvh, extract_dvh_metrics

    dose, affine, spacing = load_rtdose("1")
    mask, _ = load_gtv_mask("1")
    metrics = extract_dvh_metrics(dose, mask, spacing)
    print(metrics)
"""

from typing import Dict, Tuple

import numpy as np


def voxel_volume_cc(spacing_mm: Tuple[float, float, float]) -> float:
    """
    Compute voxel volume in cubic centimetres from voxel spacing in mm.

    Parameters
    ----------
    spacing_mm : tuple of float
        Voxel dimensions (dx, dy, dz) in mm.

    Returns
    -------
    float
        Voxel volume in cm³.
    """
    return (spacing_mm[0] * spacing_mm[1] * spacing_mm[2]) / 1000.0


def compute_dvh(
    dose_array: np.ndarray,
    mask_array: np.ndarray,
    voxel_vol_cc: float,
    n_bins: int = 1000,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the cumulative DVH for a structure defined by a binary mask.

    Returns dose bins and cumulative volume fractions (% of structure volume
    receiving at least that dose).

    Parameters
    ----------
    dose_array : np.ndarray
        3D dose array in Gy. Must have the same shape as mask_array.
    mask_array : np.ndarray
        3D binary mask (True = inside structure).
    voxel_vol_cc : float
        Volume of a single voxel in cm³.
    n_bins : int
        Number of dose bins (default 1000).

    Returns
    -------
    dose_bins : np.ndarray, shape (n_bins,)
        Dose axis in Gy (left edge of each bin).
    volume_pct : np.ndarray, shape (n_bins,)
        Cumulative volume fraction in % receiving >= dose_bins[i].

    Raises
    ------
    ValueError
        If mask contains no voxels (empty structure).
    """
    dose_in_structure = dose_array[mask_array]

    if dose_in_structure.size == 0:
        raise ValueError("GTV mask is empty — no voxels inside the structure.")

    d_min = float(dose_in_structure.min())
    d_max = float(dose_in_structure.max())

    dose_bins = np.linspace(d_min, d_max, n_bins)
    volume_pct = np.array(
        [100.0 * np.mean(dose_in_structure >= d) for d in dose_bins],
        dtype=np.float32,
    )

    return dose_bins, volume_pct


def _dx(dose_in_structure: np.ndarray, x: float) -> float:
    """
    Compute Dx: minimum dose received by at least x% of the structure volume.

    Parameters
    ----------
    dose_in_structure : np.ndarray
        1D array of dose values within the structure in Gy.
    x : float
        Volume fraction in % (e.g. 95.0 for D95).

    Returns
    -------
    float
        Dx in Gy.
    """
    return float(np.percentile(dose_in_structure, 100.0 - x))


def extract_dvh_metrics(
    dose_array: np.ndarray,
    mask_array: np.ndarray,
    voxel_spacing_mm: Tuple[float, float, float],
) -> Dict[str, float]:
    """
    Extract standard DVH metrics for the GTV structure.

    Parameters
    ----------
    dose_array : np.ndarray
        3D dose array in Gy.
    mask_array : np.ndarray
        3D binary GTV mask.
    voxel_spacing_mm : tuple of float
        Voxel dimensions (dx, dy, dz) in mm.

    Returns
    -------
    dict with keys:
        D95_gy    : dose covering ≥95% of GTV volume (Gy)
        D98_gy    : dose covering ≥98% of GTV volume (Gy)
        D50_gy    : dose covering ≥50% of GTV volume — median dose (Gy)
        D2_gy     : dose covering ≥2% of GTV volume — near-max dose (Gy)
        Dmean_gy  : mean dose within GTV (Gy)
        Dmax_gy   : maximum dose within GTV (Gy)
        Dmin_gy   : minimum dose within GTV (Gy)
        volume_cc : GTV volume in cm³

    Raises
    ------
    ValueError
        If GTV mask is empty.
    """
    dose_in_structure = dose_array[mask_array]

    if dose_in_structure.size == 0:
        raise ValueError("GTV mask is empty — no voxels inside the structure.")

    vol_cc = voxel_volume_cc(voxel_spacing_mm) * float(mask_array.sum())

    return {
        "D95_gy":   _dx(dose_in_structure, 95.0),
        "D98_gy":   _dx(dose_in_structure, 98.0),
        "D50_gy":   _dx(dose_in_structure, 50.0),
        "D2_gy":    _dx(dose_in_structure, 2.0),
        "Dmean_gy": float(dose_in_structure.mean()),
        "Dmax_gy":  float(dose_in_structure.max()),
        "Dmin_gy":  float(dose_in_structure.min()),
        "volume_cc": vol_cc,
    }
