"""
DVH calculator: compute cumulative dose-volume histogram and extract DVH metrics
from a 3D dose array and a binary structure mask.
"""

from typing import Dict

import numpy as np


def compute_dvh(
    dose_array: np.ndarray,
    mask_array: np.ndarray,
    voxel_volume_cc: float,
    n_bins: int = 1000,
) -> Dict[str, np.ndarray]:
    """
    Compute the cumulative DVH for a structure defined by a binary mask.

    Parameters
    ----------
    dose_array : np.ndarray
        3D dose array in Gy. Must have the same shape as mask_array.
    mask_array : np.ndarray
        3D binary mask (1 = inside structure, 0 = outside).
    voxel_volume_cc : float
        Volume of a single voxel in cubic centimetres (cm³).
    n_bins : int
        Number of dose bins for the histogram (default 1000).

    Returns
    -------
    dict with keys:
        dose_bins  : np.ndarray, shape (n_bins,) — dose axis in Gy
        volume_pct : np.ndarray, shape (n_bins,) — cumulative volume in %
    """
    raise NotImplementedError("Task EXTRACT-02: implement dvh_calculator.compute_dvh")


def extract_dvh_metrics(
    dose_array: np.ndarray,
    mask_array: np.ndarray,
    voxel_spacing_mm: tuple,
) -> Dict[str, float]:
    """
    Extract standard DVH metrics for a structure.

    Parameters
    ----------
    dose_array : np.ndarray
        3D dose array in Gy.
    mask_array : np.ndarray
        3D binary mask.
    voxel_spacing_mm : tuple of float
        Voxel dimensions (dx, dy, dz) in mm.

    Returns
    -------
    dict with keys:
        D95_gy     : dose covering 95% of the structure volume (Gy)
        D98_gy     : dose covering 98% of the structure volume (Gy)
        D50_gy     : dose covering 50% of the structure volume (Gy)
        D2_gy      : dose covering 2% of the structure volume (Gy)
        Dmean_gy   : mean dose within the structure (Gy)
        Dmax_gy    : maximum dose within the structure (Gy)
        volume_cc  : structure volume in cm³
    """
    raise NotImplementedError("Task EXTRACT-02: implement dvh_calculator.extract_dvh_metrics")
