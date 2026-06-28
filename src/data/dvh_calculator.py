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

from typing import Dict, List, Tuple

import numpy as np

# Volume-percentile dose metrics (Dx = dose to x% of volume)
DX_VOLUME_PCT: Tuple[float, ...] = (2.0, 5.0, 10.0, 25.0, 50.0, 90.0, 95.0, 98.0)

# Volume receiving at least X Gy (% of structure)
VX_DOSE_GY: Tuple[float, ...] = (20.0, 40.0, 50.0, 60.0)

# gEUD exponent values (Niemierko: a=10 tumor, a=-10 serial, a=1 mean dose)
GEUD_A_VALUES: Tuple[float, ...] = (-10.0, 1.0, 10.0)

# Common dose grid (Gy) for exported cumulative DVH curves
DVH_DOSE_GRID_GY: np.ndarray = np.linspace(0.0, 70.0, 141, dtype=np.float32)

SCALAR_METRIC_KEYS: Tuple[str, ...] = (
    "D2_gy",
    "D5_gy",
    "D10_gy",
    "D25_gy",
    "D50_gy",
    "D90_gy",
    "D95_gy",
    "D98_gy",
    "V20_pct",
    "V40_pct",
    "V50_pct",
    "V60_pct",
    "Dmean_gy",
    "Dmax_gy",
    "Dmin_gy",
    "Dstd_gy",
    "HI_gy",
    "volume_cc",
    "gEUD_am10_gy",
    "gEUD_a1_gy",
    "gEUD_a10_gy",
)


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
    sorted_doses = np.sort(dose_in_structure)
    n_voxels = sorted_doses.size
    counts_ge = n_voxels - np.searchsorted(sorted_doses, dose_bins, side="left")
    volume_pct = (100.0 * counts_ge / n_voxels).astype(np.float32)

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


def compute_geud(dose_in_structure: np.ndarray, a: float) -> float:
    """
    Generalized equivalent uniform dose (Niemierko).

    EUD = (mean(D^a))^(1/a); for a→0 uses geometric mean.
    """
    doses = np.asarray(dose_in_structure, dtype=np.float64)
    doses = np.clip(doses, 0.0, None)
    if doses.size == 0:
        raise ValueError("Cannot compute gEUD on empty dose array.")
    if abs(a) < 1e-6:
        positive = doses[doses > 0]
        if positive.size == 0:
            return 0.0
        return float(np.exp(np.mean(np.log(positive))))
    if a < 0:
        positive = doses[doses > 0]
        if positive.size == 0:
            return 0.0
        doses = positive
    return float(np.power(np.mean(np.power(doses, a)), 1.0 / a))


def _geud_key(a: float) -> str:
    if a == -10.0:
        return "gEUD_am10_gy"
    if a == 1.0:
        return "gEUD_a1_gy"
    if a == 10.0:
        return "gEUD_a10_gy"
    label = str(a).replace("-", "m").replace(".", "p")
    return f"gEUD_a{label}_gy"


def interpolate_dvh_on_grid(
    dose_bins: np.ndarray,
    volume_pct: np.ndarray,
    dose_grid: np.ndarray = DVH_DOSE_GRID_GY,
) -> np.ndarray:
    """Interpolate cumulative DVH (% volume >= dose) onto a fixed dose grid."""
    return np.interp(
        dose_grid,
        dose_bins,
        volume_pct,
        left=100.0,
        right=0.0,
    ).astype(np.float32)


def extract_dvh_metrics(
    dose_array: np.ndarray,
    mask_array: np.ndarray,
    voxel_spacing_mm: Tuple[float, float, float],
) -> Dict[str, object]:
    """
    Extract DVH scalar metrics, cumulative DVH curve, and mid-axial slices.

    Returns
    -------
    dict with scalar metric keys (see SCALAR_METRIC_KEYS) plus:
        dose_bins, volume_pct — cumulative DVH arrays
        dose_slice, mask_slice — mid-axial 2D arrays for visualization
        slice_index — axial index used
    """
    dose_in_structure = dose_array[mask_array]

    if dose_in_structure.size == 0:
        raise ValueError("GTV mask is empty — no voxels inside the structure.")

    vol_cc = voxel_volume_cc(voxel_spacing_mm) * float(mask_array.sum())
    dose_bins, volume_pct = compute_dvh(dose_array, mask_array, voxel_volume_cc(voxel_spacing_mm))

    metrics: Dict[str, float] = {
        "D2_gy": _dx(dose_in_structure, 2.0),
        "D5_gy": _dx(dose_in_structure, 5.0),
        "D10_gy": _dx(dose_in_structure, 10.0),
        "D25_gy": _dx(dose_in_structure, 25.0),
        "D50_gy": _dx(dose_in_structure, 50.0),
        "D90_gy": _dx(dose_in_structure, 90.0),
        "D95_gy": _dx(dose_in_structure, 95.0),
        "D98_gy": _dx(dose_in_structure, 98.0),
        "Dmean_gy": float(dose_in_structure.mean()),
        "Dmax_gy": float(dose_in_structure.max()),
        "Dmin_gy": float(dose_in_structure.min()),
        "Dstd_gy": float(dose_in_structure.std()),
        "volume_cc": vol_cc,
    }

    for dose_gy in VX_DOSE_GY:
        metrics[f"V{int(dose_gy)}_pct"] = float(100.0 * np.mean(dose_in_structure >= dose_gy))

    d50 = metrics["D50_gy"]
    metrics["HI_gy"] = float((metrics["D2_gy"] - metrics["D98_gy"]) / d50) if d50 > 0 else 0.0

    for a in GEUD_A_VALUES:
        metrics[_geud_key(a)] = compute_geud(dose_in_structure, a)

    # Mid-axial slice through GTV centroid (for notebooks without full 3D NIfTI)
    z_indices = np.where(mask_array.any(axis=(0, 1)))[0]
    if z_indices.size == 0:
        slice_index = int(mask_array.shape[2] // 2)
    else:
        slice_index = int(z_indices[len(z_indices) // 2])

    metrics["dose_bins"] = dose_bins
    metrics["volume_pct"] = volume_pct
    metrics["dose_slice"] = np.asarray(dose_array[:, :, slice_index], dtype=np.float32)
    metrics["mask_slice"] = np.asarray(mask_array[:, :, slice_index], dtype=np.uint8)
    metrics["slice_index"] = float(slice_index)

    return metrics
