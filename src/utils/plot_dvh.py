"""
Publication-style DVH overlay plots for CFB-GBM cohort.

Loads cumulative DVH curves from per-patient NPZ files in
``data/processed/dvh_curves/`` or from a pre-built dictionary.
"""

from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from src.config import DVH_CURVES_DIR, FIGURES_DIR

DvhCurve = Tuple[np.ndarray, np.ndarray]

DEFAULT_RC_PARAMS = {
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
}

OVERLAY_COLORS = (
    "#2c7bb6",
    "#d7191c",
    "#fdae61",
    "#1a9641",
    "#762a83",
    "#9970ab",
    "#35978f",
    "#80cdc1",
)


def load_dvh_curve(patient_id: str, curves_dir: Path = DVH_CURVES_DIR) -> DvhCurve:
    """
    Load cumulative DVH for one patient from ``{patient_id}_dvh.npz``.

    Parameters
    ----------
    patient_id : str
        CFB-GBM patient identifier.
    curves_dir : pathlib.Path, optional
        Directory containing per-patient DVH NPZ files.

    Returns
    -------
    tuple of (np.ndarray, np.ndarray)
        ``(dose_gy, volume_pct)`` on a common 0–70 Gy grid when available;
        otherwise raw ``(dose_bins, volume_pct)`` from extraction.

    Raises
    ------
    FileNotFoundError
        If the NPZ file does not exist.
    KeyError
        If required arrays are missing from the archive.
    """
    path = curves_dir / f"{patient_id}_dvh.npz"
    if not path.exists():
        raise FileNotFoundError(f"DVH curve not found: {path}")

    data = np.load(path)
    if "dose_grid" in data.files and "volume_pct_grid" in data.files:
        return data["dose_grid"], data["volume_pct_grid"]
    if "dose_bins" in data.files and "volume_pct" in data.files:
        return data["dose_bins"], data["volume_pct"]
    raise KeyError(f"Unrecognised DVH NPZ layout in {path}: {data.files}")


def load_dvh_curves(
    patient_ids: Sequence[str],
    curves_dir: Path = DVH_CURVES_DIR,
) -> Dict[str, DvhCurve]:
    """
    Load cumulative DVH curves for multiple patients.

    Parameters
    ----------
    patient_ids : sequence of str
        Patient identifiers.
    curves_dir : pathlib.Path, optional
        Directory containing per-patient DVH NPZ files.

    Returns
    -------
    dict
        Maps each patient ID to ``(dose_gy, volume_pct)``.
    """
    return {pid: load_dvh_curve(pid, curves_dir=curves_dir) for pid in patient_ids}


def plot_dvh_overlay(
    patient_ids: Sequence[str],
    dvh_data: Optional[Mapping[str, DvhCurve]] = None,
    save_path: Optional[Union[str, Path]] = None,
    curves_dir: Path = DVH_CURVES_DIR,
    title: str = "Cumulative DVH overlay",
    xlabel: str = "Dose (Gy)",
    ylabel: str = "Volume receiving ≥ dose (%)",
    figsize: Tuple[float, float] = (7.0, 4.5),
    ax: Optional[plt.Axes] = None,
) -> Figure:
    """
    Plot cumulative DVH curves for selected patients.

    Parameters
    ----------
    patient_ids : sequence of str
        Patients to overlay (order preserved in legend).
    dvh_data : mapping, optional
        Pre-loaded ``{patient_id: (dose_gy, volume_pct)}``. When ``None``,
        curves are read from ``curves_dir``.
    save_path : str or pathlib.Path, optional
        If provided, figure is saved at 300 dpi (PNG).
    curves_dir : pathlib.Path, optional
        Source directory when ``dvh_data`` is ``None``.
    title : str, optional
        Axes title.
    xlabel, ylabel : str, optional
        Axis labels.
    figsize : tuple of float, optional
        Figure size in inches when creating a new axes.
    ax : matplotlib.axes.Axes, optional
        Existing axes; when ``None``, a new figure is created.

    Returns
    -------
    matplotlib.figure.Figure
        Parent figure (new or existing).

    Raises
    ------
    ValueError
        If ``patient_ids`` is empty.
    FileNotFoundError
        If a curve file is missing and ``dvh_data`` is not supplied.
    """
    if not patient_ids:
        raise ValueError("patient_ids must contain at least one ID")

    if dvh_data is None:
        dvh_data = load_dvh_curves(patient_ids, curves_dir=curves_dir)

    created_fig = ax is None
    if created_fig:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    with plt.rc_context(DEFAULT_RC_PARAMS):
        for index, patient_id in enumerate(patient_ids):
            if patient_id not in dvh_data:
                raise KeyError(f"patient_id {patient_id!r} missing from dvh_data")
            dose_gy, volume_pct = dvh_data[patient_id]
            color = OVERLAY_COLORS[index % len(OVERLAY_COLORS)]
            ax.plot(
                dose_gy,
                volume_pct,
                lw=1.8,
                color=color,
                label=f"Patient {patient_id}",
            )

        ax.set_xlim(left=0.0)
        ax.set_ylim(0.0, 100.0)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=9, loc="lower left", frameon=False)
        ax.grid(True, axis="both", alpha=0.25, lw=0.6)

        if save_path is not None:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path)

    return fig


def main() -> None:
    """Generate a sample DVH overlay for verification (standard vs hypofractionated)."""
    import pandas as pd

    from src.config import DATA_PROCESSED

    modeling = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    std_id = str(modeling.loc[modeling["rt_dose_gy"] == 60.0, "patient_id"].iloc[0])
    hypo_id = str(modeling.loc[modeling["rt_dose_gy"] == 40.05, "patient_id"].iloc[0])

    out = FIGURES_DIR / "02_dvh_overlay_sample.png"
    plot_dvh_overlay(
        [std_id, hypo_id],
        title=f"Sample DVH overlay (60 Gy patient {std_id} vs 40.05 Gy patient {hypo_id})",
        save_path=out,
    )
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
