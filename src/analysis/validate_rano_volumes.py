"""
Validate RANO TSV volumes against NIfTI GTV segmentations (t0 and t1).
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
from scipy import stats

from src.config import DATA_PROCESSED, DATA_RAW, REPORTS_DIR
from src.data.dvh_calculator import voxel_volume_cc
from src.data.nifti_paths import gtv_path, gtv_t1_path

import nibabel as nib


def mask_volume_cc(mask_path: Path) -> float:
    """Compute structure volume in cm³ from a binary GTV NIfTI."""
    nii = nib.load(str(mask_path))
    data = np.asarray(nii.dataobj)
    spacing = nii.header.get_zooms()[:3]
    vox_cc = voxel_volume_cc((float(spacing[0]), float(spacing[1]), float(spacing[2])))
    return float((data > 0).sum() * vox_cc)


def validate_nifti_vs_rano(
    frame: pd.DataFrame,
    data_dir: Path = DATA_RAW,
    temporality: str = "t0",
) -> pd.DataFrame:
    """Per-patient NIfTI GTV volume vs RANO size column."""
    size_col = f"size_{temporality}_cm3"
    path_fn = gtv_path if temporality == "t0" else gtv_t1_path
    rows: List[dict] = []
    for _, row in frame.iterrows():
        pid = str(int(row["patient_id"]))
        rano_vol = row.get(size_col)
        if pd.isna(rano_vol):
            continue
        nii_path = path_fn(pid, data_dir)
        if not nii_path.exists():
            continue
        nifti_vol = mask_volume_cc(nii_path)
        rows.append(
            {
                "patient_id": pid,
                "temporality": temporality,
                "nifti_volume_cc": nifti_vol,
                "rano_size_cm3": float(rano_vol),
                "abs_diff_cc": abs(nifti_vol - float(rano_vol)),
                "pct_diff": abs(nifti_vol - float(rano_vol)) / max(float(rano_vol), 0.1) * 100,
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    rho, p = stats.spearmanr(df["nifti_volume_cc"], df["rano_size_cm3"])
    df.attrs["spearman_rho"] = float(rho)
    df.attrs["spearman_p"] = float(p)
    df.attrs["n"] = len(df)
    return df


def run_volume_validation(output_dir: Optional[Path] = None) -> pd.DataFrame:
    """Validate t0/t1 NIfTI GTV volumes vs RANO TSV where files exist."""
    frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    out_dir = output_dir or REPORTS_DIR / "metrics"
    out_dir.mkdir(parents=True, exist_ok=True)

    summaries = []
    for tp in ("t0", "t1"):
        detail = validate_nifti_vs_rano(frame, temporality=tp)
        if detail.empty:
            summaries.append({"temporality": tp, "n": 0, "spearman_rho": np.nan, "spearman_p": np.nan, "median_pct_diff": np.nan})
            continue
        detail.to_csv(out_dir / f"rano_nifti_volume_detail_{tp}.csv", index=False)
        summaries.append(
            {
                "temporality": tp,
                "n": detail.attrs.get("n", len(detail)),
                "spearman_rho": detail.attrs.get("spearman_rho"),
                "spearman_p": detail.attrs.get("spearman_p"),
                "median_pct_diff": float(detail["pct_diff"].median()),
            }
        )
    summary = pd.DataFrame(summaries)
    summary.to_csv(out_dir / "rano_nifti_volume_summary.csv", index=False)
    return summary


def main() -> None:
    summary = run_volume_validation()
    print(summary.to_string(index=False))
    t1_n = int(summary.loc[summary["temporality"] == "t1", "n"].iloc[0]) if len(summary) > 1 else 0
    if t1_n == 0:
        print("\nNo t1 GTV NIfTI found. Download with: python -m src.data.download_rt_files --include-t1-gtv")


if __name__ == "__main__":
    main()
