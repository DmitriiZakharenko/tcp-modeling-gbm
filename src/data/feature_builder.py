"""
Feature builder: extract DVH metrics for all included patients.

Computes extended dosimetric scalars (Dx, Vx, gEUD, HI), saves per-patient
cumulative DVH curves and mid-axial dose/GTV slices for notebooks without raw NIfTI.

Usage
-----
    python -m src.data.feature_builder
    python -m src.data.feature_builder --workers 4 --force
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

from src.config import (
    COHORT_CSV,
    DATA_RAW,
    DOSE_SLICES_DIR,
    DVH_CACHE_DIR,
    DVH_CURVES_DIR,
    DVH_CURVES_NPZ,
    FEATURES_CSV,
)
from src.data.dvh_calculator import (
    DVH_DOSE_GRID_GY,
    SCALAR_METRIC_KEYS,
    extract_dvh_metrics,
    interpolate_dvh_on_grid,
)
from src.data.nifti_loader import check_shape_match, load_gtv_mask, load_rtdose

CACHE_VERSION = "v2"


def _cache_path(patient_id: str, cache_dir: Path) -> Path:
    return cache_dir / f"{patient_id}_dvh_{CACHE_VERSION}.npz"


def _split_bundle(bundle: Dict[str, Any]) -> tuple[Dict[str, float], Dict[str, Any]]:
    """Separate scalar metrics from arrays for CSV vs NPZ storage."""
    scalars = {key: float(bundle[key]) for key in SCALAR_METRIC_KEYS if key in bundle}
    extras = {
        key: bundle[key]
        for key in ("dose_bins", "volume_pct", "dose_slice", "mask_slice", "slice_index")
        if key in bundle
    }
    return scalars, extras


def _save_patient_artifacts(
    patient_id: str,
    scalars: Dict[str, float],
    extras: Dict[str, Any],
    cache_dir: Path,
    curves_dir: Path,
    slices_dir: Path,
) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    curves_dir.mkdir(parents=True, exist_ok=True)
    slices_dir.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        _cache_path(patient_id, cache_dir),
        **scalars,
        dose_bins=extras["dose_bins"],
        volume_pct=extras["volume_pct"],
        dose_grid=DVH_DOSE_GRID_GY,
        volume_pct_grid=interpolate_dvh_on_grid(extras["dose_bins"], extras["volume_pct"]),
    )

    np.savez_compressed(
        curves_dir / f"{patient_id}_dvh.npz",
        dose_bins=extras["dose_bins"],
        volume_pct=extras["volume_pct"],
        dose_grid=DVH_DOSE_GRID_GY,
        volume_pct_grid=interpolate_dvh_on_grid(extras["dose_bins"], extras["volume_pct"]),
    )

    np.savez_compressed(
        slices_dir / f"{patient_id}_axial_mid.npz",
        dose_slice=extras["dose_slice"],
        mask_slice=extras["mask_slice"],
        slice_index=int(extras["slice_index"]),
    )


def _load_cached_scalars(patient_id: str, cache_dir: Path) -> Optional[Dict[str, float]]:
    path = _cache_path(patient_id, cache_dir)
    if not path.exists():
        return None
    data = np.load(path)
    if not all(key in data.files for key in SCALAR_METRIC_KEYS):
        return None
    return {key: float(data[key]) for key in SCALAR_METRIC_KEYS}


def extract_patient_features(
    patient_id: str,
    data_dir: Path,
    cache_dir: Path,
    curves_dir: Path,
    slices_dir: Path,
    use_cache: bool = True,
) -> Dict[str, float]:
    """Load NIfTI files and extract full DVH feature set for one patient."""
    if use_cache:
        cached = _load_cached_scalars(patient_id, cache_dir)
        if cached is not None:
            return {"patient_id": patient_id, **cached}

    dose, _, spacing = load_rtdose(patient_id, data_dir=data_dir, mmap=True)
    mask, _ = load_gtv_mask(patient_id, data_dir=data_dir, mmap=True)
    check_shape_match(dose, mask, patient_id)

    bundle = extract_dvh_metrics(dose, mask, spacing)
    scalars, extras = _split_bundle(bundle)
    _save_patient_artifacts(patient_id, scalars, extras, cache_dir, curves_dir, slices_dir)

    return {"patient_id": patient_id, **scalars}


def _process_patient_task(args: tuple) -> Dict[str, float]:
    """Wrapper for ProcessPoolExecutor (must be top-level for pickling)."""
    return extract_patient_features(*args)


def _load_existing_features(features_path: Path) -> pd.DataFrame:
    if not features_path.exists():
        return pd.DataFrame(columns=["patient_id", *SCALAR_METRIC_KEYS])
    frame = pd.read_csv(features_path)
    if all(key in frame.columns for key in SCALAR_METRIC_KEYS):
        return frame
    return pd.DataFrame(columns=["patient_id", *SCALAR_METRIC_KEYS])


def _bundle_dvh_curves(patient_ids: List[str], curves_dir: Path, output_path: Path) -> None:
    """Pack all patient DVH curves on a common dose grid into one NPZ file."""
    dose_grid = DVH_DOSE_GRID_GY
    volume_matrix = np.zeros((len(patient_ids), dose_grid.size), dtype=np.float32)

    for index, patient_id in enumerate(patient_ids):
        curve_file = curves_dir / f"{patient_id}_dvh.npz"
        if not curve_file.exists():
            cache_file = _cache_path(patient_id, DVH_CACHE_DIR)
            if not cache_file.exists():
                raise FileNotFoundError(f"Missing DVH curve for patient {patient_id}")
            data = np.load(cache_file)
        else:
            data = np.load(curve_file)
        volume_matrix[index] = data["volume_pct_grid"]

    np.savez_compressed(
        output_path,
        patient_id=np.array(patient_ids, dtype=str),
        dose_grid_gy=dose_grid,
        volume_pct_ge=volume_matrix,
    )
    print(f"Combined DVH curves: {output_path} ({len(patient_ids)} patients)")


def build_features(
    cohort_path: Path = COHORT_CSV,
    data_dir: Path = DATA_RAW,
    output_path: Path = FEATURES_CSV,
    cache_dir: Path = DVH_CACHE_DIR,
    curves_dir: Path = DVH_CURVES_DIR,
    slices_dir: Path = DOSE_SLICES_DIR,
    curves_bundle_path: Path = DVH_CURVES_NPZ,
    workers: int = 1,
    use_cache: bool = True,
    force: bool = False,
) -> pd.DataFrame:
    """Extract full DVH feature set for all included patients."""
    if not cohort_path.exists():
        raise FileNotFoundError(
            f"Cohort file not found: {cohort_path}\n"
            "Run: python -m src.data.cohort_builder"
        )

    cohort = pd.read_csv(cohort_path)
    patient_ids = cohort.loc[cohort["included"] == True, "patient_id"].astype(str).tolist()  # noqa: E712

    existing = _load_existing_features(output_path)
    if force:
        pending_ids = patient_ids
    else:
        done_ids = set(existing["patient_id"].astype(str))
        pending_ids = [pid for pid in patient_ids if pid not in done_ids]

    print(f"Feature extraction: {len(patient_ids)} included patients")
    print(f"  Scalar metrics: {len(SCALAR_METRIC_KEYS)} per patient")
    print(f"  Already in {output_path.name}: {len(patient_ids) - len(pending_ids)}")
    print(f"  To process: {len(pending_ids)}")
    print(f"  Workers: {workers}")

    if pending_ids:
        new_rows: List[Dict[str, float]] = []
        errors: List[str] = []

        task_args = [
            (pid, data_dir, cache_dir, curves_dir, slices_dir, use_cache and not force)
            for pid in pending_ids
        ]

        if workers <= 1:
            for args in tqdm(task_args, total=len(pending_ids), desc="DVH extraction"):
                patient_id = args[0]
                try:
                    new_rows.append(_process_patient_task(args))
                except Exception as exc:
                    errors.append(f"patient {patient_id}: {exc}")
        else:
            with ProcessPoolExecutor(max_workers=workers) as executor:
                futures = {
                    executor.submit(_process_patient_task, args): args[0]
                    for args in task_args
                }
                for future in tqdm(as_completed(futures), total=len(futures), desc="DVH extraction"):
                    patient_id = futures[future]
                    try:
                        new_rows.append(future.result())
                    except Exception as exc:
                        errors.append(f"patient {patient_id}: {exc}")

        if new_rows:
            new_df = pd.DataFrame(new_rows)
            new_df["patient_id"] = new_df["patient_id"].astype(str)
            if not existing.empty and not force:
                existing = existing.copy()
                existing["patient_id"] = existing["patient_id"].astype(str)
                features = pd.concat([existing, new_df], ignore_index=True)
            else:
                features = new_df
            features["patient_id"] = features["patient_id"].astype(str)
            features = features.drop_duplicates(subset=["patient_id"], keep="last")
            features = features.sort_values("patient_id", key=lambda s: s.astype(int)).reset_index(drop=True)
            features.to_csv(output_path, index=False)
            print(f"\nFeatures saved to: {output_path} ({len(features)} patients)")
        else:
            features = existing
            print("\nNo new features extracted.")

        if errors:
            print(f"\n{len(errors)} error(s):")
            for msg in errors[:20]:
                print(f"  {msg}")
            if len(errors) > 20:
                print(f"  ... and {len(errors) - 20} more")
            raise RuntimeError(f"{len(errors)} patient(s) failed during feature extraction")
    else:
        features = existing
        print("All patients already in features.csv — rebuilding DVH bundle only.")

    processed_ids = features["patient_id"].astype(str).tolist()
    _bundle_dvh_curves(processed_ids, curves_dir, curves_bundle_path)
    return features


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract DVH features for included patients")
    parser.add_argument("--cohort", type=Path, default=COHORT_CSV)
    parser.add_argument("--data-dir", type=Path, default=DATA_RAW)
    parser.add_argument("--output", type=Path, default=FEATURES_CSV)
    parser.add_argument("--cache-dir", type=Path, default=DVH_CACHE_DIR)
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel worker processes (default: 1)",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Disable per-patient DVH metric cache",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Reprocess all patients",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    try:
        build_features(
            cohort_path=args.cohort,
            data_dir=args.data_dir,
            output_path=args.output,
            cache_dir=args.cache_dir,
            workers=args.workers,
            use_cache=not args.no_cache,
            force=args.force,
        )
    except (FileNotFoundError, RuntimeError) as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)
