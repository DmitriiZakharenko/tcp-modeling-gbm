"""
Verify completeness of downloaded RT NIfTI files against the cohort table.

Cross-references cohort.csv with files on disk and optionally writes a manifest.

Usage
-----
    python -m src.data.verify_raw_data
    python -m src.data.verify_raw_data --write-manifest
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

from src.config import COHORT_CSV, DATA_RAW, RAW_MANIFEST_CSV
from src.data.nifti_paths import expected_nifti_paths
from src.data.organize_raw_data import find_unorganized_nifti


def verify_raw_data(
    cohort_path: Path = COHORT_CSV,
    data_dir: Path = DATA_RAW,
    included_only: bool = True,
    write_manifest: bool = False,
    manifest_path: Path = RAW_MANIFEST_CSV,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Check which patients have RTDOSE and GTV NIfTI files on disk.

    Parameters
    ----------
    cohort_path : Path
        Path to cohort.csv.
    data_dir : Path
        Root directory containing per-patient NIfTI subdirectories.
    included_only : bool
        If True, verify only patients with included == True.
    write_manifest : bool
        If True, write a manifest CSV to manifest_path.
    manifest_path : Path
        Output path for the manifest CSV.

    Returns
    -------
    pd.DataFrame
        One row per patient with columns:
        patient_id, has_rtdose_file, has_gtv_file, complete, included.

    Raises
    ------
    FileNotFoundError
        If cohort_path does not exist.
    """
    if not cohort_path.exists():
        raise FileNotFoundError(
            f"Cohort file not found: {cohort_path}\n"
            "Run: python -m src.data.cohort_builder"
        )

    cohort = pd.read_csv(cohort_path)
    if included_only:
        cohort = cohort[cohort["included"] == True]  # noqa: E712

    rows = []
    for _, patient_row in cohort.iterrows():
        patient_id = str(patient_row["patient_id"])
        paths = expected_nifti_paths(patient_id, data_dir)
        has_rtdose = paths["rtdose"].exists()
        has_gtv = paths["gtv"].exists()
        rows.append({
            "patient_id": patient_id,
            "has_rtdose_file": has_rtdose,
            "has_gtv_file": has_gtv,
            "complete": has_rtdose and has_gtv,
            "included": bool(patient_row.get("included", True)),
            "rtdose_path": str(paths["rtdose"]),
            "gtv_path": str(paths["gtv"]),
        })

    manifest = pd.DataFrame(rows)
    n_complete = int(manifest["complete"].sum())
    n_total = len(manifest)

    if verbose:
        print(f"Raw data verification ({data_dir})")
        print(f"  Patients checked : {n_total}")
        print(f"  Complete (RTDOSE + GTV) : {n_complete}")
        print(f"  Missing files    : {n_total - n_complete}")

        if n_complete < n_total:
            missing = manifest[~manifest["complete"]]
            print("\nMissing files:")
            for _, row in missing.head(20).iterrows():
                missing_parts = []
                if not row["has_rtdose_file"]:
                    missing_parts.append("RTDOSE")
                if not row["has_gtv_file"]:
                    missing_parts.append("GTV")
                print(f"  patient {row['patient_id']}: {', '.join(missing_parts)}")
            if len(missing) > 20:
                print(f"  ... and {len(missing) - 20} more")

            nested = find_unorganized_nifti(data_dir)
            if nested:
                print(
                    f"\nFound {len(nested)} RT NIfTI file(s) in nested Aspera paths "
                    f"(not under {{patient_id}}/t0/)."
                )
                print("Flatten them, then re-run verify:")
                print("  python -m src.data.organize_raw_data")
                print("  python -m src.data.verify_raw_data")

    if write_manifest:
        manifest.to_csv(manifest_path, index=False)
        if verbose:
            print(f"\nManifest written to: {manifest_path}")

    return manifest


def is_raw_data_complete(
    cohort_path: Path = COHORT_CSV,
    data_dir: Path = DATA_RAW,
) -> bool:
    """Return True if all included patients have both RTDOSE and GTV files."""
    manifest = verify_raw_data(
        cohort_path=cohort_path,
        data_dir=data_dir,
        included_only=True,
        write_manifest=False,
        verbose=False,
    )
    return bool(manifest["complete"].all())


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Verify RT NIfTI files against cohort.csv")
    parser.add_argument(
        "--cohort",
        type=Path,
        default=COHORT_CSV,
        help="Path to cohort.csv",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DATA_RAW,
        help="Root raw data directory",
    )
    parser.add_argument(
        "--all-patients",
        action="store_true",
        help="Check all patients, not only included ones",
    )
    parser.add_argument(
        "--write-manifest",
        action="store_true",
        help="Write raw_data_manifest.csv to data/processed/",
    )
    parser.add_argument(
        "--organize",
        action="store_true",
        help="Flatten nested Aspera download paths before verifying",
    )
    parser.add_argument(
        "--manifest-path",
        type=Path,
        default=RAW_MANIFEST_CSV,
        help="Output path for manifest CSV",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    try:
        if args.organize:
            from src.data.organize_raw_data import organize_raw_data

            organize_raw_data(data_dir=args.data_dir)
        manifest = verify_raw_data(
            cohort_path=args.cohort,
            data_dir=args.data_dir,
            included_only=not args.all_patients,
            write_manifest=args.write_manifest,
            manifest_path=args.manifest_path,
        )
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)

    if not manifest["complete"].all():
        sys.exit(1)
