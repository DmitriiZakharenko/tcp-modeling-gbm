"""
Unified data setup: download clinical TSVs, build cohort, fetch NIfTI files, verify.

Usage
-----
    python -m src.data.setup_data
    python -m src.data.setup_data --skip-rt
    python -m src.data.setup_data --rt-dry-run
    python -m src.data.setup_data --features --workers 4
"""

import argparse
import subprocess
import sys
from pathlib import Path

from src.config import COHORT_CSV, DATA_RAW
from src.data.cohort_builder import build_cohort
from src.data.download_clinical_data import download_clinical_data
from src.data.download_rt_files import download_rt_files, find_ascp
from src.data.verify_raw_data import is_raw_data_complete, verify_raw_data


def setup_data(
    skip_rt: bool = False,
    rt_dry_run: bool = False,
    force_clinical: bool = False,
    build_features_flag: bool = False,
    feature_workers: int = 1,
    force_features: bool = False,
) -> None:
    """
    Run the full data acquisition pipeline.

    Steps
    -----
    1. Download clinical/metadata TSV files (skip existing unless force).
    2. Build cohort.csv from TSVs.
    3. Download RT NIfTI files via Aspera if incomplete (unless --skip-rt).
    4. Verify raw data completeness and write manifest.
    5. Optionally extract DVH features (--features).
    """
    print("=" * 60)
    print("Step 1/4: Clinical TSV download")
    print("=" * 60)
    download_clinical_data(force=force_clinical)

    print("\n" + "=" * 60)
    print("Step 2/4: Build cohort")
    print("=" * 60)
    build_cohort(output_path=COHORT_CSV)

    if skip_rt:
        print("\n" + "=" * 60)
        print("Step 3/4: RT NIfTI download — SKIPPED")
        print("=" * 60)
    else:
        print("\n" + "=" * 60)
        print("Step 3/4: RT NIfTI download")
        print("=" * 60)
        if is_raw_data_complete(cohort_path=COHORT_CSV, data_dir=DATA_RAW):
            print("All included patients already have RTDOSE + GTV files. Skipping download.")
        else:
            try:
                find_ascp()
                download_rt_files(dry_run=rt_dry_run)
            except FileNotFoundError as exc:
                print(exc)
                print("\nInstall IBM Aspera Connect, then re-run:")
                print("  python -m src.data.setup_data")
                sys.exit(1)
            except (RuntimeError, subprocess.CalledProcessError) as exc:
                print(exc)
                print("\nIf the network blocks Aspera, download manually from:")
                print("  https://www.cancerimagingarchive.net/collection/cfb-gbm/")
                print("Then verify with: python -m src.data.verify_raw_data --write-manifest")
                sys.exit(1)

    print("\n" + "=" * 60)
    print("Step 4/4: Verify raw data")
    print("=" * 60)
    manifest = verify_raw_data(write_manifest=True)

    if not manifest["complete"].all() and not skip_rt and not rt_dry_run:
        print("\nWarning: raw data is incomplete. Re-run setup to resume download:")
        print("  python -m src.data.setup_data")
        sys.exit(1)

    if build_features_flag:
        from src.data.feature_builder import build_features
        from src.data.export_modeling_dataset import export_modeling_dataset

        print("\n" + "=" * 60)
        print("Step 5: Extract DVH features")
        print("=" * 60)
        build_features(workers=feature_workers, force=force_features)

        print("\n" + "=" * 60)
        print("Step 6: Export modeling table (team share)")
        print("=" * 60)
        export_modeling_dataset()

    print("\nData setup complete.")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download and prepare CFB-GBM data")
    parser.add_argument(
        "--skip-rt",
        action="store_true",
        help="Skip RT NIfTI download (clinical TSVs + cohort only)",
    )
    parser.add_argument(
        "--rt-dry-run",
        action="store_true",
        help="Print ascp command without downloading NIfTI files",
    )
    parser.add_argument(
        "--force-clinical",
        action="store_true",
        help="Re-download clinical TSV files",
    )
    parser.add_argument(
        "--features",
        action="store_true",
        help="Extract DVH features after raw data verification",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel workers for feature extraction (with --features)",
    )
    parser.add_argument(
        "--force-features",
        action="store_true",
        help="Reprocess all patients during feature extraction",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    setup_data(
        skip_rt=args.skip_rt,
        rt_dry_run=args.rt_dry_run,
        force_clinical=args.force_clinical,
        build_features_flag=args.features,
        feature_workers=args.workers,
        force_features=args.force_features,
    )
