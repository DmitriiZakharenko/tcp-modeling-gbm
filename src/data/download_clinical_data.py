"""
Download CFB-GBM clinical and metadata TSV files from TCIA.

Skips files that already exist locally unless --force is passed.

Usage
-----
    python -m src.data.download_clinical_data
    python -m src.data.download_clinical_data --force
"""

import argparse
import sys
import urllib.error
import urllib.request
from pathlib import Path

from src.config import CLINICAL_FILES, DATA_PROCESSED


def download_file(url: str, dest: Path, force: bool = False) -> str:
    """
    Download a single file if it does not already exist.

    Returns
    -------
    str
        One of: "downloaded", "skipped", "failed".
    """
    if dest.exists() and not force:
        print(f"  skip  {dest.name} (already exists)")
        return "skipped"

    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  fetch {dest.name} ...", end=" ", flush=True)
    try:
        urllib.request.urlretrieve(url, dest)
    except urllib.error.URLError as exc:
        print(f"FAILED\n    {exc}")
        return "failed"

    size_kb = dest.stat().st_size / 1024
    print(f"OK ({size_kb:.1f} KB)")
    return "downloaded"


def download_clinical_data(
    output_dir: Path = DATA_PROCESSED,
    force: bool = False,
) -> dict[str, int]:
    """
    Download all clinical/metadata TSV files to data/processed/.

    Parameters
    ----------
    output_dir : Path
        Destination directory for TSV files.
    force : bool
        Re-download even if files already exist.

    Returns
    -------
    dict
        Counts keyed by status: downloaded, skipped, failed.
    """
    print(f"Clinical TSV download → {output_dir}\n")
    counts = {"downloaded": 0, "skipped": 0, "failed": 0}

    for clinical_file in CLINICAL_FILES:
        dest = output_dir / clinical_file.path.name
        status = download_file(clinical_file.url, dest, force=force)
        counts[status] += 1
        if status == "failed" and clinical_file.required:
            pass  # counted; raise below if any required failed

    print(
        f"\nDone: {counts['downloaded']} downloaded, "
        f"{counts['skipped']} skipped, {counts['failed']} failed"
    )
    if counts["failed"]:
        failed_required = sum(
            1 for f in CLINICAL_FILES if f.required and not (output_dir / f.path.name).exists()
        )
        if failed_required:
            raise RuntimeError(f"{failed_required} required file(s) failed to download")

    return counts


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download CFB-GBM clinical TSV files from TCIA")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DATA_PROCESSED,
        help="Directory for downloaded TSV files",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download files even if they already exist",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    try:
        download_clinical_data(output_dir=args.output_dir, force=args.force)
    except RuntimeError:
        sys.exit(1)
