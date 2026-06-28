"""
Move RT NIfTI files from nested Aspera/Faspex paths into the expected layout.

Connect preserves the package directory tree under destination_root, e.g.::

    data/raw/CFB-GBM version 1 - .../CFB-GBM/6/t0/6_t0_gtv.nii.gz

The pipeline expects::

    data/raw/6/t0/6_t0_gtv.nii.gz

Usage
-----
    python -m src.data.organize_raw_data
    python -m src.data.organize_raw_data --dry-run
    python -m src.data.organize_raw_data --recover-partials
    python -m src.data.organize_raw_data --import-from ~/Downloads/.../data
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

from src.config import DATA_RAW
from src.data.nifti_paths import expected_nifti_paths

_NIFTI_RE = re.compile(r"^(\d+)_t0_(rtdose|gtv)\.nii\.gz$")
_ASPERA_STAGING_SUFFIXES = (".partial", ".aspera-ckpt")


def _expected_path(data_dir: Path, basename: str) -> Path | None:
    match = _NIFTI_RE.match(basename)
    if not match:
        return None
    patient_id = match.group(1)
    return expected_nifti_paths(patient_id, data_dir)[match.group(2)]


def find_unorganized_nifti(data_dir: Path = DATA_RAW) -> list[Path]:
    """Return NIfTI files not already at the expected patient/t0 location."""
    if not data_dir.is_dir():
        return []

    unorganized: list[Path] = []
    for path in data_dir.rglob("*_t0_*.nii.gz"):
        expected = _expected_path(data_dir, path.name)
        if expected is None:
            continue
        if path.resolve() != expected.resolve():
            unorganized.append(path)
    return unorganized


def find_misplaced_aspera_staging(data_dir: Path = DATA_RAW) -> list[Path]:
    """Return Aspera resume files (.partial, .aspera-ckpt) sitting in data_dir root."""
    if not data_dir.is_dir():
        return []

    staging: list[Path] = []
    for path in data_dir.iterdir():
        if not path.is_file():
            continue
        for suffix in _ASPERA_STAGING_SUFFIXES:
            if path.name.endswith(suffix):
                basename = path.name[: -len(suffix)]
                if _expected_path(data_dir, basename) is not None:
                    staging.append(path)
                break
    return staging


def recover_aspera_staging(
    data_dir: Path = DATA_RAW,
    dry_run: bool = False,
    verbose: bool = True,
) -> dict[str, int]:
    """
    Move Connect resume files from ``data/raw/`` root into ``{patient_id}/t0/``.

    Connect sometimes writes ``220_t0_rtdose.nii.gz.partial`` at the destination
    root even when the transfer spec destination is ``220/t0/220_t0_rtdose.nii.gz``.
    Resume only works when staging files sit next to the final path.
    """
    stats = {"moved": 0, "skipped": 0}

    for path in find_misplaced_aspera_staging(data_dir):
        for suffix in _ASPERA_STAGING_SUFFIXES:
            if path.name.endswith(suffix):
                basename = path.name[: -len(suffix)]
                break
        else:
            continue

        expected = _expected_path(data_dir, basename)
        if expected is None:
            continue
        target = expected.with_name(expected.name + path.name[len(basename) :])

        if target.exists():
            stats["skipped"] += 1
            if verbose:
                print(f"  skip (already exists): {target}")
            continue

        if verbose:
            print(f"  {path.name} -> {target.relative_to(data_dir)}")

        if not dry_run:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(path), str(target))
        stats["moved"] += 1

    if verbose:
        print(
            f"\nRecover staging summary: moved={stats['moved']}, skipped={stats['skipped']}"
        )

    return stats


def import_nifti_from(
    source_dir: Path,
    data_dir: Path = DATA_RAW,
    dry_run: bool = False,
    move: bool = False,
    verbose: bool = True,
) -> dict[str, int]:
    """
    Copy or move RT NIfTI files from an external tree (e.g. Downloads) into ``{patient_id}/t0/``.

    Skips ``.partial`` files. Does not overwrite a larger existing destination file.
    Use ``move=True`` to avoid doubling disk use (same volume: near-instant rename).
    """
    if not source_dir.is_dir():
        raise FileNotFoundError(f"Source directory not found: {source_dir}")

    stats = {"copied": 0, "skipped_existing": 0, "ignored": 0, "moved": 0}

    # Same basename may appear in multiple nested download paths — keep the largest.
    best_sources: dict[str, Path] = {}
    for path in source_dir.rglob("*_t0_*.nii.gz"):
        if path.name.endswith(".partial") or ".partial" in path.name:
            stats["ignored"] += 1
            continue
        if not _NIFTI_RE.match(path.name):
            stats["ignored"] += 1
            continue
        previous = best_sources.get(path.name)
        if previous is None or path.stat().st_size > previous.stat().st_size:
            best_sources[path.name] = path

    for path in best_sources.values():
        expected = _expected_path(data_dir, path.name)
        if expected is None:
            stats["ignored"] += 1
            continue

        if expected.exists() and expected.stat().st_size >= path.stat().st_size:
            stats["skipped_existing"] += 1
            continue

        if verbose:
            print(f"  {path} -> {expected}")

        if not dry_run:
            expected.parent.mkdir(parents=True, exist_ok=True)
            if move:
                shutil.move(str(path), str(expected))
                stats["moved"] += 1
            else:
                shutil.copy2(path, expected)
                stats["copied"] += 1
        elif move:
            stats["moved"] += 1
        else:
            stats["copied"] += 1

    if verbose:
        action = "moved" if move else "copied"
        print(
            f"\nImport summary: {action}={stats[action]}, "
            f"skipped_existing={stats['skipped_existing']}, ignored={stats['ignored']}"
        )

    return stats


def organize_raw_data(
    data_dir: Path = DATA_RAW,
    dry_run: bool = False,
    verbose: bool = True,
) -> dict[str, int]:
    """
    Move RT NIfTI files into ``{data_dir}/{patient_id}/t0/``.

    Returns counts: moved, skipped_existing, removed_duplicates.
    """
    stats = {"moved": 0, "skipped_existing": 0, "removed_duplicates": 0}

    for path in find_unorganized_nifti(data_dir):
        expected = _expected_path(data_dir, path.name)
        if expected is None:
            continue

        if expected.exists():
            if expected.stat().st_size >= path.stat().st_size:
                stats["removed_duplicates"] += 1
                if verbose:
                    print(f"  remove duplicate: {path}")
                if not dry_run:
                    path.unlink()
                continue

            stats["skipped_existing"] += 1
            if verbose:
                print(f"  skip (destination larger): {expected}")
            continue

        if verbose:
            print(f"  {path} -> {expected}")

        if not dry_run:
            expected.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(path), str(expected))
        stats["moved"] += 1

    if verbose:
        print(
            f"\nOrganize summary: moved={stats['moved']}, "
            f"skipped_existing={stats['skipped_existing']}, "
            f"removed_duplicates={stats['removed_duplicates']}"
        )

    return stats


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Flatten nested Aspera downloads into data/raw/{patient_id}/t0/"
    )
    parser.add_argument("--data-dir", type=Path, default=DATA_RAW)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--recover-partials",
        action="store_true",
        help="Move .partial/.aspera-ckpt files from data/raw/ root into patient/t0/",
    )
    parser.add_argument(
        "--import-from",
        type=Path,
        metavar="DIR",
        help="Copy RT NIfTI files from an external directory (e.g. Downloads mirror)",
    )
    parser.add_argument(
        "--move",
        action="store_true",
        help="With --import-from: move files instead of copy (saves disk space)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()

    if args.import_from:
        print(f"Importing RT NIfTI from {args.import_from} -> {args.data_dir}")
        import_nifti_from(
            source_dir=args.import_from,
            data_dir=args.data_dir,
            dry_run=args.dry_run,
            move=args.move,
        )

    if args.recover_partials:
        staging = find_misplaced_aspera_staging(args.data_dir)
        if not staging:
            print(f"No misplaced Aspera staging files under {args.data_dir}")
        else:
            print(f"Found {len(staging)} staging file(s) under {args.data_dir}")
            recover_aspera_staging(data_dir=args.data_dir, dry_run=args.dry_run)

    unorganized = find_unorganized_nifti(args.data_dir)
    if unorganized:
        print(f"Found {len(unorganized)} nested RT NIfTI file(s) under {args.data_dir}")
        organize_raw_data(data_dir=args.data_dir, dry_run=args.dry_run)
    elif not args.recover_partials and not args.import_from:
        print(f"No nested RT NIfTI files found under {args.data_dir}")
