"""
Selective download of RTDOSE and GTV NIfTI files from the CFB-GBM TCIA package.

Uses the ascp CLI (ships with IBM Aspera Connect) with include/exclude patterns
to download only *_t0_rtdose.nii.gz and *_t0_gtv.nii.gz files (~1–1.5 GB)
instead of the full 208 GB archive.

Usage
-----
    python -m src.data.download_rt_files

Requirements
------------
    IBM Aspera Connect must be installed:
    https://www.ibm.com/products/aspera/downloads

    The ascp binary is located at:
    macOS : /Applications/Aspera Connect.app/Contents/Resources/ascp
    Linux : ~/.aspera/connect/bin/ascp
    Windows: C:/Users/<user>/AppData/Local/Programs/Aspera/Aspera Connect/bin/ascp.exe
"""

import subprocess
import sys
from pathlib import Path

from src.config import DATA_RAW

# --- Aspera / TCIA package credentials (public, embedded in the TCIA download URL) ---
FASPEX_HOST = "faspex.cancerimagingarchive.net"
FASPEX_PACKAGE_PATH = "/packages/1196/CFB-GBM"
FASPEX_PASSCODE = "77a6f7a7258b1a4ef9bf2d13fcd3f1ec87b9c313"
FASPEX_USER = "anon"

# --- ascp binary path (macOS default; adjust if on Linux/Windows) ---
ASCP_PATHS = [
    Path("/Applications/IBM Aspera.app/Contents/Resources/transferd/bin/ascp"),
    Path("/Applications/Aspera Connect.app/Contents/Resources/ascp"),
    Path.home() / ".aspera/connect/bin/ascp",
    Path("C:/Users") / Path.home().name / "AppData/Local/Programs/Aspera/Aspera Connect/bin/ascp.exe",
]


def find_ascp() -> Path:
    """
    Locate the ascp binary from Aspera Connect installation.

    Returns
    -------
    Path
        Absolute path to the ascp binary.

    Raises
    ------
    FileNotFoundError
        If ascp is not found at any known location.
    """
    for path in ASCP_PATHS:
        if path.exists():
            return path
    raise FileNotFoundError(
        "ascp binary not found. Ensure IBM Aspera Connect is installed.\n"
        "Download: https://www.ibm.com/products/aspera/downloads\n"
        f"Searched: {[str(p) for p in ASCP_PATHS]}"
    )


def download_rt_files(
    output_dir: Path = DATA_RAW,
    max_rate_mbps: int = 200,
    dry_run: bool = False,
) -> None:
    """
    Download only RTDOSE and GTV NIfTI files for all patients.

    Downloads files matching:
        *_t0_rtdose.nii.gz
        *_t0_gtv.nii.gz

    All other files (MRI sequences, CT) are excluded.

    Parameters
    ----------
    output_dir : Path
        Local directory to download files into. Created if it does not exist.
    max_rate_mbps : int
        Maximum transfer rate in Mbps (default 200).
    dry_run : bool
        If True, print the ascp command without executing it.

    Raises
    ------
    FileNotFoundError
        If the ascp binary is not found.
    subprocess.CalledProcessError
        If ascp exits with a non-zero return code.
    """
    ascp = find_ascp()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Password passed via environment variable, not command-line flag
    import os
    env = os.environ.copy()
    env["ASPERA_SCP_PASS"] = FASPEX_PASSCODE

    cmd = [
        str(ascp),
        "-T",                        # Disable encryption (faster, data is public)
        "-l", f"{max_rate_mbps}m",   # Max transfer rate
        "-P", "33001",               # Faspex SSH port
        "--mode=recv",
        "-N", "*_t0_rtdose.nii.gz",  # Include RTDOSE files only
        "-N", "*_t0_gtv.nii.gz",     # Include GTV segmentation files only
        "-E", "*",                   # Exclude everything else
        f"{FASPEX_USER}@{FASPEX_HOST}:{FASPEX_PACKAGE_PATH}",
        str(output_dir),
    ]

    print("ascp command:")
    print(" \\\n  ".join(cmd))
    print()

    if dry_run:
        print("Dry run — not executing.")
        return

    print(f"Downloading RTDOSE and GTV files to: {output_dir}")
    print("Expected size: ~1–1.5 GB\n")

    result = subprocess.run(cmd, check=True, env=env)
    print(f"\nDownload complete. Return code: {result.returncode}")


if __name__ == "__main__":
    dry_run = "--dry-run" in sys.argv
    download_rt_files(dry_run=dry_run)
