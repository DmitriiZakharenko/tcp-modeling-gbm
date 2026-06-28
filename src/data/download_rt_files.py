"""
Selective download of RTDOSE and GTV NIfTI files from the CFB-GBM TCIA package.

Uses the ascp CLI (ships with IBM Aspera Connect) with include/exclude patterns
to download only *_t0_rtdose.nii.gz and *_t0_gtv.nii.gz files (~1–1.5 GB)
instead of the full 208 GB archive.

Usage
-----
    python -m src.data.download_rt_files
    python -m src.data.download_rt_files --dry-run
    python -m src.data.download_rt_files --port 33001
    python -m src.data.download_rt_files --max-rate 500 --no-resume

Requirements
------------
    IBM Aspera Connect must be installed:
    https://www.ibm.com/products/aspera/downloads

    Firewall must allow outbound TCP and UDP on port 33001 (Aspera default).
"""

import argparse
import os
import socket
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Optional

from src.config import (
    ASPERA_MAX_RATE_MBPS,
    DATA_RAW,
    FASPEX_FALLBACK_PORTS,
    FASPEX_HOST,
    FASPEX_PACKAGE_PATH,
    FASPEX_PASSCODE,
    FASPEX_PORT,
    FASPEX_UDP_PORT,
    FASPEX_USER,
    GTV_GLOB,
    RTDOSE_GLOB,
)

# --- ascp binary path (macOS default; adjust if on Linux/Windows) ---
ASCP_PATHS = [
    Path("/Applications/IBM Aspera.app/Contents/Resources/transferd/bin/ascp"),
    Path("/Applications/Aspera Connect.app/Contents/Resources/ascp"),
    Path.home() / ".aspera/connect/bin/ascp",
    Path("C:/Users") / Path.home().name / "AppData/Local/Programs/Aspera/Aspera Connect/bin/ascp.exe",
]

SSH_CONNECT_ERROR_MARKERS = (
    "failed to open tcp connection for ssh",
    "connection refused",
    "connection timed out",
    "no route to host",
    "network is unreachable",
)


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


def check_tcp_port(host: str, port: int, timeout: float = 5.0) -> bool:
    """Return True if a TCP connection to host:port succeeds."""
    try:
        with socket.create_connection((host, port), timeout=timeout):
            return True
    except OSError:
        return False


def resolve_ssh_port(preferred_port: Optional[int] = None) -> int:
    """Return SSH port for TCIA Faspex (33001 by default)."""
    if preferred_port is not None:
        return preferred_port
    return FASPEX_PORT


def build_ascp_command(
    ascp: Path,
    output_dir: Path,
    ssh_port: int = FASPEX_PORT,
    udp_port: int = FASPEX_UDP_PORT,
    max_rate_mbps: int = ASPERA_MAX_RATE_MBPS,
    resume: bool = True,
) -> list[str]:
    """Build the ascp command for selective RT NIfTI download."""
    cmd = [
        str(ascp),
        "-T",                        # Disable encryption (faster, data is public)
        "-l", f"{max_rate_mbps}m",   # Max transfer rate
        "-P", str(ssh_port),         # SSH / session initiation port
        "-O", str(udp_port),         # FASP UDP data port
        "--mode=recv",
        "--host", FASPEX_HOST,
        "-N", RTDOSE_GLOB.replace("**/", "*"),  # Include RTDOSE files only
        "-N", GTV_GLOB.replace("**/", "*"),     # Include GTV segmentation files only
        "-E", "*",                   # Exclude everything else
        f"{FASPEX_USER}@{FASPEX_HOST}:{FASPEX_PACKAGE_PATH}",
        str(output_dir),
    ]
    if resume:
        cmd[1:1] = ["-k", "1"]  # Resume partially transferred files
    return cmd


def build_ascp_env() -> dict[str, str]:
    """Environment for TCIA Faspex anonymous package download."""
    env = os.environ.copy()
    # TCIA embeds the Faspex passcode as the transfer password (not -W token).
    env["ASPERA_SCP_PASS"] = FASPEX_PASSCODE
    return env


def _is_ssh_connect_error(stderr: str) -> bool:
    text = stderr.lower()
    return any(marker in text for marker in SSH_CONNECT_ERROR_MARKERS)


def _run_ascp(cmd: list[str], env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        check=False,
        env=env,
        capture_output=True,
        text=True,
    )


def _print_ascp_failure(result: subprocess.CompletedProcess[str]) -> None:
    if result.stdout:
        print(result.stdout, end="" if result.stdout.endswith("\n") else "\n")
    if result.stderr:
        print(result.stderr, end="" if result.stderr.endswith("\n") else "\n", file=sys.stderr)


def _raise_connectivity_help(host: str, ports: Iterable[int]) -> None:
    ports_str = ", ".join(str(p) for p in ports)
    raise RuntimeError(
        f"Cannot reach {host} on TCP port(s) {ports_str}.\n\n"
        "TCIA downloads require IBM Aspera with outbound access to:\n"
        "  - TCP 33001 (SSH session)\n"
        "  - UDP 33001 (FASP data transfer)\n\n"
        "Try:\n"
        "  1. Install IBM Aspera Connect: https://www.ibm.com/products/aspera/downloads\n"
        "  2. Test from browser: https://www.cancerimagingarchive.net/collection/cfb-gbm/\n"
        "  3. If on university/corporate network, ask IT to allow TCP+UDP 33001\n"
        "  4. Retry with explicit port: python -m src.data.download_rt_files --port 33001\n"
        "  5. Download manually via browser Aspera dialog → data/raw/"
    )


def download_rt_files(
    output_dir: Path = DATA_RAW,
    max_rate_mbps: int = ASPERA_MAX_RATE_MBPS,
    ssh_port: Optional[int] = None,
    udp_port: int = FASPEX_UDP_PORT,
    resume: bool = True,
    dry_run: bool = False,
    auto_port: bool = True,
    prefer_connect: bool = True,
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
    ssh_port : int, optional
        SSH port for ascp (-P). Defaults to first reachable port in
        FASPEX_FALLBACK_PORTS, usually 33001 for TCIA.
    udp_port : int
        FASP UDP port for ascp (-O). Default 33001.
    resume : bool
        If True, pass -k 1 to ascp to resume interrupted transfers.
    dry_run : bool
        If True, print the ascp command without executing it.
    auto_port : bool
        If True and ssh_port is None, try fallback ports after SSH failures.
    prefer_connect : bool
        If True, try IBM Aspera Desktop via HTTPS first (same path as browser plugin).

    Raises
    ------
    FileNotFoundError
        If the ascp binary is not found.
    RuntimeError
        If ascp cannot connect to the Faspex server.
    subprocess.CalledProcessError
        If ascp exits with a non-zero return code for other reasons.
    """
    if prefer_connect:
        try:
            from src.data.download_rt_connect import download_rt_connect

            print("Attempting download via IBM Aspera Desktop (HTTPS / Connect agent)...")
            download_rt_connect(
                output_dir=output_dir,
                max_rate_mbps=max_rate_mbps,
                dry_run=dry_run,
            )
            return
        except (FileNotFoundError, RuntimeError, subprocess.CalledProcessError) as exc:
            print(f"Connect/Desktop download unavailable: {exc}")
            print("Falling back to direct ascp (requires TCP 33001)...\n")

    ascp = find_ascp()
    output_dir.mkdir(parents=True, exist_ok=True)

    if ssh_port is None and not check_tcp_port(FASPEX_HOST, FASPEX_PORT):
        print(
            f"Warning: TCP preflight to {FASPEX_HOST}:{FASPEX_PORT} failed. "
            "Aspera may still work, but your network might block TCIA downloads."
        )

    resolved_port = resolve_ssh_port(ssh_port)
    env = build_ascp_env()
    ports_to_try = [resolved_port]
    if auto_port and ssh_port is None:
        ports_to_try = list(dict.fromkeys(FASPEX_FALLBACK_PORTS))

    last_result: Optional[subprocess.CompletedProcess[str]] = None

    for attempt_port in ports_to_try:
        cmd = build_ascp_command(
            ascp, output_dir, attempt_port, udp_port, max_rate_mbps, resume
        )

        print("ascp command:")
        print(" \\\n  ".join(cmd))
        print()

        if dry_run:
            print("Dry run — not executing.")
            return

        print(f"Downloading RTDOSE and GTV files to: {output_dir}")
        print(f"SSH port: {attempt_port} | FASP UDP port: {udp_port}")
        if resume:
            print("Resume enabled (-k 1): partially downloaded files will continue.")
        print("Expected size: ~1–1.5 GB\n")

        result = _run_ascp(cmd, env)
        if result.returncode == 0:
            print(f"\nDownload complete. Return code: {result.returncode}")
            return

        last_result = result
        combined = f"{result.stdout}\n{result.stderr}"
        if auto_port and ssh_port is None and _is_ssh_connect_error(combined):
            print(f"SSH connection failed on port {attempt_port}, trying next port...\n")
            continue

        _print_ascp_failure(result)
        raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)

    if last_result is not None:
        _print_ascp_failure(last_result)
        if _is_ssh_connect_error(f"{last_result.stdout}\n{last_result.stderr}"):
            _raise_connectivity_help(FASPEX_HOST, ports_to_try)
        raise subprocess.CalledProcessError(
            last_result.returncode,
            last_result.args,
            last_result.stdout,
            last_result.stderr,
        )

    _raise_connectivity_help(FASPEX_HOST, ports_to_try)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download RTDOSE and GTV NIfTI files from CFB-GBM via Aspera"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DATA_RAW,
        help="Local directory for NIfTI files",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=None,
        metavar="PORT",
        help=f"SSH port for ascp (-P). Default: auto-detect ({FASPEX_PORT} for TCIA)",
    )
    parser.add_argument(
        "--udp-port",
        type=int,
        default=FASPEX_UDP_PORT,
        metavar="PORT",
        help="FASP UDP port for ascp (-O)",
    )
    parser.add_argument(
        "--max-rate",
        type=int,
        default=ASPERA_MAX_RATE_MBPS,
        metavar="MBPS",
        help="Maximum transfer rate in Mbps",
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="Disable ascp resume (-k 1)",
    )
    parser.add_argument(
        "--no-auto-port",
        action="store_true",
        help="Do not retry alternate SSH ports on connection failure",
    )
    parser.add_argument(
        "--ascp-only",
        action="store_true",
        help="Skip Connect/Desktop and use direct ascp only",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the ascp command without executing",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    try:
        download_rt_files(
            output_dir=args.output_dir,
            max_rate_mbps=args.max_rate,
            ssh_port=args.port,
            udp_port=args.udp_port,
            resume=not args.no_resume,
            dry_run=args.dry_run,
            auto_port=not args.no_auto_port,
            prefer_connect=not args.ascp_only,
        )
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)
    except (RuntimeError, subprocess.CalledProcessError) as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)
