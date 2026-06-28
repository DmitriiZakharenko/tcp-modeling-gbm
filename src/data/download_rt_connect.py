"""
Download RT NIfTI files via IBM Aspera Connect (same path as the browser plugin).

When direct `ascp` SSH to port 33001 is blocked, the browser plugin still works because it
routes through the locally installed IBM Aspera Connect app over HTTPS (port 443).

Flow:
  1. OAuth bearer token from the TCIA package passcode (HTTPS 443)
  2. Browse package files and collect exact paths for RT NIfTI files
  3. Faspex transfer_spec/download with that file list (token must match paths)
  4. Start transfer through local Connect REST API without modifying paths

Usage
-----
    python -m src.data.download_rt_connect
"""

from __future__ import annotations

import base64
import fnmatch
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
from pathlib import Path
from typing import Any, Optional

from src.config import (
    ASPERA_MAX_RATE_MBPS,
    COHORT_CSV,
    DATA_RAW,
    FASPEX_API_BASE,
    FASPEX_CLIENT_ID,
    FASPEX_PACKAGE_ID,
    FASPEX_PASSCODE,
    FASPEX_REDIRECT_URI,
    GTV_GLOB,
    RTDOSE_GLOB,
)
from src.data.nifti_paths import expected_nifti_paths

CONNECT_RUN_DIRS = (
    Path.home() / "Library/Application Support/Aspera/Aspera Connect/var/run",
    Path.home() / ".aspera/connect/var/run",
)


def build_download_context() -> str:
    """Build the base64 OAuth context Faspex expects for external package downloads."""
    payload = {
        "resource": "packages",
        "type": "external_download_package",
        "id": FASPEX_PACKAGE_ID,
        "passcode": FASPEX_PASSCODE,
        "package_id": FASPEX_PACKAGE_ID,
    }
    return base64.b64encode(json.dumps(payload, separators=(",", ":")).encode()).decode()


def _matches_nifti_basename(basename: str) -> bool:
    """Return True if basename matches RTDOSE or GTV NIfTI patterns."""
    return fnmatch.fnmatch(basename, "*_t0_rtdose.nii.gz") or fnmatch.fnmatch(
        basename, "*_t0_gtv.nii.gz"
    )


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):  # noqa: ARG002
        return None


def _http_json(
    url: str,
    method: str = "GET",
    body: Optional[dict[str, Any]] = None,
    headers: Optional[dict[str, str]] = None,
) -> dict[str, Any]:
    """Perform HTTPS JSON request and return parsed response."""
    data, _ = _http_json_with_headers(url, method=method, body=body, headers=headers)
    return data


def _http_json_with_headers(
    url: str,
    method: str = "GET",
    body: Optional[dict[str, Any]] = None,
    headers: Optional[dict[str, str]] = None,
) -> tuple[dict[str, Any], dict[str, str]]:
    """Perform HTTPS JSON request and return parsed response plus response headers."""
    req_headers = {"Accept": "application/json", "Content-Type": "application/json"}
    if headers:
        req_headers.update(headers)

    payload = None
    if body is not None:
        payload = json.dumps(body).encode("utf-8")

    request = urllib.request.Request(url, data=payload, headers=req_headers, method=method)
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            raw = response.read().decode("utf-8")
            data = json.loads(raw) if raw else {}
            response_headers = {k: v for k, v in response.headers.items()}
            return data, response_headers
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {exc.code} from {url}: {detail}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Network error contacting {url}: {exc.reason}") from exc


def _http_form(url: str, form: dict[str, str]) -> dict[str, Any]:
    """POST application/x-www-form-urlencoded and parse JSON response."""
    payload = urllib.parse.urlencode(form).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=payload,
        headers={
            "Accept": "application/json",
            "Content-Type": "application/x-www-form-urlencoded",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            raw = response.read().decode("utf-8")
            return json.loads(raw) if raw else {}
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {exc.code} from {url}: {detail}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Network error contacting {url}: {exc.reason}") from exc


def _http_authorize_code(context: str) -> str:
    """Exchange a Faspex public-link context for an OAuth authorization code."""
    query = urllib.parse.urlencode(
        {
            "response_type": "code",
            "state": context,
            "client_id": FASPEX_CLIENT_ID,
            "redirect_uri": FASPEX_REDIRECT_URI,
        }
    )
    authorize_url = f"{FASPEX_API_BASE}/auth/authorize_public_link?{query}"
    opener = urllib.request.build_opener(_NoRedirect())
    request = urllib.request.Request(authorize_url, method="GET")
    try:
        opener.open(request, timeout=60)
        raise RuntimeError(f"No redirect from Faspex authorize endpoint: {authorize_url}")
    except urllib.error.HTTPError as exc:
        if exc.code not in (301, 302, 303, 307, 308):
            detail = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"HTTP {exc.code} from {authorize_url}: {detail}") from exc
        location = exc.headers.get("Location")
        if not location:
            raise RuntimeError(f"Missing Location header from {authorize_url}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Network error contacting {authorize_url}: {exc.reason}") from exc

    redirect = urllib.parse.urljoin(authorize_url, location)
    params = urllib.parse.parse_qs(urllib.parse.urlparse(redirect).query)
    if params.get("action_message"):
        raise RuntimeError(params["action_message"][0])
    code = params.get("code", [None])[0]
    if not code:
        raise RuntimeError(f"No authorization code in redirect: {redirect}")
    return code


def fetch_faspex_bearer_token() -> str:
    """Obtain a Faspex bearer token from the TCIA package passcode over HTTPS."""
    code = _http_authorize_code(build_download_context())
    token_url = f"{FASPEX_API_BASE}/auth/token"
    body = {
        "grant_type": "authorization_code",
        "code": code,
        "redirect_uri": FASPEX_REDIRECT_URI,
        "client_id": FASPEX_CLIENT_ID,
    }

    errors: list[str] = []
    for label, request in (
        ("form", lambda: _http_form(token_url, body)),
        ("json", lambda: _http_json(token_url, method="POST", body=body)),
    ):
        try:
            data = request()
            token = data.get("access_token") or data.get("token")
            if token:
                return str(token)
            errors.append(f"{label}: no access_token in response")
        except RuntimeError as exc:
            errors.append(f"{label}: {exc}")

    raise RuntimeError(
        "Could not obtain Faspex bearer token from authorization code.\n"
        + "\n".join(f"  - {e}" for e in errors[:5])
    )


def fetch_package_info_bearer(bearer_token: str) -> dict[str, Any]:
    """Fetch package metadata needed for transfer_spec query parameters."""
    url = f"{FASPEX_API_BASE}/api/v5/packages/{FASPEX_PACKAGE_ID}"
    return _http_json(url, headers={"Authorization": f"Bearer {bearer_token}"})


def _recipient_query_params(package_info: dict[str, Any]) -> dict[str, str]:
    recipients = package_info.get("recipients") or []
    if not recipients:
        return {}
    recipient = recipients[0]
    rtype = recipient.get("recipient_type")
    rid = recipient.get("id")
    if rtype == "workgroup" and rid:
        return {"recipient_workgroup_id": str(rid)}
    if rtype == "user" and rid:
        return {"recipient_user_id": str(rid)}
    return {}


def _browse_package_page(
    bearer_token: str,
    recipient_query: dict[str, str],
    path: str,
    iteration_token: Optional[str] = None,
) -> tuple[dict[str, Any], Optional[str]]:
    """Browse one page of package contents; returns items and next iteration token."""
    query_params = {"per_page": 1000, **recipient_query}
    if iteration_token:
        query_params["iteration_token"] = iteration_token
    query = urllib.parse.urlencode(query_params)
    url = (
        f"{FASPEX_API_BASE}/api/v5/packages/{FASPEX_PACKAGE_ID}/files/received/page"
        f"?{query}"
    )
    data, headers = _http_json_with_headers(
        url,
        method="POST",
        body={"path": path},
        headers={"Authorization": f"Bearer {bearer_token}"},
    )
    next_token = headers.get("X-Aspera-Next-Iteration-Token") or headers.get(
        "x-aspera-next-iteration-token"
    )
    if next_token == "":
        next_token = None
    return data, next_token


def _is_downloadable_item(item: dict[str, Any]) -> bool:
    """Return True for package entries that represent downloadable files."""
    return item.get("type") in ("file", "symbolic_link")


def _browse_all_items(
    bearer_token: str,
    recipient_query: dict[str, str],
    path: str,
) -> list[dict[str, Any]]:
    """Return all items in a package directory, following pagination."""
    items: list[dict[str, Any]] = []
    iteration_token: Optional[str] = None
    while True:
        page, iteration_token = _browse_package_page(
            bearer_token, recipient_query, path, iteration_token
        )
        items.extend(page.get("items") or [])
        if not iteration_token:
            break
    return items


def browse_nifti_files(
    bearer_token: str,
    package_info: dict[str, Any],
) -> list[dict[str, Any]]:
    """
    Collect RT NIfTI files from the CFB-GBM package layout.

    Package structure is ``/CFB-GBM/{patient_id}/t0/*_t0_rtdose.nii.gz`` and
    ``*_t0_gtv.nii.gz``. Faspex issues transfer tokens for a fixed file list,
    so we must resolve exact paths before requesting the transfer spec.
    """
    recipient_query = _recipient_query_params(package_info)
    matched: list[dict[str, Any]] = []

    root_items = _browse_all_items(bearer_token, recipient_query, "/")
    package_dirs = [
        str(item["path"])
        for item in root_items
        if item.get("type") == "directory" and item.get("path")
    ]
    if not package_dirs:
        raise RuntimeError("Package root is empty — cannot locate CFB-GBM directory.")

    patient_dirs: list[str] = []
    for package_dir in package_dirs:
        for item in _browse_all_items(bearer_token, recipient_query, package_dir):
            if item.get("type") == "directory" and item.get("path"):
                patient_dirs.append(str(item["path"]))

    if not patient_dirs:
        raise RuntimeError(f"No patient directories found under {package_dirs[0]!r}.")

    print(f"  Scanning {len(patient_dirs)} patient folders for RT NIfTI files...")
    for index, patient_dir in enumerate(patient_dirs, start=1):
        t0_path = f"{patient_dir.rstrip('/')}/t0"
        for item in _browse_all_items(bearer_token, recipient_query, t0_path):
            if _is_downloadable_item(item) and _matches_nifti_basename(
                str(item.get("basename", ""))
            ):
                matched.append(item)
        if index % 50 == 0 or index == len(patient_dirs):
            print(f"  ... {index}/{len(patient_dirs)} patients, {len(matched)} NIfTI files")

    return matched


def _faspex_api_paths(files: list[dict[str, Any]]) -> list[dict[str, str]]:
    """Convert browse results to PostDownloadTransferSpec path entries."""
    api_paths: list[dict[str, str]] = []
    for item in files:
        item_type = str(item.get("type", "file"))
        if item_type == "symbolic_link":
            item_type = "file"
        api_paths.append(
            {
                "path": str(item["path"]),
                "basename": str(item["basename"]),
                "type": item_type,
            }
        )
    return api_paths


def _file_needed_on_disk(item: dict[str, Any], output_dir: Path) -> bool:
    """Return True if this package file is missing or empty locally."""
    basename = str(item.get("basename", ""))
    match = re.match(r"^(\d+)_t0_(rtdose|gtv)\.nii\.gz$", basename)
    if not match:
        return True
    expected = expected_nifti_paths(match.group(1), output_dir)[match.group(2)]
    return not expected.is_file() or expected.stat().st_size == 0


def _filter_files_needed(
    nifti_files: list[dict[str, Any]],
    output_dir: Path,
) -> list[dict[str, Any]]:
    return [item for item in nifti_files if _file_needed_on_disk(item, output_dir)]


def fetch_transfer_spec_bearer(
    bearer_token: str,
    api_paths: list[dict[str, str]],
    package_info: Optional[dict[str, Any]] = None,
    transfer_type: str = "connect",
) -> dict[str, Any]:
    """
    Fetch a Faspex transfer spec for an explicit file list.

    The server binds the transfer token to these paths; do not replace them later.
    """
    if package_info is None:
        package_info = fetch_package_info_bearer(bearer_token)
    if not api_paths:
        raise RuntimeError("No RT NIfTI files found in package — nothing to download.")

    query = urllib.parse.urlencode(
        {
            "transfer_type": transfer_type,
            "type": "received",
            **_recipient_query_params(package_info),
        }
    )
    url = (
        f"{FASPEX_API_BASE}/api/v5/packages/{FASPEX_PACKAGE_ID}/transfer_spec/download"
        f"?{query}"
    )
    data = _http_json(
        url,
        method="POST",
        body={"paths": api_paths},
        headers={"Authorization": f"Bearer {bearer_token}"},
    )
    if "data" in data and isinstance(data["data"], dict):
        return data["data"]
    return data


def _connect_api_base() -> str:
    for run_dir in CONNECT_RUN_DIRS:
        for name in ("http.uri", "https.uri"):
            uri_file = run_dir / name
            if uri_file.is_file():
                base = uri_file.read_text(encoding="utf-8").splitlines()[0].strip()
                return f"{base.rstrip('/')}/v5/connect"
    raise RuntimeError(
        "IBM Aspera Connect is not running.\n"
        'Start it with: open -a "IBM Aspera Connect"'
    )


def _connect_post(path: str, body: dict[str, Any]) -> dict[str, Any]:
    url = f"{_connect_api_base()}/{path.lstrip('/')}"
    request = urllib.request.Request(
        url,
        data=json.dumps(body).encode("utf-8"),
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
            "Origin": "tcp-modeling-gbm",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            raw = response.read().decode("utf-8")
            return json.loads(raw) if raw else {}
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"Connect API HTTP {exc.code} from {url}: {detail}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Connect API error at {url}: {exc.reason}") from exc


def _connect_settings() -> dict[str, str]:
    return {"app_id": str(uuid.uuid4())}


def _extract_xfer_id(start_response: dict[str, Any]) -> str:
    try:
        return start_response["transfer_specs"][0]["transfer_spec"]["tags"]["aspera"]["xfer_id"]
    except (KeyError, IndexError, TypeError) as exc:
        raise RuntimeError(f"Unexpected Connect start response: {start_response}") from exc


def _format_bytes(num_bytes: int) -> str:
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"


def _flatten_transfer_spec_paths(transfer_spec: dict[str, Any]) -> None:
    """
    Map remote package paths to ``{patient_id}/t0/{basename}`` under destination_root.

    Connect otherwise preserves the full Faspex package tree and verify cannot find files.
    """
    paths = transfer_spec.get("paths")
    if not isinstance(paths, list):
        return

    flat_paths: list[dict[str, Any]] = []
    for entry in paths:
        if not isinstance(entry, dict):
            continue
        source = str(entry.get("source") or entry.get("path") or "")
        basename = source.rsplit("/", 1)[-1]
        match = re.match(r"^(\d+)_t0_(rtdose|gtv)\.nii\.gz$", basename)
        destination = f"{match.group(1)}/t0/{basename}" if match else basename
        flat_paths.append({**entry, "source": source, "destination": destination})
    transfer_spec["paths"] = flat_paths


def start_connect_transfer(
    spec: dict[str, Any],
    output_dir: Path,
    max_rate_mbps: int,
) -> tuple[str, dict[str, str]]:
    """Start the transfer via Connect REST API using the Faspex-issued spec unchanged."""
    settings = _connect_settings()
    transfer_spec = dict(spec)
    transfer_spec.pop("authentication", None)
    transfer_spec["destination_root"] = str(output_dir.resolve())
    _flatten_transfer_spec_paths(transfer_spec)
    transfer_spec["resume_policy"] = "sparse_csum"
    transfer_spec["target_rate_kbps"] = max_rate_mbps * 1000
    if transfer_spec.get("token"):
        transfer_spec["authentication"] = "token"

    body = {
        "aspera_connect_settings": {
            **settings,
            "request_id": str(uuid.uuid4()),
            "allow_dialogs": True,
        },
        "transfer_specs": [{"transfer_spec": transfer_spec}],
    }
    response = _connect_post("transfers/start", body)
    return _extract_xfer_id(response), settings


def wait_connect_transfer(xfer_id: str, settings: dict[str, str], poll_seconds: float = 2.0) -> None:
    """Poll Connect until the transfer completes or fails."""
    body = {"aspera_connect_settings": settings}
    last_written = -1
    while True:
        info = _connect_post(f"transfers/info/{xfer_id}", body)
        transfer = info.get("transfer_info")
        if not isinstance(transfer, dict):
            time.sleep(poll_seconds)
            continue

        status = transfer.get("status", "")
        if status in ("initiating", "queued"):
            print(f"Transfer {status}...")
        elif status == "running":
            written = int(transfer.get("bytes_written") or 0)
            expected = int(transfer.get("bytes_expected") or 0)
            if written != last_written:
                if expected > 0:
                    pct = 100.0 * written / expected
                    print(
                        f"Downloading: {_format_bytes(written)} / {_format_bytes(expected)} "
                        f"({pct:.1f}%)"
                    )
                else:
                    print(f"Downloading: {_format_bytes(written)}")
                last_written = written
        elif status == "completed":
            print("Transfer completed.")
            return
        elif status == "failed":
            raise RuntimeError(transfer.get("error_desc") or "Connect transfer failed")
        elif status == "cancelled":
            raise RuntimeError("Connect transfer cancelled")
        else:
            raise RuntimeError(f"Unknown Connect transfer status: {status}")

        time.sleep(poll_seconds)


def download_rt_connect(
    output_dir: Path = DATA_RAW,
    max_rate_mbps: int = ASPERA_MAX_RATE_MBPS,
    dry_run: bool = False,
    missing_only: bool = False,
    cohort_path: Path = COHORT_CSV,
) -> None:
    """Download filtered RT NIfTI files via Faspex HTTPS + local IBM Aspera Connect."""
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Using IBM Aspera Connect (HTTPS 443 + local Connect API)")
    print(f"Destination: {output_dir}")
    print(f"Files: {RTDOSE_GLOB}, {GTV_GLOB}")
    if missing_only:
        print("Mode: missing files only (skip paths already on disk)\n")
    else:
        print()

    if dry_run:
        print("Would:")
        print("  1. Obtain OAuth bearer token from TCIA passcode")
        print("  2. Browse package and collect RT NIfTI file paths")
        if missing_only:
            print("  3. Keep only files missing under destination")
        print("  4. Fetch Faspex transfer_spec/download for those paths")
        print("  5. Start Connect transfer")
        print(f"Connect API base: {_connect_api_base()}")
        return

    print("Requesting OAuth bearer token...")
    token = fetch_faspex_bearer_token()

    print("Fetching package metadata...")
    package_info = fetch_package_info_bearer(token)

    print("Browsing package for RT NIfTI files...")
    nifti_files = browse_nifti_files(token, package_info)
    print(f"Found {len(nifti_files)} files matching {RTDOSE_GLOB} and {GTV_GLOB}")

    if missing_only:
        if cohort_path.is_file():
            import pandas as pd

            cohort = pd.read_csv(cohort_path)
            included_ids = {
                str(pid) for pid in cohort.loc[cohort["included"], "patient_id"].astype(str)
            }
            nifti_files = [
                item
                for item in nifti_files
                if str(item.get("basename", "")).split("_", 1)[0] in included_ids
            ]
        nifti_files = _filter_files_needed(nifti_files, output_dir)
        print(f"Missing locally: {len(nifti_files)} file(s) to download")
        if not nifti_files:
            print("All requested files already present.")
            return

    print("Fetching transfer specification (transfer_type=connect)...")
    spec = fetch_transfer_spec_bearer(
        token,
        _faspex_api_paths(nifti_files),
        package_info=package_info,
        transfer_type="connect",
    )

    print("Starting Connect transfer...")
    xfer_id, settings = start_connect_transfer(spec, output_dir, max_rate_mbps)
    print(f"Transfer id: {xfer_id}")

    wait_connect_transfer(xfer_id, settings)
    print(f"\nFiles saved under: {output_dir}")


def _parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description="Download RT NIfTI via IBM Aspera Connect (browser plugin path)"
    )
    parser.add_argument("--output-dir", type=Path, default=DATA_RAW)
    parser.add_argument("--max-rate", type=int, default=ASPERA_MAX_RATE_MBPS, metavar="MBPS")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--missing-only",
        action="store_true",
        help="Download only RT NIfTI files missing under output-dir (included cohort)",
    )
    parser.add_argument("--cohort", type=Path, default=COHORT_CSV)
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    try:
        download_rt_connect(
            output_dir=args.output_dir,
            max_rate_mbps=args.max_rate,
            dry_run=args.dry_run,
            missing_only=args.missing_only,
            cohort_path=args.cohort,
        )
    except (FileNotFoundError, RuntimeError) as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)
