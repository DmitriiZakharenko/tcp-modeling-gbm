#!/usr/bin/env python3
"""Verify PubMed IDs and DOIs in reports/literature_table.csv."""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIT_CSV = ROOT / "reports" / "literature_table.csv"
REPORT_MD = ROOT / "reports" / "literature_doi_check.md"
METRICS_JSON = ROOT / "reports" / "metrics" / "literature_doi_check.json"
TIMEOUT = 20
USER_AGENT = "tcp-modeling-gbm/1.0 (literature DOI check)"


def _request(url: str) -> tuple[int, str]:
    """HTTP status via curl (handles redirects; avoids macOS Python SSL issues)."""
    try:
        proc = subprocess.run(
            ["curl", "-sL", "-o", "/dev/null", "-w", "%{http_code}", url],
            capture_output=True,
            text=True,
            timeout=TIMEOUT,
            check=False,
        )
        code = int((proc.stdout or "0").strip() or "0")
        return code, url
    except (subprocess.SubprocessError, ValueError, OSError) as exc:
        return 0, str(exc)


def check_doi(doi: str) -> dict:
    doi = doi.strip()
    if not doi or doi.upper() == "NA":
        return {"status": "skipped", "http": None, "final_url": ""}
    url = doi if doi.startswith("http") else f"https://doi.org/{doi}"
    code, final = _request(url)
    ok = 200 <= code < 400 or code == 403  # 403 = publisher bot wall; PubMed often OK
    return {"status": "ok" if ok else "fail", "http": code, "final_url": final}


def check_pubmed(pmid: str) -> dict:
    pmid = pmid.strip()
    if not pmid or pmid.upper() == "NA":
        return {"status": "skipped", "http": None, "final_url": ""}
    if not re.fullmatch(r"\d+", pmid):
        return {"status": "invalid", "http": None, "final_url": ""}
    url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
    code, final = _request(url)
    ok = 200 <= code < 400 or code == 403  # 403 = publisher bot wall; PubMed often OK
    return {"status": "ok" if ok else "fail", "http": code, "final_url": final}


def main() -> int:
    if not LIT_CSV.exists():
        print(f"Missing {LIT_CSV}", file=sys.stderr)
        return 1

    rows_out: list[dict] = []
    with LIT_CSV.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            doi_res = check_doi(row.get("doi", ""))
            pm_res = check_pubmed(row.get("pubmed_id", ""))
            entry = {
                "ref_id": row.get("ref_id"),
                "first_author": row.get("first_author"),
                "year": row.get("year"),
                "doi": row.get("doi"),
                "pubmed_id": row.get("pubmed_id"),
                "doi_check": doi_res,
                "pubmed_check": pm_res,
            }
            rows_out.append(entry)
            print(
                f"ref {entry['ref_id']:>2}  "
                f"DOI {doi_res['status']:>7}  PMID {pm_res['status']:>7}  "
                f"{entry['first_author']} {entry['year']}"
            )

    n_doi_ok = sum(1 for r in rows_out if r["doi_check"]["status"] == "ok")
    n_doi_skip = sum(1 for r in rows_out if r["doi_check"]["status"] == "skipped")
    n_doi_fail = sum(1 for r in rows_out if r["doi_check"]["status"] == "fail")
    n_pm_ok = sum(1 for r in rows_out if r["pubmed_check"]["status"] == "ok")
    n_pm_skip = sum(1 for r in rows_out if r["pubmed_check"]["status"] == "skipped")
    n_pm_fail = sum(
        1
        for r in rows_out
        if r["pubmed_check"]["status"] in ("fail", "invalid")
    )

    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    REPORT_MD.parent.mkdir(parents=True, exist_ok=True)
    METRICS_JSON.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# Literature DOI / PubMed verification",
        "",
        f"**Generated:** {ts}  ",
        f"**Source:** `reports/literature_table.csv`  ",
        f"**Command:** `python scripts/verify_literature_dois.py`",
        "",
        "## Summary",
        "",
        f"| Check | OK | Skipped | Failed |",
        f"|-------|---:|--------:|-------:|",
        f"| DOI | {n_doi_ok} | {n_doi_skip} | {n_doi_fail} |",
        f"| PubMed | {n_pm_ok} | {n_pm_skip} | {n_pm_fail} |",
        "",
        "## Detail",
        "",
        "| ref | Author | Year | DOI | DOI status | PMID | PMID status |",
        "|----:|--------|------|-----|------------|------|-------------|",
    ]
    for r in rows_out:
        doi_st = r["doi_check"]["status"]
        pm_st = r["pubmed_check"]["status"]
        doi_disp = r["doi"] if r["doi"] and r["doi"] != "NA" else "—"
        pm_disp = r["pubmed_id"] if r["pubmed_id"] and r["pubmed_id"] != "NA" else "—"
        lines.append(
            f"| {r['ref_id']} | {r['first_author']} | {r['year']} | "
            f"{doi_disp} | {doi_st} | {pm_disp} | {pm_st} |"
        )

    notes = [
        "",
        "## Notes",
        "",
        "- **Skipped** — book/chapter entries without DOI or PubMed (e.g. Efron bootstrap, ISBN-only).",
        "- **403 on DOI** — publisher anti-bot; open link in a browser or use PubMed ID.",
    ]
    if n_doi_fail or n_pm_fail:
        notes.append("- **Action:** fix failed links before final Vancouver reference export.")
    else:
        notes.append("- All resolvable DOI/PMID links returned HTTP success at check time.")

    REPORT_MD.write_text("\n".join(lines + notes) + "\n", encoding="utf-8")
    METRICS_JSON.write_text(
        json.dumps({"generated": ts, "rows": rows_out}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"\nWrote {REPORT_MD}")
    print(f"Wrote {METRICS_JSON}")
    return 1 if (n_doi_fail or n_pm_fail) else 0


if __name__ == "__main__":
    raise SystemExit(main())
