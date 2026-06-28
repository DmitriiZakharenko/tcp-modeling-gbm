"""
Build literature TCP parameter comparison table (assignment Part VI).

Combines our bootstrap Poisson D50/gamma50 with published GBM / TCP references.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.config import REPORTS_DIR


def build_literature_tcp_d50_table(metrics_dir: Path | None = None) -> pd.DataFrame:
    """Assemble comparison table: our estimates vs published TCP/DVH literature."""
    metrics_dir = metrics_dir or REPORTS_DIR / "metrics"
    boot_path = metrics_dir / "poisson_tcp_bootstrap_ci.csv"

    our_d50 = our_gamma = d50_ci = gamma_ci = ""
    if boot_path.exists():
        boot = pd.read_csv(boot_path).set_index("parameter")
        our_d50 = f"{boot.loc['D50_gy', 'estimate']:.1f}"
        d50_ci = f"[{boot.loc['D50_gy', 'ci_lower']:.1f}, {boot.loc['D50_gy', 'ci_upper']:.1f}]"
        our_gamma = f"{boot.loc['gamma50', 'estimate']:.2f}"
        gamma_ci = f"[{boot.loc['gamma50', 'ci_lower']:.2f}, {boot.loc['gamma50', 'ci_upper']:.2f}]"

    rows = [
        {
            "source": "This study (CFB-GBM)",
            "ref_id": "—",
            "model": "Poisson TCP (EQD2)",
            "endpoint": "OS ≥ median (exploratory proxy)",
            "n": 190,
            "D50_gy": our_d50,
            "D50_ci_95": d50_ci,
            "gamma50": our_gamma,
            "gamma50_ci_95": gamma_ci,
            "comparable_to_ours": "reference",
            "notes": "Confounded by fractionation scheme; not formal LC",
        },
        {
            "source": "Maitre et al. 2020",
            "ref_id": "4",
            "model": "Poisson / LQ TCP (review)",
            "endpoint": "Various (LC, NTCP)",
            "n": "review",
            "D50_gy": "40–80",
            "D50_ci_95": "—",
            "gamma50": "1–5",
            "gamma50_ci_95": "—",
            "comparable_to_ours": "partial",
            "notes": "Site-specific; GBM LC data sparse in routine cohorts",
        },
        {
            "source": "Ohri et al. 2017",
            "ref_id": "7",
            "model": "TCP/NTCP (review)",
            "endpoint": "LC preferred over OS",
            "n": "review",
            "D50_gy": "—",
            "D50_ci_95": "—",
            "gamma50": "—",
            "gamma50_ci_95": "—",
            "comparable_to_ours": "no",
            "notes": "Highlights need for dose heterogeneity ≥1 Gy and LC endpoints",
        },
        {
            "source": "Embring et al. 2020",
            "ref_id": "16",
            "model": "DVH metrics (not TCP fit)",
            "endpoint": "OS",
            "n": 120,
            "D50_gy": "—",
            "D50_ci_95": "—",
            "gamma50": "—",
            "gamma50_ci_95": "—",
            "comparable_to_ours": "partial",
            "notes": "GTV Dmean/volume prognostic; no Poisson D50 reported",
        },
        {
            "source": "Gardner et al. 2024",
            "ref_id": "5",
            "model": "Radiobiology review",
            "endpoint": "Multi-scale modelling",
            "n": "review",
            "D50_gy": "context-dependent",
            "D50_ci_95": "—",
            "gamma50": "context-dependent",
            "gamma50_ci_95": "—",
            "comparable_to_ours": "partial",
            "notes": "Modern TCP modelling requires appropriate endpoints and dose spread",
        },
        {
            "source": "Okunieff et al. 1995",
            "ref_id": "assignment",
            "model": "Poisson TCP (meta-analysis)",
            "endpoint": "Tumour control (mixed sites)",
            "n": "review",
            "D50_gy": "~50–70",
            "D50_ci_95": "—",
            "gamma50": "~1–4",
            "gamma50_ci_95": "—",
            "comparable_to_ours": "partial",
            "notes": "Course reference; not GBM-specific; LC-based trials",
        },
    ]
    return pd.DataFrame(rows)


def run_literature_tcp_comparison(output_dir: Path | None = None) -> pd.DataFrame:
    """Write literature_tcp_d50_comparison.csv."""
    out_dir = output_dir or REPORTS_DIR / "metrics"
    out_dir.mkdir(parents=True, exist_ok=True)
    table = build_literature_tcp_d50_table(out_dir)
    table.to_csv(out_dir / "literature_tcp_d50_comparison.csv", index=False)
    return table


def main() -> None:
    table = run_literature_tcp_comparison()
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
