"""
Compare TCP model quality: OS median-split vs RANO non-PD endpoint.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from src.config import FIGURES_DIR, RANDOM_SEED
from src.models.poisson_tcp import PoissonTCPModel


def _evaluate_models(doses, os_outcomes, rano_outcomes, os_def, rano_def):
    from src.reporting.update_results import evaluate_logistic, evaluate_poisson, evaluate_probit

    rows: List[Dict[str, object]] = []
    for model_fn, name in (
        (evaluate_poisson, "poisson_tcp"),
        (evaluate_logistic, "logistic_tcp"),
        (evaluate_probit, "probit_tcp"),
    ):
        os_m = model_fn(doses, os_outcomes, "eqd2_gy", os_def)
        rano_m = model_fn(doses, rano_outcomes, "eqd2_gy", rano_def)
        for endpoint, m, outcomes in (
            ("os_median_proxy", os_m, os_outcomes),
            ("rano_non_pd_t1", rano_m, rano_outcomes),
        ):
            rows.append(
                {
                    "model": name,
                    "endpoint": endpoint,
                    "outcome_definition": m["outcome_definition"],
                    "n": len(outcomes),
                    "event_rate": float(outcomes.mean()),
                    "roc_auc_insample": m["roc_auc_insample"],
                    "roc_auc_cv_mean": m["roc_auc_cv_mean"],
                    "roc_auc_cv_std": m["roc_auc_cv_std"],
                    "aic": m["aic"],
                    "lr_p_value": m["lr_p_value"],
                    "D50_gy": m["D50_gy"],
                }
            )
    return rows

DEFAULT_RC_PARAMS = {
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
}


def rano_evaluable(frame: pd.DataFrame) -> pd.DataFrame:
    """Patients with RANO t0→t1 label and EQD2."""
    if "rano_controlled_t1" not in frame.columns:
        return frame.iloc[0:0].copy()
    sub = frame[frame["rano_controlled_t1"].notna() & frame["eqd2_gy"].notna()].copy()
    return sub


def compare_outcomes_on_rano_subset(frame: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Fit Poisson/Logistic/Probit on RANO-evaluable patients for both endpoints.

    Returns
    -------
    tuple
        (comparison_rows, rano_category_counts)
    """
    sub = rano_evaluable(frame)
    if sub.empty:
        return pd.DataFrame(), pd.DataFrame()

    median_os = float(sub["survival_weeks"].median())
    os_outcomes = (sub["survival_weeks"] >= median_os).astype(float).to_numpy()
    rano_outcomes = sub["rano_controlled_t1"].astype(float).to_numpy()
    doses = sub["eqd2_gy"].to_numpy()
    n = len(sub)
    pd_rate = float(1.0 - rano_outcomes.mean())

    os_def = f"OS >= median ({median_os:.0f} wk), n={n}"
    rano_def = f"RANO non-PD at t1, n={n}, PD rate={pd_rate:.1%}"

    rows: List[Dict[str, object]] = []
    rows = _evaluate_models(doses, os_outcomes, rano_outcomes, os_def, rano_def)

    # Dose–outcome association on same subset
    r_eqd2_os, p_eqd2_os = stats.pointbiserialr(os_outcomes, doses)
    r_eqd2_rano, p_eqd2_rano = stats.pointbiserialr(rano_outcomes, doses)

    assoc = pd.DataFrame(
        [
            {
                "pair": "EQD2_vs_os_median_proxy",
                "n": n,
                "point_biserial_r": float(r_eqd2_os),
                "p_value": float(p_eqd2_os),
            },
            {
                "pair": "EQD2_vs_rano_non_pd",
                "n": n,
                "point_biserial_r": float(r_eqd2_rano),
                "p_value": float(p_eqd2_rano),
            },
        ]
    )

    if "rano_t0_t1" in sub.columns:
        rano_counts = sub["rano_t0_t1"].value_counts(dropna=True).rename_axis("category").reset_index(name="n")
    else:
        rano_counts = pd.DataFrame()

    comparison = pd.DataFrame(rows)
    comparison.attrs["association"] = assoc
    comparison.attrs["pd_rate"] = pd_rate
    comparison.attrs["n_rano"] = n
    return comparison, rano_counts


def plot_auc_comparison(comparison: pd.DataFrame, save_path: Path | None = None) -> plt.Figure:
    """Bar chart: in-sample AUC by model and endpoint."""
    if comparison.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No RANO data", ha="center", va="center")
        return fig

    pivot = comparison.pivot(index="model", columns="endpoint", values="roc_auc_insample")
    pivot = pivot.reindex(columns=["os_median_proxy", "rano_non_pd_t1"])
    pivot.columns = ["OS median proxy", "RANO non-PD (t1)"]

    with plt.rc_context(DEFAULT_RC_PARAMS):
        fig, ax = plt.subplots(figsize=(7, 4.5))
        x = np.arange(len(pivot))
        width = 0.35
        ax.bar(x - width / 2, pivot.iloc[:, 0], width, label=pivot.columns[0], color="#2c7bb6")
        ax.bar(x + width / 2, pivot.iloc[:, 1], width, label=pivot.columns[1], color="#d7191c")
        ax.axhline(0.5, color="gray", ls="--", lw=0.8, alpha=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels([m.replace("_tcp", "") for m in pivot.index], rotation=0)
        ax.set_ylabel("ROC AUC (in-sample)")
        ax.set_ylim(0.4, max(0.85, pivot.values.max() + 0.05))
        ax.set_title(f"TCP discrimination on RANO-evaluable subset (n={int(comparison['n'].iloc[0])})")
        ax.legend(frameon=False, fontsize=9)
        plt.tight_layout()
        if save_path is not None:
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path)

    return fig


def run_rano_tcp_comparison(
    frame: pd.DataFrame,
    metrics_dir: Path,
    figure_path: Path | None = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Run comparison and save CSVs + figure.

    Returns
    -------
    tuple
        comparison, rano_counts, association
    """
    metrics_dir.mkdir(parents=True, exist_ok=True)
    comparison, rano_counts = compare_outcomes_on_rano_subset(frame)
    comparison.to_csv(metrics_dir / "tcp_outcome_comparison_rano_subset.csv", index=False)
    rano_counts.to_csv(metrics_dir / "rano_category_counts_modeling.csv", index=False)

    assoc = comparison.attrs.get("association", pd.DataFrame())
    if isinstance(assoc, pd.DataFrame) and not assoc.empty:
        assoc.to_csv(metrics_dir / "rano_subset_dose_outcome_association.csv", index=False)

    plot_auc_comparison(
        comparison,
        save_path=figure_path or FIGURES_DIR / "05_rano_vs_os_tcp_auc.png",
    )
    return comparison, rano_counts, assoc if isinstance(assoc, pd.DataFrame) else pd.DataFrame()
