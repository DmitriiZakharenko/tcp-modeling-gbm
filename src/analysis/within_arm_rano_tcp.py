"""
Within-arm RANO / DVH TCP and Cox analyses.

Tests whether DVH metrics predict RANO non-PD separately within each
fractionation scheme (where pooled EQD2 is constant or confounded).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy import stats
from sklearn.metrics import roc_auc_score

from src.config import FIGURES_DIR
from src.models.poisson_tcp import PoissonTCPModel

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

PRIMARY_SCHEMES: Tuple[Tuple[float, str], ...] = (
    (60.0, "60Gy_30fr"),
    (40.05, "40Gy_15fr"),
)
DVH_METRICS: Tuple[str, ...] = ("Dmean_gy", "D95_gy", "volume_cc", "gEUD_a10_gy", "HI_gy")
MIN_METRIC_STD: Dict[str, float] = {
    "Dmean_gy": 0.05,
    "D95_gy": 0.05,
    "volume_cc": 0.5,
    "gEUD_a10_gy": 0.05,
    "HI_gy": 0.05,
}


def _fit_poisson_auc(doses: np.ndarray, outcomes: np.ndarray) -> Tuple[float, float, float]:
    """Return in-sample AUC, LR p-value, D50; NaNs if fit fails."""
    if len(np.unique(outcomes)) < 2 or np.std(doses) < 1e-9:
        return np.nan, np.nan, np.nan
    p0 = float(outcomes.mean())
    nll_null = float(-np.sum(outcomes * np.log(p0) + (1.0 - outcomes) * np.log(1.0 - p0)))
    try:
        model = PoissonTCPModel(d50_init=float(np.median(doses)), gamma50_init=1.5)
        model.fit(doses, outcomes)
        preds = model.predict(doses)
        auc = float(roc_auc_score(outcomes, preds))
        nll = model.nll_
        lr_stat = 2 * ((-nll) - (-nll_null))
        lr_p = float(1.0 - stats.chi2.cdf(lr_stat, 2))
        return auc, lr_p, float(model.params_["D50_gy"])
    except (ValueError, RuntimeError):
        return np.nan, np.nan, np.nan


def within_arm_rano_tcp_table(frame: pd.DataFrame) -> pd.DataFrame:
    """
    Per scheme × DVH metric: Spearman and Poisson TCP AUC vs RANO non-PD.
    """
    rows: List[Dict[str, object]] = []
    for dose, label in PRIMARY_SCHEMES:
        sub = frame[(frame["rt_dose_gy"] == dose) & frame["rano_controlled_t1"].notna()].copy()
        if sub.empty:
            continue
        y = sub["rano_controlled_t1"].astype(float).to_numpy()
        pd_rate = float(1.0 - y.mean())
        for metric in DVH_METRICS:
            if metric not in sub.columns:
                continue
            x = sub[metric].astype(float)
            valid = x.notna()
            if valid.sum() < 15:
                continue
            xv = x[valid].to_numpy()
            yv = y[valid.to_numpy()]
            std = float(np.std(xv))
            min_std = MIN_METRIC_STD.get(metric, 0.01)
            rho, p_spearman = stats.spearmanr(xv, yv)
            auc, lr_p, d50 = (
                _fit_poisson_auc(xv, yv) if std >= min_std else (np.nan, np.nan, np.nan)
            )
            rows.append(
                {
                    "scheme": label,
                    "rt_dose_gy": dose,
                    "n_rano": int(valid.sum()),
                    "pd_rate": pd_rate,
                    "metric": metric,
                    "metric_std": std,
                    "metric_feasible": std >= min_std,
                    "endpoint": "rano_non_pd_t1",
                    "spearman_rho": float(rho),
                    "spearman_p": float(p_spearman),
                    "poisson_auc": auc,
                    "poisson_lr_p": lr_p,
                    "poisson_d50": d50,
                }
            )
    return pd.DataFrame(rows)


def cox_os_with_rano(frame: pd.DataFrame) -> pd.DataFrame:
    """
    Cox on RANO-evaluable patients: OS ~ age + sex + WHO PS + EQD2 + RANO non-PD.
    """
    if "rano_controlled_t1" not in frame.columns:
        return pd.DataFrame()
    sub = frame[frame["rano_controlled_t1"].notna()].copy()
    sub["event"] = 1
    sub["sex_M"] = (sub["sex"] == "M").astype(int)
    cols = ["survival_weeks", "event", "age", "sex_M", "who_status", "eqd2_gy", "rano_controlled_t1"]
    data = sub[cols].dropna()
    if len(data) < 20:
        return pd.DataFrame()

    cph = CoxPHFitter()
    cph.fit(data, duration_col="survival_weeks", event_col="event")
    summary = cph.summary.reset_index()
    if "covariate" in summary.columns:
        summary = summary.rename(columns={"covariate": "term"})
    elif "index" in summary.columns:
        summary = summary.rename(columns={"index": "term"})
    out = summary[["term", "coef", "exp(coef)", "p"]].copy()
    out.columns = ["term", "coef", "hazard_ratio", "p_value"]
    out["concordance_index"] = float(cph.concordance_index_)
    out["n_patients"] = len(data)
    return out


def plot_within_arm_rano_auc(table: pd.DataFrame, save_path: Optional[Path] = None) -> plt.Figure:
    """Grouped bar chart of Poisson AUC by scheme and DVH metric (RANO endpoint)."""
    plot_df = table[table["metric_feasible"]].dropna(subset=["poisson_auc"]).copy()
    if plot_df.empty:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.text(0.5, 0.5, "No feasible within-arm RANO TCP fits", ha="center", va="center")
        return fig

    metrics_order = [m for m in DVH_METRICS if m in plot_df["metric"].unique()]
    scheme_labels = [lab for _, lab in PRIMARY_SCHEMES]

    with plt.rc_context(DEFAULT_RC_PARAMS):
        fig, ax = plt.subplots(figsize=(8, 4.5))
        x = np.arange(len(metrics_order))
        width = 0.35
        colors = ["#2c7bb6", "#d7191c"]
        for i, label in enumerate(scheme_labels):
            sub = plot_df[plot_df["scheme"] == label].set_index("metric")
            vals = [sub.loc[m, "poisson_auc"] if m in sub.index else np.nan for m in metrics_order]
            ns = [int(sub.loc[m, "n_rano"]) if m in sub.index else 0 for m in metrics_order]
            bars = ax.bar(x + (i - 0.5) * width, vals, width, label=label, color=colors[i])
            for bar, n in zip(bars, ns):
                if not np.isnan(bar.get_height()):
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.02,
                        f"n={n}",
                        ha="center",
                        va="bottom",
                        fontsize=7,
                    )
        ax.axhline(0.5, color="gray", ls="--", lw=0.8, alpha=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels([m.replace("_gy", "").replace("_cc", "") for m in metrics_order], rotation=25, ha="right")
        ax.set_ylabel("Poisson TCP AUC (RANO non-PD)")
        ax.set_ylim(0.4, 1.05)
        ax.set_title("Within-arm DVH → RANO (Poisson TCP)")
        ax.legend(frameon=False, fontsize=9)
        plt.tight_layout()
        if save_path is not None:
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path)
    return fig


def run_within_arm_rano_analysis(
    frame: pd.DataFrame,
    metrics_dir: Path,
    figure_path: Optional[Path] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run within-arm RANO TCP table, Cox model, save CSVs and figure.

    Returns
    -------
    tuple
        (within_arm_table, cox_summary)
    """
    metrics_dir.mkdir(parents=True, exist_ok=True)
    table = within_arm_rano_tcp_table(frame)
    cox = cox_os_with_rano(frame)
    table.to_csv(metrics_dir / "within_arm_rano_tcp.csv", index=False)
    if not cox.empty:
        cox.to_csv(metrics_dir / "cox_os_with_rano.csv", index=False)
    plot_within_arm_rano_auc(table, save_path=figure_path or FIGURES_DIR / "06_within_arm_rano_tcp.png")
    return table, cox


def main() -> None:
    from src.config import DATA_PROCESSED, REPORTS_DIR

    frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    table, cox = run_within_arm_rano_analysis(frame, REPORTS_DIR / "metrics")
    print("Within-arm RANO TCP:")
    print(table.to_string(index=False))
    if not cox.empty:
        print("\nCox OS ~ covariates + RANO:")
        print(cox.to_string(index=False))


if __name__ == "__main__":
    main()
