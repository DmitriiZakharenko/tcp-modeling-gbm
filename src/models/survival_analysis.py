"""
Survival analysis: Kaplan-Meier curves and Cox proportional hazards.

Uses overall survival in weeks from the modeling table.
"""

from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

from src.config import FIGURES_DIR, RANDOM_SEED

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

COX_COVARIATES: Tuple[str, ...] = ("eqd2_gy", "Dmean_gy", "age", "sex")
COLORS = ["#2c7bb6", "#d7191c", "#fdae61", "#1a9641"]


def prepare_survival_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare lifelines-compatible survival DataFrame.

    Parameters
    ----------
    frame : pd.DataFrame
        Modeling table with ``survival_weeks``, covariates, and ``sex``.

    Returns
    -------
    pd.DataFrame
        Columns: ``duration_weeks``, ``event_observed``, covariates, ``sex_M``.
    """
    data = frame.copy()
    data["duration_weeks"] = data["survival_weeks"].astype(float)
    data["event_observed"] = 1
    data["sex_M"] = (data["sex"] == "M").astype(int)
    return data


def kaplan_meier_by_dose(
    frame: pd.DataFrame,
    dose_threshold: float = 50.0,
    save_path: Optional[Path] = None,
) -> Tuple[pd.DataFrame, plt.Figure]:
    """
    Plot Kaplan-Meier curves stratified by EQD2 threshold.

    Parameters
    ----------
    frame : pd.DataFrame
        Modeling table.
    dose_threshold : float, optional
        Split EQD2 at this value (Gy).
    save_path : pathlib.Path, optional
        Output PNG path.

    Returns
    -------
    tuple
        ``(logrank_summary_df, figure)``
    """
    data = prepare_survival_frame(frame)
    high = data["eqd2_gy"] >= dose_threshold
    low = ~high

    kmf = KaplanMeierFitter()
    with plt.rc_context(DEFAULT_RC_PARAMS):
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for mask, label, color in [
            (high, f"EQD2 ≥ {dose_threshold} Gy (n={high.sum()})", COLORS[0]),
            (low, f"EQD2 < {dose_threshold} Gy (n={low.sum()})", COLORS[1]),
        ]:
            kmf.fit(
                data.loc[mask, "duration_weeks"],
                event_observed=data.loc[mask, "event_observed"],
                label=label,
            )
            kmf.plot_survival_function(ax=ax, ci_show=True, color=color)

        ax.set_xlabel("Overall survival (weeks)")
        ax.set_ylabel("Survival probability")
        ax.set_title("Kaplan-Meier by EQD2 group")
        ax.set_ylim(0, 1.05)
        ax.legend(fontsize=9, frameon=False)
        plt.tight_layout()
        if save_path is not None:
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path)

    lr = logrank_test(
        data.loc[high, "duration_weeks"],
        data.loc[low, "duration_weeks"],
        event_observed_A=data.loc[high, "event_observed"],
        event_observed_B=data.loc[low, "event_observed"],
    )
    summary = pd.DataFrame(
        [
            {
                "comparison": f"EQD2>={dose_threshold}_vs_<{dose_threshold}",
                "n_high": int(high.sum()),
                "n_low": int(low.sum()),
                "logrank_chi2": float(lr.test_statistic),
                "logrank_p_value": float(lr.p_value),
            }
        ]
    )
    return summary, fig


def fit_cox_model(
    frame: pd.DataFrame,
    covariates: Sequence[str] = COX_COVARIATES,
) -> Tuple[CoxPHFitter, pd.DataFrame]:
    """
    Fit Cox PH model and return hazard ratio summary.

    Parameters
    ----------
    frame : pd.DataFrame
        Modeling table.
    covariates : sequence of str
        Covariate column names; ``sex`` is encoded as ``sex_M``.

    Returns
    -------
    tuple
        ``(fitted CoxPHFitter, summary DataFrame with HR, CI, p)``

    Raises
    ------
    ValueError
        If a covariate column is missing.
    """
    data = prepare_survival_frame(frame)
    cox_cols: List[str] = []
    for col in covariates:
        if col == "sex":
            cox_cols.append("sex_M")
        elif col in data.columns:
            cox_cols.append(col)
        else:
            raise ValueError(f"Missing covariate column: {col}")

    cox_df = data[["duration_weeks", "event_observed", *cox_cols]].dropna()
    cph = CoxPHFitter()
    cph.fit(cox_df, duration_col="duration_weeks", event_col="event_observed")

    summary = cph.summary.copy()
    summary = summary.reset_index().rename(columns={"index": "covariate", "covariate": "term"})
    if "term" not in summary.columns:
        summary = summary.rename(columns={"covariate": "term"})
    summary["hazard_ratio"] = np.exp(summary["coef"])
    summary["hr_ci_lower"] = np.exp(summary["coef lower 95%"])
    summary["hr_ci_upper"] = np.exp(summary["coef upper 95%"])
    out = summary[["term", "hazard_ratio", "hr_ci_lower", "hr_ci_upper", "p"]].copy()
    out["concordance_index"] = float(cph.concordance_index_)
    return cph, out


def plot_cox_forest(
    cox_summary: pd.DataFrame,
    save_path: Optional[Path] = None,
) -> plt.Figure:
    """
    Plot Cox hazard ratio forest plot.

    Parameters
    ----------
    cox_summary : pd.DataFrame
        Output of ``fit_cox_model`` (HR summary rows only).
    hr_rows = cox_summary.dropna(subset=["hazard_ratio"])
    save_path : pathlib.Path, optional
        Output PNG path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    hr_rows = cox_summary[cox_summary["term"].isin(COX_COVARIATES + ("sex_M",))].copy()
    if hr_rows.empty:
        hr_rows = cox_summary.dropna(subset=["hazard_ratio"]).copy()

    labels = hr_rows["term"].tolist()
    hrs = hr_rows["hazard_ratio"].to_numpy()
    lo = hr_rows["hr_ci_lower"].to_numpy()
    hi = hr_rows["hr_ci_upper"].to_numpy()
    y_pos = np.arange(len(labels))

    with plt.rc_context(DEFAULT_RC_PARAMS):
        fig, ax = plt.subplots(figsize=(7, 0.8 * len(labels) + 1.5))
        ax.errorbar(
            hrs,
            y_pos,
            xerr=[hrs - lo, hi - hrs],
            fmt="o",
            color=COLORS[0],
            ecolor=COLORS[0],
            capsize=4,
        )
        ax.axvline(1.0, color="0.4", ls="--", lw=1)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.set_xlabel("Hazard ratio (95% CI)")
        c_index = float(hr_rows["concordance_index"].iloc[0]) if "concordance_index" in hr_rows.columns else np.nan
        ax.set_title(f"Cox PH forest plot (C-index = {c_index:.3f})")
        plt.tight_layout()
        if save_path is not None:
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path)

    return fig


def run_survival_analysis(
    frame: pd.DataFrame,
    km_path: Optional[Path] = None,
    forest_path: Optional[Path] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, plt.Figure, plt.Figure]:
    """
    Run KM and Cox analysis; save figures.

    Parameters
    ----------
    frame : pd.DataFrame
        Modeling table.
    km_path, forest_path : pathlib.Path, optional
        Figure output paths.

    Returns
    -------
    tuple
        ``(km_logrank_df, cox_summary_df, km_fig, forest_fig)``
    """
    km_summary, km_fig = kaplan_meier_by_dose(
        frame,
        save_path=km_path or FIGURES_DIR / "03_kaplan_meier_eqd2.png",
    )
    _, cox_summary = fit_cox_model(frame)
    forest_fig = plot_cox_forest(
        cox_summary,
        save_path=forest_path or FIGURES_DIR / "03_cox_forest.png",
    )
    return km_summary, cox_summary, km_fig, forest_fig


def main() -> None:
    """Run survival analysis on modeling cohort."""
    from src.config import DATA_PROCESSED

    frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    km_summary, cox_summary, _, _ = run_survival_analysis(frame)

    print(f"Cohort: n={len(frame)}, all events observed (OS)")
    print("\nLog-rank (EQD2 groups):")
    print(km_summary.to_string(index=False))
    print("\nCox PH summary:")
    print(cox_summary[["term", "hazard_ratio", "hr_ci_lower", "hr_ci_upper", "p"]].to_string(index=False))
    print(f"\nRandom seed (project default): {RANDOM_SEED}")


if __name__ == "__main__":
    main()
