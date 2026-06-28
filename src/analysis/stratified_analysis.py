"""
Stratified clinical and DVH survival analyses for CFB-GBM.

Focuses on verified prognostic signals (scheme, WHO PS) and hypothesis-driven
within-arm DVH tests where dose heterogeneity allows.
"""

from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats

from src.config import DATA_PROCESSED, FIGURES_DIR, REPORTS_DIR, TREATMENT_TSV

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
COLORS = ["#2c7bb6", "#d7191c", "#fdae61", "#1a9641", "#762a83"]


def _prep(frame: pd.DataFrame) -> pd.DataFrame:
    data = frame.copy()
    data["event"] = 1
    data["sex_M"] = (data["sex"] == "M").astype(int)
    data["scheme_60gy"] = (data["rt_dose_gy"] == 60.0).astype(int)
    return data


def clinical_cox_summary(frame: pd.DataFrame) -> pd.DataFrame:
    """
    Cox model: OS ~ age + sex + WHO PS + scheme (60 Gy indicator).

    Parameters
    ----------
    frame : pd.DataFrame
        Modeling table.

    Returns
    -------
    pd.DataFrame
        Hazard ratios and p-values.
    """
    data = _prep(frame)
    cph = CoxPHFitter()
    cph.fit(
        data[["survival_weeks", "event", "age", "sex_M", "who_status", "scheme_60gy"]].dropna(),
        duration_col="survival_weeks",
        event_col="event",
    )
    summary = cph.summary.reset_index()
    if "covariate" in summary.columns:
        summary = summary.rename(columns={"covariate": "term"})
    elif "index" in summary.columns:
        summary = summary.rename(columns={"index": "term"})
    out = summary[["term", "coef", "exp(coef)", "p"]].copy()
    out.columns = ["term", "coef", "hazard_ratio", "p_value"]
    out["concordance_index"] = float(cph.concordance_index_)
    return out


def who_ps_os_table(frame: pd.DataFrame) -> pd.DataFrame:
    """Median OS by WHO PS with Kruskal-Wallis p-value."""
    rows = []
    for ps in sorted(frame["who_status"].dropna().unique()):
        sub = frame[frame["who_status"] == ps]
        rows.append(
            {
                "who_status": int(ps),
                "n": len(sub),
                "os_median_wk": float(sub["survival_weeks"].median()),
                "os_q25_wk": float(sub["survival_weeks"].quantile(0.25)),
                "os_q75_wk": float(sub["survival_weeks"].quantile(0.75)),
            }
        )
    table = pd.DataFrame(rows)
    groups = [frame.loc[frame["who_status"] == ps, "survival_weeks"].values for ps in sorted(frame["who_status"].unique())]
    table.attrs["kruskal_p"] = float(stats.kruskal(*groups).pvalue)
    return table


def within_arm_dvh_tests(frame: pd.DataFrame) -> pd.DataFrame:
    """
    Spearman correlation OS vs DVH metrics within 60 Gy and 40 Gy arms.

    Parameters
    ----------
    frame : pd.DataFrame
        Modeling table.

    Returns
    -------
    pd.DataFrame
        One row per scheme × metric.
    """
    metrics = ["Dmean_gy", "D95_gy", "volume_cc", "HI_gy", "gEUD_a10_gy"]
    rows = []
    for dose, label in [(60.0, "60Gy_30fr"), (40.05, "40Gy_15fr")]:
        sub = frame[frame["rt_dose_gy"] == dose]
        for metric in metrics:
            r, p = stats.spearmanr(sub["survival_weeks"], sub[metric])
            rows.append(
                {
                    "scheme": label,
                    "n": len(sub),
                    "metric": metric,
                    "spearman_rho": float(r),
                    "p_value": float(p),
                }
            )
    return pd.DataFrame(rows)


def hypofractionated_volume_cox(frame: pd.DataFrame) -> pd.DataFrame:
    """
    Cox within 40.05 Gy arm: OS ~ volume + age + sex + WHO PS.

    Parameters
    ----------
    frame : pd.DataFrame
        Modeling table.

    Returns
    -------
    pd.DataFrame
        Cox summary for hypofractionated subgroup (n≈61).
    """
    sub = _prep(frame[frame["rt_dose_gy"] == 40.05])
    cph = CoxPHFitter()
    cph.fit(
        sub[["survival_weeks", "event", "volume_cc", "age", "sex_M", "who_status"]].dropna(),
        duration_col="survival_weeks",
        event_col="event",
    )
    summary = cph.summary.reset_index()
    if "covariate" in summary.columns:
        summary = summary.rename(columns={"covariate": "term"})
    elif "index" in summary.columns:
        summary = summary.rename(columns={"index": "term"})
    out = summary[["term", "coef", "exp(coef)", "p"]].copy()
    out.columns = ["term", "coef", "hazard_ratio", "p_value"]
    out["concordance_index"] = float(cph.concordance_index_)
    out["n_patients"] = len(sub)
    return out


def plot_clinical_prognosis(
    frame: pd.DataFrame,
    save_path: Optional[Path] = None,
) -> plt.Figure:
    """
    Two-panel figure: KM by fractionation scheme and by WHO PS (0–2).

    Parameters
    ----------
    frame : pd.DataFrame
        Modeling table.
    save_path : pathlib.Path, optional
        PNG output path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    data = _prep(frame)
    kmf = KaplanMeierFitter()

    with plt.rc_context(DEFAULT_RC_PARAMS):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

        for dose, color, label in [
            (60.0, COLORS[0], "60 Gy / 30 fr"),
            (40.05, COLORS[1], "40.05 Gy / 15 fr"),
        ]:
            mask = data["rt_dose_gy"] == dose
            sub = data.loc[mask]
            kmf.fit(sub["survival_weeks"], event_observed=sub["event"],
                    label=f"{label} (n={mask.sum()}, med={sub['survival_weeks'].median():.0f} wk)")
            kmf.plot_survival_function(ax=axes[0], ci_show=True, color=color)
        g60 = data.loc[data["rt_dose_gy"] == 60, "survival_weeks"]
        g40 = data.loc[data["rt_dose_gy"] == 40.05, "survival_weeks"]
        lr = logrank_test(g60, g40, event_observed_A=np.ones(len(g60)), event_observed_B=np.ones(len(g40)))
        axes[0].set_title(f"OS by fractionation (log-rank p={lr.p_value:.1e})")
        axes[0].set_xlabel("Survival (weeks)")
        axes[0].set_ylabel("Survival probability")
        axes[0].legend(fontsize=8, frameon=False)

        ps_data = data[data["who_status"].isin([0, 1, 2])]
        for i, ps in enumerate([0, 1, 2]):
            mask = ps_data["who_status"] == ps
            sub = ps_data.loc[mask]
            kmf.fit(sub["survival_weeks"], event_observed=sub["event"],
                    label=f"PS {int(ps)} (n={mask.sum()}, med={sub['survival_weeks'].median():.0f} wk)")
            kmf.plot_survival_function(ax=axes[1], ci_show=True, color=COLORS[i])
        ps_table = who_ps_os_table(frame)
        axes[1].set_title(f"OS by WHO PS (KW p={ps_table.attrs['kruskal_p']:.1e})")
        axes[1].set_xlabel("Survival (weeks)")
        axes[1].set_ylabel("Survival probability")
        axes[1].legend(fontsize=8, frameon=False)

        plt.tight_layout()
        if save_path is not None:
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path)

    return fig


def run_stratified_analysis(
    frame: pd.DataFrame,
    metrics_dir: Path,
    figure_path: Optional[Path] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Run all stratified tests and save CSVs + figure.

    Returns
    -------
    tuple
        clinical_cox, who_ps, within_arm, hypo_volume_cox
    """
    metrics_dir.mkdir(parents=True, exist_ok=True)
    clinical = clinical_cox_summary(frame)
    who = who_ps_os_table(frame)
    within = within_arm_dvh_tests(frame)
    hypo = hypofractionated_volume_cox(frame)

    clinical.to_csv(metrics_dir / "clinical_cox_summary.csv", index=False)
    who.to_csv(metrics_dir / "who_ps_os.csv", index=False)
    within.to_csv(metrics_dir / "within_arm_dvh_spearman.csv", index=False)
    hypo.to_csv(metrics_dir / "hypofractionated_volume_cox.csv", index=False)

    plot_clinical_prognosis(frame, save_path=figure_path or FIGURES_DIR / "04_clinical_prognosis.png")
    return clinical, who, within, hypo


def main() -> None:
    """Run stratified analysis on modeling cohort."""
    frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    clinical, who, within, hypo = run_stratified_analysis(frame, REPORTS_DIR / "metrics")
    print("Clinical Cox (OS ~ age + sex + PS + scheme):")
    print(clinical.to_string(index=False))
    print("\nWHO PS OS:")
    print(who.to_string(index=False))
    print("\nWithin-arm DVH Spearman (40 Gy volume p):",
          within.loc[(within.scheme == "40Gy_15fr") & (within.metric == "volume_cc"), "p_value"].iloc[0])
    print("\n40 Gy volume Cox:")
    print(hypo[hypo["term"] == "volume_cc"].to_string(index=False))


if __name__ == "__main__":
    main()
