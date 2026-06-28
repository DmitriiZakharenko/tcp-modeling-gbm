"""
Regenerate project results report from verified data files.

Run before each commit that changes analysis outputs::

    python -m src.reporting.update_results

Writes CSV metrics to ``reports/metrics/`` and a human-readable
``reports/RESULTS.md`` with figure links.
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import accuracy_score, brier_score_loss, confusion_matrix, roc_auc_score
from sklearn.model_selection import StratifiedKFold

from src.config import DATA_PROCESSED, FIGURES_DIR, RANDOM_SEED, REPORTS_DIR
from src.data.dvh_calculator import SCALAR_METRIC_KEYS
from src.analysis.rano_tcp_comparison import run_rano_tcp_comparison
from src.analysis.rano_multivariable import run_rano_multivariable_40gy
from src.analysis.rano_prediction_suite import run_rano_prediction_suite
from src.analysis.validate_rano_volumes import run_volume_validation
from src.analysis.within_arm_rano_tcp import run_within_arm_rano_analysis
from src.analysis.confounding_audit import run_confounding_audit, tcp_feasibility_summary
from src.analysis.stratified_analysis import run_stratified_analysis
from src.models.bootstrap_ci import bootstrap_tcp_params
from src.models.eud_tcp import EUDTCPModel
from src.models.model_comparison import run_model_comparison
from src.models.survival_analysis import run_survival_analysis
from src.models.logistic_tcp import LogisticTCPModel
from src.models.poisson_tcp import PoissonTCPModel
from src.models.probit_tcp import ProbitTCPModel

METRICS_DIR = REPORTS_DIR / "metrics"
RESULTS_MD = REPORTS_DIR / "RESULTS.md"
GIT_HEAD = Path(__file__).resolve().parents[2]


def _git_short_hash() -> str:
    import subprocess

    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=GIT_HEAD,
            text=True,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def _save_csv(frame: pd.DataFrame, name: str) -> Path:
    METRICS_DIR.mkdir(parents=True, exist_ok=True)
    path = METRICS_DIR / name
    frame.to_csv(path, index=False)
    return path


def cohort_metrics(frame: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    """Build cohort summary table."""
    included = cohort[cohort["included"]]
    rows = [
        {"metric": "total_patients_tsv", "value": len(cohort), "unit": "n"},
        {"metric": "included_cohort_rules", "value": int(included["included"].sum()), "unit": "n"},
        {"metric": "modeling_table_rows", "value": len(frame), "unit": "n"},
        {"metric": "excluded", "value": int((~cohort["included"]).sum()), "unit": "n"},
        {"metric": "fractionation_60gy_30fr", "value": int((frame["rt_dose_gy"] == 60.0).sum()), "unit": "n"},
        {"metric": "fractionation_40gy_15fr", "value": int((frame["rt_dose_gy"] == 40.05).sum()), "unit": "n"},
        {"metric": "age_median_yr", "value": float(frame["age"].median()), "unit": "yr"},
        {"metric": "os_median_wk", "value": float(frame["survival_weeks"].median()), "unit": "wk"},
        {"metric": "sex_male", "value": int((frame["sex"] == "M").sum()), "unit": "n"},
        {"metric": "sex_female", "value": int((frame["sex"] == "F").sum()), "unit": "n"},
        {"metric": "dmean_median_gy", "value": float(frame["Dmean_gy"].median()), "unit": "Gy"},
        {"metric": "gtv_volume_median_cc", "value": float(frame["volume_cc"].median()), "unit": "cc"},
    ]
    if "rano_controlled_t1" in frame.columns:
        n_rano = int(frame["rano_controlled_t1"].notna().sum())
        pd_rate = 1.0 - float(frame["rano_controlled_t1"].dropna().mean())
        rows.extend(
            [
                {"metric": "with_rano_t0_t1", "value": n_rano, "unit": "n"},
                {"metric": "rano_pd_rate_t1", "value": pd_rate, "unit": "proportion"},
            ]
        )
    return pd.DataFrame(rows)


def survival_by_fractionation(frame: pd.DataFrame) -> pd.DataFrame:
    """Descriptive OS by primary fractionation schemes."""
    rows: List[Dict[str, Any]] = []
    for dose, label in [(60.0, "60Gy_30fr"), (40.05, "40Gy_15fr")]:
        sub = frame.loc[frame["rt_dose_gy"] == dose, "survival_weeks"]
        rows.append(
            {
                "scheme": label,
                "rt_dose_gy": dose,
                "n": len(sub),
                "os_median_wk": float(sub.median()),
                "os_q25_wk": float(sub.quantile(0.25)),
                "os_q75_wk": float(sub.quantile(0.75)),
                "p_value": np.nan,
            }
        )
    s60 = frame.loc[frame["rt_dose_gy"] == 60.0, "survival_weeks"]
    s40 = frame.loc[frame["rt_dose_gy"] == 40.05, "survival_weeks"]
    _, p_value = stats.mannwhitneyu(s60, s40, alternative="two-sided")
    rows.append(
        {
            "scheme": "mann_whitney_p",
            "rt_dose_gy": np.nan,
            "n": np.nan,
            "os_median_wk": np.nan,
            "os_q25_wk": np.nan,
            "os_q75_wk": np.nan,
            "p_value": float(p_value),
        }
    )
    return pd.DataFrame(rows)


def association_metrics(frame: pd.DataFrame, outcomes: np.ndarray) -> pd.DataFrame:
    """Correlations between dose and OS / outcome proxy."""
    rows = []
    for col, label in [("eqd2_gy", "EQD2"), ("Dmean_gy", "Dmean"), ("survival_weeks", "OS_weeks")]:
        if col == "survival_weeks":
            r, p = stats.pearsonr(frame["eqd2_gy"], frame["survival_weeks"])
            rows.append({"pair": "EQD2_vs_OS", "statistic": "pearson_r", "value": float(r), "p_value": float(p)})
            r2, p2 = stats.pearsonr(frame["Dmean_gy"], frame["survival_weeks"])
            rows.append({"pair": "Dmean_vs_OS", "statistic": "pearson_r", "value": float(r2), "p_value": float(p2)})
        else:
            r, p = stats.pointbiserialr(outcomes, frame[col].to_numpy())
            rows.append({"pair": f"{label}_vs_outcome_proxy", "statistic": "point_biserial_r",
                         "value": float(r), "p_value": float(p)})
    return pd.DataFrame(rows)


def evaluate_tcp_model(
    model_name: str,
    model_factory,
    doses: np.ndarray,
    outcomes: np.ndarray,
    dose_label: str,
    outcome_definition: str,
    k_params: int,
    param2_name: str,
    n_folds: int = 5,
) -> Dict[str, Any]:
    """Fit any binary TCP model and return quality metrics."""
    n = len(outcomes)
    p0 = float(outcomes.mean())
    nll_null = float(-np.sum(outcomes * np.log(p0) + (1.0 - outcomes) * np.log(1.0 - p0)))

    model = model_factory()
    model.fit(doses, outcomes)
    preds = model.predict(doses)
    nll = model.nll_
    aic = 2 * k_params + 2 * nll
    bic = k_params * np.log(n) + 2 * nll
    lr_stat = 2 * ((-nll) - (-nll_null))
    lr_p = float(1.0 - stats.chi2.cdf(lr_stat, k_params))
    mcfadden = 1.0 - nll / nll_null

    yhat = (preds >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(outcomes, yhat).ravel()

    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=RANDOM_SEED)
    cv_aucs = []
    for train_idx, test_idx in skf.split(doses, outcomes):
        cv_model = model_factory()
        cv_model.fit(doses[train_idx], outcomes[train_idx])
        cv_preds = cv_model.predict(doses[test_idx])
        cv_aucs.append(float(roc_auc_score(outcomes[test_idx], cv_preds)))

    return {
        "model": model_name,
        "dose_variable": dose_label,
        "n": n,
        "outcome_definition": outcome_definition,
        "D50_gy": float(model.params_["D50_gy"]),
        param2_name: float(model.params_[param2_name]),
        "param2_name": param2_name,
        "nll": nll,
        "nll_null": nll_null,
        "aic": aic,
        "bic": bic,
        "aic_null": 2 * nll_null,
        "bic_null": 2 * nll_null,
        "lr_chi2": lr_stat,
        "lr_p_value": lr_p,
        "mcfadden_pseudo_r2": mcfadden,
        "roc_auc_insample": float(roc_auc_score(outcomes, preds)),
        "roc_auc_cv_mean": float(np.mean(cv_aucs)),
        "roc_auc_cv_std": float(np.std(cv_aucs)),
        "brier_score": float(brier_score_loss(outcomes, preds)),
        "brier_null": float(brier_score_loss(outcomes, np.full(n, p0))),
        "accuracy_at_0p5": float(accuracy_score(outcomes, yhat)),
        "sensitivity": float(tp / (tp + fn)),
        "specificity": float(tn / (tn + fp)),
        "pred_tcp_min": float(preds.min()),
        "pred_tcp_max": float(preds.max()),
        "cv_fold_aucs": cv_aucs,
    }


def evaluate_poisson(
    doses: np.ndarray,
    outcomes: np.ndarray,
    dose_label: str,
    outcome_definition: str,
    n_folds: int = 5,
) -> Dict[str, Any]:
    """Fit Poisson TCP and return quality metrics."""
    return evaluate_tcp_model(
        "poisson_tcp",
        lambda: PoissonTCPModel(d50_init=55.0, gamma50_init=1.5),
        doses,
        outcomes,
        dose_label,
        outcome_definition,
        k_params=2,
        param2_name="gamma50",
        n_folds=n_folds,
    )


def evaluate_probit(
    doses: np.ndarray,
    outcomes: np.ndarray,
    dose_label: str,
    outcome_definition: str,
    n_folds: int = 5,
) -> Dict[str, Any]:
    """Fit Probit TCP and return quality metrics."""
    return evaluate_tcp_model(
        "probit_tcp",
        lambda: ProbitTCPModel(d50_init=53.0, sigma_init=10.0),
        doses,
        outcomes,
        dose_label,
        outcome_definition,
        k_params=2,
        param2_name="sigma_gy",
        n_folds=n_folds,
    )


def evaluate_logistic(
    doses: np.ndarray,
    outcomes: np.ndarray,
    dose_label: str,
    outcome_definition: str,
    n_folds: int = 5,
) -> Dict[str, Any]:
    """Fit Logistic TCP and return quality metrics."""
    return evaluate_tcp_model(
        "logistic_tcp",
        lambda: LogisticTCPModel(d50_init=53.0, k_init=0.1),
        doses,
        outcomes,
        dose_label,
        outcome_definition,
        k_params=2,
        param2_name="k",
        n_folds=n_folds,
    )


def evaluate_eud(
    frame: pd.DataFrame,
    outcomes: np.ndarray,
    outcome_definition: str,
    n_folds: int = 5,
) -> Dict[str, Any]:
    """Fit gEUD TCP with best-a selection and return quality metrics."""
    geud_cols = EUDTCPModel.geud_columns_from_frame(frame)
    n = len(outcomes)
    p0 = float(outcomes.mean())
    nll_null = float(-np.sum(outcomes * np.log(p0) + (1.0 - outcomes) * np.log(1.0 - p0)))

    model = EUDTCPModel.fit_select_a(geud_cols, outcomes)
    selected_a = float(model.params_["geud_a"])
    doses = geud_cols[selected_a]
    preds = model.predict(doses)
    nll = model.nll_
    k_params = 3
    aic = 2 * k_params + 2 * nll
    bic = k_params * np.log(n) + 2 * nll
    lr_stat = 2 * ((-nll) - (-nll_null))
    lr_p = float(1.0 - stats.chi2.cdf(lr_stat, k_params))
    mcfadden = 1.0 - nll / nll_null

    yhat = (preds >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(outcomes, yhat).ravel()

    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=RANDOM_SEED)
    cv_aucs = []
    for train_idx, test_idx in skf.split(doses, outcomes):
        train_geud = {a: values[train_idx] for a, values in geud_cols.items()}
        cv_model = EUDTCPModel.fit_select_a(train_geud, outcomes[train_idx])
        test_doses = geud_cols[float(cv_model.params_["geud_a"])][test_idx]
        cv_preds = cv_model.predict(test_doses)
        cv_aucs.append(float(roc_auc_score(outcomes[test_idx], cv_preds)))

    return {
        "model": "eud_tcp",
        "dose_variable": f"gEUD_a{int(selected_a)}",
        "n": n,
        "outcome_definition": outcome_definition,
        "D50_gy": float(model.params_["D50_gy"]),
        "gamma50": float(model.params_["gamma50"]),
        "geud_a": selected_a,
        "param2_name": "gamma50",
        "nll": nll,
        "nll_null": nll_null,
        "aic": aic,
        "bic": bic,
        "aic_null": 2 * nll_null,
        "bic_null": 2 * nll_null,
        "lr_chi2": lr_stat,
        "lr_p_value": lr_p,
        "mcfadden_pseudo_r2": mcfadden,
        "roc_auc_insample": float(roc_auc_score(outcomes, preds)),
        "roc_auc_cv_mean": float(np.mean(cv_aucs)),
        "roc_auc_cv_std": float(np.std(cv_aucs)),
        "brier_score": float(brier_score_loss(outcomes, preds)),
        "brier_null": float(brier_score_loss(outcomes, np.full(n, p0))),
        "accuracy_at_0p5": float(accuracy_score(outcomes, yhat)),
        "sensitivity": float(tp / (tp + fn)),
        "specificity": float(tn / (tn + fp)),
        "pred_tcp_min": float(preds.min()),
        "pred_tcp_max": float(preds.max()),
        "cv_fold_aucs": cv_aucs,
    }


def list_figures() -> List[str]:
    """Return sorted figure filenames in figures/."""
    if not FIGURES_DIR.exists():
        return []
    return sorted(p.name for p in FIGURES_DIR.glob("*.png"))


def render_results_md(
    cohort_df: pd.DataFrame,
    survival_df: pd.DataFrame,
    assoc_df: pd.DataFrame,
    model_rows: List[Dict[str, Any]],
    figures: List[str],
    bootstrap_ci: Optional[pd.DataFrame] = None,
    model_comparison: Optional[pd.DataFrame] = None,
    cox_summary: Optional[pd.DataFrame] = None,
    km_logrank: Optional[pd.DataFrame] = None,
    c_index: Optional[float] = None,
    clinical_cox: Optional[pd.DataFrame] = None,
    who_ps: Optional[pd.DataFrame] = None,
    within_arm: Optional[pd.DataFrame] = None,
    hypo_volume: Optional[pd.DataFrame] = None,
    dose_heterogeneity: Optional[pd.DataFrame] = None,
    confounding_corr: Optional[pd.DataFrame] = None,
    unused_fields: Optional[pd.DataFrame] = None,
    tcp_verdict: Optional[Dict[str, str]] = None,
    rano_comparison: Optional[pd.DataFrame] = None,
    rano_assoc: Optional[pd.DataFrame] = None,
    rano_counts: Optional[pd.DataFrame] = None,
    within_arm_rano: Optional[pd.DataFrame] = None,
    cox_rano: Optional[pd.DataFrame] = None,
    rano_logistic: Optional[pd.DataFrame] = None,
    rano_logistic_boot: Optional[pd.DataFrame] = None,
    rano_pooled: Optional[pd.DataFrame] = None,
    rano_loocv_40: Optional[pd.DataFrame] = None,
    rano_loocv_pooled: Optional[pd.DataFrame] = None,
    pyro_comparison: Optional[pd.DataFrame] = None,
) -> str:
    """Render human-readable results markdown."""
    git_hash = _git_short_hash()
    today = date.today().isoformat()

    lines = [
        "# Current Results (auto-generated)",
        "",
        f"**Last updated:** {today}  ",
        f"**Git commit:** `{git_hash}`  ",
        f"**Regenerate:** `python -m src.reporting.update_results`",
        "",
        "> **Outcome caveat:** Primary TCP models still use an exploratory OS median-split proxy.",
        "> **CFB-GBM v3** adds RANO imaging response (non-PD vs PD at t1) — see §4b below.",
        "> RANO ≠ formal local control but is closer to tumor response than OS alone.",
        "",
        "---",
        "",
        "## 1. Cohort (verified from `modeling_table.csv`)",
        "",
        "| Metric | Value |",
        "|---|---|",
    ]
    display_keys = [
        "total_patients_tsv", "included_cohort_rules", "modeling_table_rows", "excluded",
        "fractionation_60gy_30fr", "fractionation_40gy_15fr",
        "age_median_yr", "os_median_wk", "sex_male", "sex_female",
        "dmean_median_gy", "gtv_volume_median_cc",
    ]
    lookup = {r["metric"]: r for _, r in cohort_df.iterrows()}
    labels = {
        "total_patients_tsv": "Total patients (TSV)",
        "included_cohort_rules": "Included (cohort rules)",
        "modeling_table_rows": "Modeling table rows",
        "excluded": "Excluded",
        "fractionation_60gy_30fr": "60 Gy / 30 fr",
        "fractionation_40gy_15fr": "40.05 Gy / 15 fr",
        "age_median_yr": "Median age",
        "os_median_wk": "Median OS",
        "sex_male": "Sex — male",
        "sex_female": "Sex — female",
        "dmean_median_gy": "Median GTV Dmean",
        "gtv_volume_median_cc": "Median GTV volume",
    }
    for key in display_keys:
        row = lookup[key]
        val = row["value"]
        unit = row["unit"]
        if unit == "n":
            disp = f"**{int(val)}**"
        elif unit == "yr":
            disp = f"**{val:.0f}** yr"
        elif unit == "wk":
            disp = f"**{val:.0f}** wk"
        elif unit == "Gy":
            disp = f"**{val:.2f}** Gy"
        elif unit == "cc":
            disp = f"**{val:.1f}** cc"
        else:
            disp = str(val)
        lines.append(f"| {labels[key]} | {disp} |")

    lines.extend(["", "## 2. Survival by fractionation (descriptive)", "", "| Scheme | n | Median OS (wk) | IQR (wk) |",
                  "|---|---:|---:|---:|"])
    for _, row in survival_df.iterrows():
        if row["scheme"] == "mann_whitney_p":
            continue
        lines.append(
            f"| {row['scheme']} | {int(row['n'])} | {row['os_median_wk']:.0f} | "
            f"{row['os_q25_wk']:.0f}–{row['os_q75_wk']:.0f} |"
        )
    p_row = survival_df[survival_df["scheme"] == "mann_whitney_p"]
    if not p_row.empty and "p_value" in p_row.columns:
        p_val = p_row["p_value"].iloc[0]
        lines.append(f"\nMann–Whitney U (60 Gy vs 40.05 Gy): **p = {p_val:.2e}**")

    lines.extend(["", "## 3. Dose–outcome association", "", "| Pair | r | p |", "|---|---:|---:|"])
    for _, row in assoc_df.iterrows():
        lines.append(f"| {row['pair']} | {row['value']:.4f} | {row['p_value']:.2e} |")

    lines.extend(["", "## 4. TCP model quality metrics", ""])
    for m in model_rows:
        param_row = f"| {m.get('param2_name', 'gamma50')} | {m.get(m.get('param2_name', 'gamma50'), m.get('gamma50', float('nan'))):.4f} |"
        block = [
            f"### {m['model']} — dose = `{m['dose_variable']}`",
            "",
            f"Outcome proxy: {m['outcome_definition']}",
            "",
            "| Metric | Value |",
            "|---|---|",
            f"| D50 (Gy) | {m['D50_gy']:.3f} |",
            param_row,
        ]
        if m.get("model") == "eud_tcp":
            block.append(f"| geud_a | {m['geud_a']:.1f} |")
        block.extend([
            f"| NLL (model / null) | {m['nll']:.2f} / {m['nll_null']:.2f} |",
            f"| AIC (model / null) | {m['aic']:.2f} / {m['aic_null']:.2f} |",
            f"| BIC (model / null) | {m['bic']:.2f} / {m['bic_null']:.2f} |",
            f"| LR test p-value | {m['lr_p_value']:.2e} |",
            f"| McFadden pseudo-R² | {m['mcfadden_pseudo_r2']:.4f} |",
            f"| ROC AUC (in-sample) | {m['roc_auc_insample']:.4f} |",
            f"| ROC AUC (5-fold CV) | {m['roc_auc_cv_mean']:.4f} ± {m['roc_auc_cv_std']:.4f} |",
            f"| Brier (model / null) | {m['brier_score']:.4f} / {m['brier_null']:.4f} |",
            f"| Accuracy @ TCP=0.5 | {m['accuracy_at_0p5']:.4f} |",
            f"| Sensitivity / Specificity | {m['sensitivity']:.3f} / {m['specificity']:.3f} |",
            "",
        ])
        lines.extend(block)

    if rano_comparison is not None and not rano_comparison.empty:
        n_rano = int(rano_comparison["n"].iloc[0])
        lines.extend(
            [
                "",
                "## 4b. RANO endpoint vs OS proxy (same patients, EQD2 models)",
                "",
                f"Subset with RANO t0→t1 labels: **n = {n_rano}**",
                "",
                "| Model | Endpoint | Event rate | AUC (in-sample) | AUC (5-fold CV) | LR p |",
                "|---|---|---:|---:|---:|---:|",
            ]
        )
        for _, row in rano_comparison.iterrows():
            ep = "OS median" if row["endpoint"] == "os_median_proxy" else "RANO non-PD"
            lines.append(
                f"| {row['model']} | {ep} | {row['event_rate']:.2f} | "
                f"{row['roc_auc_insample']:.4f} | {row['roc_auc_cv_mean']:.4f} ± {row['roc_auc_cv_std']:.4f} | "
                f"{row['lr_p_value']:.2e} |"
            )
        if rano_assoc is not None and not rano_assoc.empty:
            lines.append("")
            lines.append("Dose–outcome on RANO subset:")
            lines.append("")
            lines.append("| Pair | r | p |")
            lines.append("|---|---:|---:|")
            for _, row in rano_assoc.iterrows():
                lines.append(
                    f"| {row['pair']} | {row['point_biserial_r']:.4f} | {row['p_value']:.4f} |"
                )
        if rano_counts is not None and not rano_counts.empty:
            lines.append("")
            lines.append("RANO categories (modeling subset):")
            lines.append("")
            for _, row in rano_counts.iterrows():
                lines.append(f"- {row['category']}: {int(row['n'])}")
        lines.append("")

    if within_arm_rano is not None and not within_arm_rano.empty:
        lines.extend(
            [
                "",
                "## 4c. Within-arm DVH → RANO (Poisson TCP + Spearman)",
                "",
                "| Scheme | n | Metric | std | Spearman ρ | p | Poisson AUC | LR p |",
                "|---|---:|---|---:|---:|---:|---:|---:|",
            ]
        )
        for _, row in within_arm_rano.iterrows():
            auc = row["poisson_auc"]
            auc_str = f"{auc:.3f}" if pd.notna(auc) else "—"
            lr = row["poisson_lr_p"]
            lr_str = f"{lr:.3f}" if pd.notna(lr) else "—"
            feas = "" if row.get("metric_feasible", True) else " (low variance)"
            lines.append(
                f"| {row['scheme']} | {int(row['n_rano'])} | {row['metric']}{feas} | "
                f"{row['metric_std']:.3f} | {row['spearman_rho']:.3f} | {row['spearman_p']:.4f} | "
                f"{auc_str} | {lr_str} |"
            )
        lines.append("")

    if cox_rano is not None and not cox_rano.empty:
        c_idx = float(cox_rano["concordance_index"].iloc[0])
        n_cox = int(cox_rano["n_patients"].iloc[0])
        lines.extend(
            [
                "",
                f"## 4d. Cox OS ~ age + sex + WHO PS + EQD2 + RANO non-PD (n={n_cox})",
                "",
                f"**Concordance index:** {c_idx:.4f}",
                "",
                "| Covariate | HR | p |",
                "|---|---:|---:|",
            ]
        )
        for _, row in cox_rano.iterrows():
            lines.append(f"| {row['term']} | {row['hazard_ratio']:.3f} | {row['p_value']:.4f} |")
        lines.append("")

    if rano_logistic is not None and not rano_logistic.empty:
        lines.extend(
            [
                "",
                "## 4e. Multivariable logistic — RANO non-PD (40 Gy arm)",
                "",
                "| Model | n | AUC | Brier |",
                "|---|---:|---:|---:|",
            ]
        )
        for _, row in rano_logistic.iterrows():
            auc = row["auc"]
            auc_s = f"{auc:.3f}" if pd.notna(auc) else "—"
            lines.append(f"| {row['model']} | {int(row['n'])} | {auc_s} | {row['brier']:.3f} |")
        if rano_logistic_boot is not None and not rano_logistic_boot.empty:
            b = rano_logistic_boot.iloc[0]
            lines.append(
                f"\nBootstrap AUC (volume+age+PS): **{b['auc_mean']:.3f}** "
                f"[{b['auc_ci_lower']:.3f}, {b['auc_ci_upper']:.3f}] (n={int(b['n_bootstrap'])} resamples)"
            )
        lines.append("")

    if rano_pooled is not None and not rano_pooled.empty:
        lines.extend(
            [
                "",
                "## 4f. Pooled logistic — RANO non-PD (n≈137, volume + clinical + scheme)",
                "",
                "| Model | n | AUC | Brier |",
                "|---|---:|---:|---:|",
            ]
        )
        for _, row in rano_pooled.iterrows():
            lines.append(
                f"| {row['model']} | {int(row['n'])} | {row['auc_insample']:.3f} | {row['brier']:.3f} |"
            )
        lines.append("")

    if rano_loocv_40 is not None and not rano_loocv_40.empty:
        lines.extend(
            [
                "",
                "## 4g. LOOCV — 40 Gy arm (honest AUC)",
                "",
                "| Model | n | LOOCV AUC | LOOCV Brier |",
                "|---|---:|---:|---:|",
            ]
        )
        for _, row in rano_loocv_40.iterrows():
            lines.append(
                f"| {row['model']} | {int(row['n'])} | {row['loocv_auc']:.3f} | {row['loocv_brier']:.3f} |"
            )
        lines.append("")

    if rano_loocv_pooled is not None and not rano_loocv_pooled.empty:
        lines.extend(
            [
                "",
                "## 4h. LOOCV — pooled RANO models",
                "",
                "| Model | n | LOOCV AUC | LOOCV Brier |",
                "|---|---:|---:|---:|",
            ]
        )
        for _, row in rano_loocv_pooled.iterrows():
            lines.append(
                f"| {row['model']} | {int(row['n'])} | {row['loocv_auc']:.3f} | {row['loocv_brier']:.3f} |"
            )
        lines.append("")

    if pyro_comparison is not None and not pyro_comparison.empty:
        rano_pyro = pyro_comparison[pyro_comparison["endpoint"] == "RANO_non_PD"]
        lines.extend(
            [
                "",
                "## 4i. PyRadiomics vs DVH volume — RANO endpoint (t1gd GTV t0)",
                "",
                "| Model | n | AUC | Brier |",
                "|---|---:|---:|---:|",
            ]
        )
        for _, row in rano_pyro.iterrows():
            lines.append(
                f"| {row['model']} | {int(row['n'])} | {row['auc_insample']:.3f} | {row['brier']:.3f} |"
            )
        lines.append("")

    if bootstrap_ci is not None and not bootstrap_ci.empty:
        lines.extend(["", "## 5. Bootstrap 95% CI (Poisson TCP, EQD2)", ""])
        lines.append("| Parameter | Estimate | 95% CI | Bootstrap SD |")
        lines.append("|---|---:|---|---:|")
        for _, row in bootstrap_ci.iterrows():
            lines.append(
                f"| {row['parameter']} | {row['estimate']:.3f} | "
                f"[{row['ci_lower']:.3f}, {row['ci_upper']:.3f}] | {row['bootstrap_std']:.3f} |"
            )
        lines.append("")

        lines.append("")

    if model_comparison is not None and not model_comparison.empty:
        lines.extend(["", "## 6. Four-model comparison (EQD2 / gEUD, sorted by AIC)", ""])
        lines.append("| Model | k | AIC | BIC | ROC AUC | Brier | HL p-value |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|")
        for _, row in model_comparison.iterrows():
            lines.append(
                f"| {row['model']} | {int(row['k_params'])} | {row['aic']:.2f} | "
                f"{row['bic']:.2f} | {row['roc_auc']:.4f} | {row['brier_score']:.4f} | "
                f"{row['hl_p_value']:.4f} |"
            )
        lines.append("")

        lines.append("")

    section_num = 7
    if cox_summary is not None and not cox_summary.empty:
        lines.extend(["", f"## {section_num}. Cox PH regression (OS ~ EQD2 + Dmean + age + sex)", ""])
        if c_index is not None:
            lines.append(f"**Concordance index:** {c_index:.4f}")
        lines.append("")
        lines.append("| Covariate | HR | 95% CI | p |")
        lines.append("|---|---:|---|---:|")
        for _, row in cox_summary.iterrows():
            lines.append(
                f"| {row['term']} | {row['hazard_ratio']:.3f} | "
                f"[{row['hr_ci_lower']:.3f}, {row['hr_ci_upper']:.3f}] | {row['p']:.4f} |"
            )
        if km_logrank is not None and not km_logrank.empty:
            lr = km_logrank.iloc[0]
            lines.append(
                f"\nKM log-rank ({lr['comparison']}): χ² = {lr['logrank_chi2']:.2f}, "
                f"p = {lr['logrank_p_value']:.2e}"
            )
        lines.append("")
        section_num += 1

    if clinical_cox is not None and not clinical_cox.empty:
        c_idx = clinical_cox["concordance_index"].iloc[0]
        lines.extend(["", f"## {section_num}. Clinical prognosis (Cox: OS ~ age + sex + WHO PS + scheme)", ""])
        lines.append(f"**Concordance index:** {c_idx:.4f}")
        lines.append("")
        lines.append("| Covariate | HR | p |")
        lines.append("|---|---:|---:|")
        for _, row in clinical_cox.iterrows():
            lines.append(f"| {row['term']} | {row['hazard_ratio']:.3f} | {row['p_value']:.4f} |")
        lines.append("")
        section_num += 1

    if who_ps is not None and not who_ps.empty:
        lines.extend(["", f"## {section_num}. OS by WHO performance status", ""])
        lines.append("| WHO PS | n | Median OS (wk) | IQR (wk) |")
        lines.append("|---:|---:|---:|---:|")
        for _, row in who_ps.iterrows():
            lines.append(
                f"| {int(row['who_status'])} | {int(row['n'])} | {row['os_median_wk']:.0f} | "
                f"{row['os_q25_wk']:.0f}–{row['os_q75_wk']:.0f} |"
            )
        kw_p = who_ps.attrs.get("kruskal_p")
        if kw_p is not None:
            lines.append(f"\nKruskal–Wallis: **p = {kw_p:.2e}**")
        lines.append("")
        section_num += 1

    if within_arm is not None and not within_arm.empty:
        lines.extend(["", f"## {section_num}. Within-arm DVH vs OS (Spearman)", ""])
        lines.append("| Scheme | n | Metric | ρ | p |")
        lines.append("|---|---:|---|---:|---:|")
        for _, row in within_arm.iterrows():
            lines.append(
                f"| {row['scheme']} | {int(row['n'])} | {row['metric']} | "
                f"{row['spearman_rho']:.3f} | {row['p_value']:.4f} |"
            )
        lines.append("")
        section_num += 1

    if hypo_volume is not None and not hypo_volume.empty:
        vol = hypo_volume[hypo_volume["term"] == "volume_cc"]
        if not vol.empty:
            v = vol.iloc[0]
            lines.extend(["", f"## {section_num}. Hypofractionated arm: Cox OS ~ volume + covariates", ""])
            lines.append(f"n = {int(hypo_volume['n_patients'].iloc[0])}, C-index = {v.get('concordance_index', hypo_volume['concordance_index'].iloc[0]):.3f}")
            lines.append(f"\nGTV volume HR = **{v['hazard_ratio']:.4f}**/cc, p = **{v['p_value']:.4f}** (exploratory)")
            lines.append("")
            section_num += 1

    if dose_heterogeneity is not None and not dose_heterogeneity.empty:
        lines.extend(["", f"## {section_num}. TCP feasibility: dose heterogeneity within arm", ""])
        lines.append("| Scheme | n | Dmean SD (Gy) | min–max (Gy) | DVH-TCP feasible? |")
        lines.append("|---|---:|---:|---|---|")
        for _, row in dose_heterogeneity.iterrows():
            feasible = "yes" if row["tcp_dvh_feasible"] else "no"
            lines.append(
                f"| {row['rt_dose_gy']:.2f} Gy / {int(row['n_fractions'])} fr | {int(row['n'])} | "
                f"{row['dmean_std_gy']:.2f} | {row['dmean_min_gy']:.1f}–{row['dmean_max_gy']:.1f} | {feasible} |"
            )
        lines.append("")
        section_num += 1

    if confounding_corr is not None and not confounding_corr.empty:
        lines.extend(["", f"## {section_num}. Confounding correlations", ""])
        lines.append("| Pair | r | p |")
        lines.append("|---|---:|---:|")
        for _, row in confounding_corr.iterrows():
            lines.append(f"| {row['pair']} | {row['value']:.4f} | {row['p_value']:.2e} |")
        lines.append("")
        section_num += 1

    if unused_fields is not None and not unused_fields.empty:
        lines.extend(["", f"## {section_num}. Unused TSV fields (association with OS)", ""])
        lines.append("| Field | In modeling table? | ρ vs OS | p |")
        lines.append("|---|---|---:|---:|")
        for _, row in unused_fields.iterrows():
            in_table = "yes" if row["in_modeling_table"] else "no"
            lines.append(
                f"| {row['field']} | {in_table} | {row['spearman_rho_vs_os']:.3f} | {row['p_value']:.2e} |"
            )
        lines.append("")
        section_num += 1

    if tcp_verdict:
        lines.extend(["", f"## {section_num}. TCP feasibility verdict", ""])
        for key, text in tcp_verdict.items():
            lines.append(f"- **{key.replace('_', ' ')}:** {text}")
        lines.append("")
        section_num += 1

    fig_section = f"## {section_num}. Figures"
    lines.extend([fig_section, ""])
    if figures:
        for fig in figures:
            lines.append(f"- [`figures/{fig}`](../figures/{fig})")
    else:
        lines.append("_No PNG figures found in `figures/`._")

    lines.extend([
        "",
        "---",
        "",
        "## Interpretation snapshot",
        "",
        "| Question | Current answer |",
        "|---|---|",
        "| Data pipeline complete? | Yes — 190-patient modeling table with 21 DVH metrics |",
        "| Strong clinical signal? | Yes — 60 vs 40 Gy OS (p ≈ 3×10⁻⁶); WHO PS (p ≈ 1.6×10⁻⁴); Cox scheme HR≈0.54 |",
        "| DVH-TCP within standard arm? | No — 60 Gy Dmean SD = 0.28 Gy (below 1 Gy threshold) |",
        "| TCP model beats null? | Yes — LR p ≈ 3×10⁻⁶ (Poisson, EQD2) but confounded by scheme/age |",
        "| Good discrimination (AUC ≥ 0.7)? | No — in-sample AUC ≈ 0.68, CV ≈ 0.68 ± 0.10 |",
        "| True TCP validation? | Partial — RANO non-PD available (v3); still not formal LC |",
        "| RANO improves AUC vs OS on same n? | No — pooled RANO AUC ≈ 0.43 vs OS ≈ 0.62 (n=137) |",
        "| Within-arm volume → RANO (40 Gy)? | Yes — Poisson AUC ≈ 0.83, LR p ≈ 0.037 (n=34); LOOCV AUC ≈ 0.74 |",
        "| Pooled volume + scheme → RANO? | Yes — in-sample AUC ≈ 0.72 (n=137); LOOCV ≈ 0.64 |",
        "| PyRadiomics beats DVH volume for RANO? | Exploratory — top-5 radiomics AUC ≈ 0.78 vs volume 0.71 |",
        "| Within-arm volume → RANO (60 Gy)? | Exploratory — AUC ≈ 0.66, Spearman p ≈ 0.019 (n=96) |",
        "| Calibration fixes ranking? | No — Platt scaling does not change AUC on same data |",
        "",
    ])
    return "\n".join(lines) + "\n"


# Module-level frame reference removed — use update_results() locals only


def update_results() -> Path:
    """
    Recompute all metrics and write report files.

    Returns
    -------
    pathlib.Path
        Path to ``reports/RESULTS.md``.
    """
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    METRICS_DIR.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    cohort = pd.read_csv(DATA_PROCESSED / "cohort.csv")
    median_os = float(frame["survival_weeks"].median())
    os_outcome_def = f"OS >= median ({median_os:.0f} wk)"
    outcomes = (frame["survival_weeks"] >= median_os).astype(float).to_numpy()

    cohort_df = cohort_metrics(frame, cohort)
    survival_df = survival_by_fractionation(frame)
    assoc_df = association_metrics(frame, outcomes)

    doses = frame["eqd2_gy"].to_numpy()
    model_rows = [
        evaluate_poisson(doses, outcomes, "eqd2_gy", os_outcome_def),
        evaluate_logistic(doses, outcomes, "eqd2_gy", os_outcome_def),
        evaluate_probit(doses, outcomes, "eqd2_gy", os_outcome_def),
        evaluate_eud(frame, outcomes, os_outcome_def),
        evaluate_poisson(frame["Dmean_gy"].to_numpy(), outcomes, "Dmean_gy", os_outcome_def),
    ]
    model_csv_rows = [{k: v for k, v in m.items() if k != "cv_fold_aucs"} for m in model_rows]
    model_df = pd.DataFrame(model_csv_rows)

    _save_csv(cohort_df, "cohort_summary.csv")
    _save_csv(survival_df, "survival_by_fractionation.csv")
    _save_csv(assoc_df, "dose_outcome_association.csv")
    _save_csv(model_df, "tcp_model_quality.csv")

    comparison_df, _ = run_model_comparison(doses, outcomes, frame=frame)
    _save_csv(comparison_df, "tcp_model_comparison.csv")

    km_logrank, cox_summary, _, _ = run_survival_analysis(frame)
    _save_csv(km_logrank, "km_logrank_eqd2.csv")
    _save_csv(cox_summary.drop(columns=["concordance_index"], errors="ignore"), "cox_ph_summary.csv")
    c_index = float(cox_summary["concordance_index"].iloc[0])

    bootstrap_ci = bootstrap_tcp_params(
        PoissonTCPModel,
        doses,
        outcomes,
        n_bootstrap=1000,
        model_factory=lambda: PoissonTCPModel(d50_init=55.0, gamma50_init=1.5),
    )
    _save_csv(bootstrap_ci, "poisson_tcp_bootstrap_ci.csv")

    manifest = {
        "generated": date.today().isoformat(),
        "git_hash": _git_short_hash(),
        "n_modeling": len(frame),
        "median_os_wk": median_os,
        "models": model_csv_rows,
        "figures": list_figures(),
    }
    (METRICS_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    dvh_cols = [c for c in SCALAR_METRIC_KEYS if c in frame.columns]
    dvh_summary = frame[dvh_cols].describe().T[["count", "mean", "std", "min", "50%", "max"]]
    dvh_summary.columns = ["n", "mean", "std", "min", "median", "max"]
    _save_csv(dvh_summary.reset_index().rename(columns={"index": "metric"}), "dvh_scalars_summary.csv")

    clinical_cox, who_ps, within_arm, hypo_volume = run_stratified_analysis(
        frame, METRICS_DIR, figure_path=FIGURES_DIR / "04_clinical_prognosis.png"
    )
    audit_tables = run_confounding_audit(frame, METRICS_DIR)
    tcp_verdict = tcp_feasibility_summary(frame)
    rano_comparison, rano_counts, rano_assoc = run_rano_tcp_comparison(
        frame, METRICS_DIR, figure_path=FIGURES_DIR / "05_rano_vs_os_tcp_auc.png"
    )
    within_arm_rano, cox_rano = run_within_arm_rano_analysis(
        frame, METRICS_DIR, figure_path=FIGURES_DIR / "06_within_arm_rano_tcp.png"
    )
    rano_logistic, _, rano_logistic_boot, _ = run_rano_multivariable_40gy(frame, METRICS_DIR)
    rano_suite = run_rano_prediction_suite(frame, METRICS_DIR)
    run_volume_validation(METRICS_DIR)

    figures = list_figures()
    md = render_results_md(
        cohort_df,
        survival_df,
        assoc_df,
        model_rows,
        figures,
        bootstrap_ci,
        comparison_df,
        cox_summary,
        km_logrank,
        c_index,
        clinical_cox=clinical_cox,
        who_ps=who_ps,
        within_arm=within_arm,
        hypo_volume=hypo_volume,
        dose_heterogeneity=audit_tables["dose_heterogeneity"],
        confounding_corr=audit_tables["confounding_correlations"],
        unused_fields=audit_tables["unused_clinical_fields"],
        tcp_verdict=tcp_verdict,
        rano_comparison=rano_comparison,
        rano_assoc=rano_assoc,
        rano_counts=rano_counts,
        within_arm_rano=within_arm_rano,
        cox_rano=cox_rano,
        rano_logistic=rano_logistic,
        rano_logistic_boot=rano_logistic_boot,
        rano_pooled=rano_suite["pooled"],
        rano_loocv_40=rano_suite["loocv_40"],
        rano_loocv_pooled=rano_suite["loocv_pooled"],
        pyro_comparison=rano_suite["pyro_comparison"],
    )
    RESULTS_MD.write_text(md)
    print(f"Wrote {RESULTS_MD}")
    print(f"Metrics CSVs in {METRICS_DIR}/")
    return RESULTS_MD


def main() -> None:
    update_results()


if __name__ == "__main__":
    main()
