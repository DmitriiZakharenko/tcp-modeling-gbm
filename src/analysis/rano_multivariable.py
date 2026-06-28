"""
Multivariable logistic models for RANO non-PD within the 40 Gy hypofractionated arm.

Also validates DVH GTV volume against RANO segmentation volumes (TSV).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.preprocessing import StandardScaler

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

HYPO_DOSE_GY = 40.05


def hypofractionated_rano_cohort(frame: pd.DataFrame) -> pd.DataFrame:
    """40.05 Gy / 15 fr patients with RANO t0→t1 label."""
    return frame[
        (frame["rt_dose_gy"] == HYPO_DOSE_GY) & frame["rano_controlled_t1"].notna()
    ].copy()


def _prep_features(sub: pd.DataFrame, covariates: List[str]) -> Tuple[np.ndarray, np.ndarray]:
    y = sub["rano_controlled_t1"].astype(float).to_numpy()
    x = sub[covariates].astype(float).to_numpy()
    return x, y


def _fit_logistic(
    x: np.ndarray,
    y: np.ndarray,
    covariate_names: List[str],
) -> Dict[str, object]:
    scaler = StandardScaler()
    xs = scaler.fit_transform(x)
    model = LogisticRegression(max_iter=2000, random_state=RANDOM_SEED)
    model.fit(xs, y)
    probs = model.predict_proba(xs)[:, 1]
    auc = float(roc_auc_score(y, probs)) if len(np.unique(y)) > 1 else np.nan
    brier = float(brier_score_loss(y, probs))

    coef_rows = []
    for name, coef in zip(covariate_names, model.coef_.ravel()):
        coef_rows.append({"term": name, "coef_stdized": float(coef), "odds_ratio": float(np.exp(coef))})
    coef_rows.insert(
        0,
        {"term": "intercept", "coef_stdized": float(model.intercept_[0]), "odds_ratio": float(np.exp(model.intercept_[0]))},
    )

    return {
        "model": model,
        "scaler": scaler,
        "covariates": covariate_names,
        "n": len(y),
        "event_rate": float(y.mean()),
        "auc": auc,
        "brier": brier,
        "coef_table": pd.DataFrame(coef_rows),
        "probs": probs,
    }


def bootstrap_multivariable_auc(
    sub: pd.DataFrame,
    covariates: List[str],
    n_bootstrap: int = 1000,
) -> pd.DataFrame:
    """Bootstrap 95% CI for multivariable logistic AUC."""
    rng = np.random.default_rng(RANDOM_SEED)
    x, y = _prep_features(sub, covariates)
    n = len(y)
    aucs: List[float] = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            continue
        try:
            fit = _fit_logistic(x[idx], y[idx], covariates)
            if not np.isnan(fit["auc"]):
                aucs.append(float(fit["auc"]))
        except ValueError:
            continue
    aucs_arr = np.array(aucs)
    return pd.DataFrame(
        [
            {
                "n_bootstrap": len(aucs_arr),
                "auc_mean": float(aucs_arr.mean()) if len(aucs_arr) else np.nan,
                "auc_ci_lower": float(np.percentile(aucs_arr, 2.5)) if len(aucs_arr) else np.nan,
                "auc_ci_upper": float(np.percentile(aucs_arr, 97.5)) if len(aucs_arr) else np.nan,
            }
        ]
    )


def compare_logistic_models(frame: pd.DataFrame) -> pd.DataFrame:
    """Compare intercept-only, volume-only, and multivariable models (40 Gy arm)."""
    sub = hypofractionated_rano_cohort(frame)
    rows: List[Dict[str, object]] = []
    specs = [
        ("intercept_only", []),
        ("volume_only", ["volume_cc"]),
        ("volume_age_ps", ["volume_cc", "age", "who_status"]),
    ]
    fits: Dict[str, Dict] = {}
    for name, covs in specs:
        if not covs:
            p0 = float(sub["rano_controlled_t1"].mean())
            y = sub["rano_controlled_t1"].astype(float).to_numpy()
            nll = float(-np.sum(y * np.log(p0) + (1 - y) * np.log(1 - p0)))
            rows.append(
                {
                    "model": name,
                    "n": len(sub),
                    "k_params": 0,
                    "covariates": "",
                    "auc": np.nan,
                    "brier": float(brier_score_loss(y, np.full(len(y), p0))),
                    "nll": nll,
                }
            )
            continue
        data = sub[covs + ["rano_controlled_t1"]].dropna()
        x, y = _prep_features(data, covs)
        fit = _fit_logistic(x, y, covs)
        fits[name] = fit
        probs = np.clip(fit["probs"], 1e-9, 1 - 1e-9)
        nll = float(-np.sum(y * np.log(probs) + (1 - y) * np.log(1 - probs)))
        rows.append(
            {
                "model": name,
                "n": fit["n"],
                "k_params": len(covs),
                "covariates": "+".join(covs),
                "auc": fit["auc"],
                "brier": fit["brier"],
                "nll": nll,
            }
        )

    table = pd.DataFrame(rows)
    if "volume_only" in fits and "volume_age_ps" in fits:
        nll1 = float(table.loc[table["model"] == "volume_only", "nll"].iloc[0])
        nll2 = float(table.loc[table["model"] == "volume_age_ps", "nll"].iloc[0])
        lr = 2 * (nll1 - nll2)
        table.attrs["lr_volume_vs_multi"] = float(lr)
        table.attrs["lr_p_volume_vs_multi"] = float(1 - stats.chi2.cdf(lr, df=2))
    return table


def multivariable_coef_table(frame: pd.DataFrame) -> pd.DataFrame:
    """Standardized coefficients and ORs for volume + age + WHO PS."""
    sub = hypofractionated_rano_cohort(frame)
    covs = ["volume_cc", "age", "who_status"]
    data = sub[covs + ["rano_controlled_t1"]].dropna()
    x, y = _prep_features(data, covs)
    fit = _fit_logistic(x, y, covs)
    out = fit["coef_table"].copy()
    out["model"] = "volume_age_ps"
    out["n"] = fit["n"]
    out["auc"] = fit["auc"]
    return out


def validate_dvh_vs_rano_volumes(frame: pd.DataFrame) -> pd.DataFrame:
    """Compare DVH volume_cc with RANO size_t0_cm3 / size_t1_cm3."""
    rows = []
    for dvh_col, rano_col, label in [
        ("volume_cc", "size_t0_cm3", "t0_DVH_vs_RANO"),
        ("volume_cc", "size_t1_cm3", "t0_DVH_vs_RANO_t1size"),
    ]:
        if dvh_col not in frame.columns or rano_col not in frame.columns:
            continue
        sub = frame[[dvh_col, rano_col]].dropna()
        if len(sub) < 10:
            continue
        rho, p = stats.spearmanr(sub[dvh_col], sub[rano_col])
        r_pearson, p_pearson = stats.pearsonr(sub[dvh_col], sub[rano_col])
        rows.append(
            {
                "comparison": label,
                "n": len(sub),
                "spearman_rho": float(rho),
                "spearman_p": float(p),
                "pearson_r": float(r_pearson),
                "pearson_p": float(p_pearson),
                "mean_abs_pct_diff": float(
                    (np.abs(sub[dvh_col] - sub[rano_col]) / sub[rano_col].clip(lower=0.1)).mean() * 100
                ),
            }
        )
    return pd.DataFrame(rows)


def plot_volume_validation(frame: pd.DataFrame, save_path: Optional[Path] = None) -> plt.Figure:
    """Scatter DVH volume vs RANO size_t0 for 40 Gy RANO patients."""
    sub = hypofractionated_rano_cohort(frame)[["volume_cc", "size_t0_cm3"]].dropna()
    with plt.rc_context(DEFAULT_RC_PARAMS):
        fig, ax = plt.subplots(figsize=(5.5, 5))
        ax.scatter(sub["size_t0_cm3"], sub["volume_cc"], alpha=0.75, c="#2c7bb6", edgecolors="white", s=45)
        lim = max(sub["size_t0_cm3"].max(), sub["volume_cc"].max()) * 1.05
        ax.plot([0, lim], [0, lim], "k--", lw=0.8, alpha=0.6)
        rho, p = stats.spearmanr(sub["size_t0_cm3"], sub["volume_cc"])
        ax.set_xlabel("RANO size t0 (cm³)")
        ax.set_ylabel("DVH GTV volume t0 (cm³)")
        ax.set_title(f"40 Gy arm volume validation (n={len(sub)})\nSpearman ρ={rho:.2f}, p={p:.3f}")
        plt.tight_layout()
        if save_path is not None:
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path)
    return fig


def plot_multivariable_roc(frame: pd.DataFrame, save_path: Optional[Path] = None) -> plt.Figure:
    """ROC curves for univariate vs multivariable logistic (40 Gy)."""
    from sklearn.metrics import RocCurveDisplay

    sub = hypofractionated_rano_cohort(frame)
    covs = ["volume_cc", "age", "who_status"]
    data = sub[covs + ["rano_controlled_t1"]].dropna()
    y = data["rano_controlled_t1"].astype(float).to_numpy()

    with plt.rc_context(DEFAULT_RC_PARAMS):
        fig, ax = plt.subplots(figsize=(6, 5))
        for covs_use, label, color in [
            (["volume_cc"], "Volume only", "#2c7bb6"),
            (covs, "Volume + age + WHO PS", "#d7191c"),
        ]:
            xi = data[covs_use].astype(float).to_numpy()
            fit = _fit_logistic(xi, y, covs_use)
            RocCurveDisplay.from_predictions(y, fit["probs"], ax=ax, name=f"{label} (AUC={fit['auc']:.2f})")
            ax.get_lines()[-1].set_color(color)
        ax.plot([0, 1], [0, 1], "k--", lw=0.8)
        ax.set_title(f"RANO non-PD at t1 — 40 Gy arm (n={len(y)})")
        ax.legend(frameon=False, fontsize=9)
        plt.tight_layout()
        if save_path is not None:
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path)
    return fig


def run_rano_multivariable_40gy(
    frame: pd.DataFrame,
    metrics_dir: Path,
    figure_dir: Path | None = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run full 40 Gy multivariable analysis and save CSVs + figures."""
    metrics_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = figure_dir or FIGURES_DIR

    comparison = compare_logistic_models(frame)
    coefs = multivariable_coef_table(frame)
    boot = bootstrap_multivariable_auc(
        hypofractionated_rano_cohort(frame)[["volume_cc", "age", "who_status", "rano_controlled_t1"]].dropna(),
        ["volume_cc", "age", "who_status"],
    )
    validation = validate_dvh_vs_rano_volumes(frame)

    comparison.to_csv(metrics_dir / "rano_logistic_40gy_model_comparison.csv", index=False)
    coefs.to_csv(metrics_dir / "rano_logistic_40gy_coefficients.csv", index=False)
    boot.to_csv(metrics_dir / "rano_logistic_40gy_bootstrap_auc.csv", index=False)
    validation.to_csv(metrics_dir / "rano_volume_validation.csv", index=False)

    plot_volume_validation(frame, fig_dir / "07_rano_volume_validation_40gy.png")
    plot_multivariable_roc(frame, fig_dir / "07_rano_logistic_roc_40gy.png")

    return comparison, coefs, boot, validation


def main() -> None:
    from src.config import DATA_PROCESSED, REPORTS_DIR

    frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    comp, coefs, boot, val = run_rano_multivariable_40gy(frame, REPORTS_DIR / "metrics")
    print("Model comparison:\n", comp.to_string(index=False))
    print("\nCoefficients:\n", coefs.to_string(index=False))
    print("\nBootstrap AUC:\n", boot.to_string(index=False))
    print("\nVolume validation:\n", val.to_string(index=False))


if __name__ == "__main__":
    main()
