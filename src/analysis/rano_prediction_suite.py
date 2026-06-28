"""
Pooled RANO prediction, PyRadiomics comparison, LOOCV, and volume×scheme interaction.

Extends the 40 Gy arm analysis to the full RANO subset (n≈137).
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
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

from src.analysis.rano_multivariable import (
    DEFAULT_RC_PARAMS,
    HYPO_DOSE_GY,
    _fit_logistic,
    _prep_features,
    hypofractionated_rano_cohort,
)
from src.config import FIGURES_DIR, PYRADIOMICS_TSV, RANDOM_SEED

RADIOMICS_SEQUENCE = "t1gd"
PYRO_TOP_N = 5


def rano_cohort(frame: pd.DataFrame) -> pd.DataFrame:
    """All modeling patients with RANO t0→t1 label."""
    return frame[frame["rano_controlled_t1"].notna()].copy()


def _with_scheme(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["scheme_60gy"] = (out["rt_dose_gy"] >= 59.0).astype(float)
    return out


def _fit_on_subset(
    sub: pd.DataFrame,
    covariates: List[str],
    outcome_col: str = "rano_controlled_t1",
) -> Dict[str, object]:
    data = sub[covariates + [outcome_col]].dropna()
    y = data[outcome_col].astype(float).to_numpy()
    x = data[covariates].astype(float).to_numpy()
    return _fit_logistic(x, y, covariates)


def compare_pooled_rano_models(frame: pd.DataFrame) -> pd.DataFrame:
    """Logistic models for RANO non-PD on the pooled RANO subset."""
    sub = _with_scheme(rano_cohort(frame))
    specs = [
        ("eqd2_only", ["eqd2_gy"]),
        ("volume_only", ["volume_cc"]),
        ("volume_age_ps", ["volume_cc", "age", "who_status"]),
        ("volume_age_ps_scheme", ["volume_cc", "age", "who_status", "scheme_60gy"]),
        (
            "volume_age_ps_scheme_interaction",
            ["volume_cc", "age", "who_status", "scheme_60gy", "volume_x_scheme"],
        ),
    ]
    sub["volume_x_scheme"] = sub["volume_cc"] * sub["scheme_60gy"]

    rows: List[Dict[str, object]] = []
    for name, covs in specs:
        fit = _fit_on_subset(sub, covs)
        rows.append(
            {
                "model": name,
                "n": fit["n"],
                "covariates": "+".join(covs),
                "auc_insample": fit["auc"],
                "brier": fit["brier"],
                "event_rate": fit["event_rate"],
            }
        )
    return pd.DataFrame(rows)


def loocv_logistic_auc(
    sub: pd.DataFrame,
    covariates: List[str],
    outcome_col: str = "rano_controlled_t1",
) -> pd.DataFrame:
    """Leave-one-out cross-validated AUC for logistic regression."""
    data = sub[covariates + [outcome_col]].dropna().reset_index(drop=True)
    y = data[outcome_col].astype(float).to_numpy()
    x = data[covariates].astype(float).to_numpy()
    n = len(y)
    if n < 5 or len(np.unique(y)) < 2:
        return pd.DataFrame([{"n": n, "loocv_auc": np.nan, "loocv_brier": np.nan}])

    loo = LeaveOneOut()
    probs = np.zeros(n)
    for train_idx, test_idx in loo.split(x):
        scaler = StandardScaler()
        xs_train = scaler.fit_transform(x[train_idx])
        xs_test = scaler.transform(x[test_idx])
        model = LogisticRegression(max_iter=2000, random_state=RANDOM_SEED)
        model.fit(xs_train, y[train_idx])
        probs[test_idx] = model.predict_proba(xs_test)[:, 1]

    return pd.DataFrame(
        [
            {
                "n": n,
                "covariates": "+".join(covariates),
                "loocv_auc": float(roc_auc_score(y, probs)),
                "loocv_brier": float(brier_score_loss(y, probs)),
            }
        ]
    )


def loocv_40gy_models(frame: pd.DataFrame) -> pd.DataFrame:
    """LOOCV for 40 Gy arm logistic models."""
    sub = hypofractionated_rano_cohort(frame)
    rows = []
    for name, covs in [
        ("volume_only", ["volume_cc"]),
        ("volume_age_ps", ["volume_cc", "age", "who_status"]),
    ]:
        cv = loocv_logistic_auc(sub, covs)
        cv.insert(0, "model", name)
        rows.append(cv)
    return pd.concat(rows, ignore_index=True)


def load_pyradiomics_t1gd() -> pd.DataFrame:
    """Baseline GTV t1gd PyRadiomics features keyed by patient_id."""
    pyro = pd.read_csv(PYRADIOMICS_TSV, sep="\t")
    mask = (
        (pyro["Temporality"] == "t0")
        & (pyro["Label name"] == "GTV")
        & (pyro["Sequence"] == RADIOMICS_SEQUENCE)
    )
    sub = pyro.loc[mask].copy()
    sub["patient_id"] = sub["Patient"].astype(int)
    feature_cols = [c for c in sub.columns if c.startswith("original_")]
    return sub[["patient_id"] + feature_cols]


def _top_radiomics_features(
    merged: pd.DataFrame,
    feature_cols: List[str],
    top_n: int = PYRO_TOP_N,
) -> List[str]:
    """Rank radiomics features by univariate AUC for RANO non-PD."""
    y = merged["rano_controlled_t1"].astype(float)
    scores: List[Tuple[str, float]] = []
    for col in feature_cols:
        x = merged[col].astype(float)
        valid = x.notna() & y.notna()
        if valid.sum() < 20 or len(np.unique(y[valid])) < 2:
            continue
        try:
            auc = float(roc_auc_score(y[valid], x[valid]))
            auc_adj = max(auc, 1 - auc)
            scores.append((col, auc_adj))
        except ValueError:
            continue
    scores.sort(key=lambda t: t[1], reverse=True)
    return [s[0] for s in scores[:top_n]]


def compare_pyradiomics_vs_volume(frame: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compare DVH volume vs PyRadiomics for RANO and OS median-split endpoints.

    Returns
    -------
    model_comparison, top_features_table
    """
    pyro = load_pyradiomics_t1gd()
    sub = rano_cohort(frame)
    merged = _with_scheme(sub.merge(pyro, on="patient_id", how="inner"))
    feature_cols = [c for c in pyro.columns if c.startswith("original_")]
    top_feats = _top_radiomics_features(merged, feature_cols)

    median_os = float(frame["survival_weeks"].median())
    merged["os_median_proxy"] = (merged["survival_weeks"] >= median_os).astype(float)

    rows: List[Dict[str, object]] = []
    for endpoint, outcome_col in [
        ("RANO_non_PD", "rano_controlled_t1"),
        ("OS_median_proxy", "os_median_proxy"),
    ]:
        base_covs = ["volume_cc", "age", "who_status", "scheme_60gy"]
        data = merged[base_covs + [outcome_col]].dropna()
        if len(data) < 20:
            continue

        for model_name, covs in [
            ("dvh_volume_only", ["volume_cc"]),
            ("dvh_volume_clinical", base_covs),
            ("pyro_mesh_volume", ["original_shape_MeshVolume"]),
            (f"pyro_top{len(top_feats)}", top_feats),
            (f"pyro_top{len(top_feats)}_clinical", top_feats + ["age", "who_status", "scheme_60gy"]),
            ("dvh_volume_plus_pyro_top1", ["volume_cc"] + top_feats[:1]),
        ]:
            avail = [c for c in covs if c in merged.columns]
            fit = _fit_on_subset(merged, avail, outcome_col)
            rows.append(
                {
                    "endpoint": endpoint,
                    "model": model_name,
                    "n": fit["n"],
                    "n_features": len(avail),
                    "auc_insample": fit["auc"],
                    "brier": fit["brier"],
                }
            )

    top_table = pd.DataFrame(
        [{"rank": i + 1, "feature": f} for i, f in enumerate(top_feats)]
    )
    return pd.DataFrame(rows), top_table


def plot_pooled_rano_roc(frame: pd.DataFrame, save_path: Optional[Path] = None) -> plt.Figure:
    """ROC for pooled volume models vs EQD2 baseline."""
    from sklearn.metrics import RocCurveDisplay

    sub = _with_scheme(rano_cohort(frame))
    sub["volume_x_scheme"] = sub["volume_cc"] * sub["scheme_60gy"]
    data = sub[
        ["volume_cc", "age", "who_status", "scheme_60gy", "volume_x_scheme", "eqd2_gy", "rano_controlled_t1"]
    ].dropna()
    y = data["rano_controlled_t1"].astype(float).to_numpy()

    with plt.rc_context(DEFAULT_RC_PARAMS):
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        specs = [
            (["eqd2_gy"], "EQD2 only", "#969696"),
            (["volume_cc", "age", "who_status", "scheme_60gy"], "Volume + clinical + scheme", "#2c7bb6"),
            (
                ["volume_cc", "age", "who_status", "scheme_60gy", "volume_x_scheme"],
                "Volume + interaction",
                "#d7191c",
            ),
        ]
        for covs, label, color in specs:
            fit = _fit_on_subset(data.assign(rano_controlled_t1=y), covs)
            RocCurveDisplay.from_predictions(
                y, fit["probs"], ax=ax, name=f"{label} (AUC={fit['auc']:.2f})"
            )
            ax.get_lines()[-1].set_color(color)
        ax.plot([0, 1], [0, 1], "k--", lw=0.8)
        ax.set_title(f"Pooled RANO non-PD models (n={len(y)})")
        ax.legend(frameon=False, fontsize=8, loc="lower right")
        plt.tight_layout()
        if save_path is not None:
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path)
    return fig


def plot_pyradiomics_comparison(
    comparison: pd.DataFrame,
    save_path: Optional[Path] = None,
) -> plt.Figure:
    """Bar chart of AUC by model for RANO endpoint."""
    sub = comparison[comparison["endpoint"] == "RANO_non_PD"].copy()
    with plt.rc_context(DEFAULT_RC_PARAMS):
        fig, ax = plt.subplots(figsize=(8, 4.5))
        labels = sub["model"].str.replace("_", " ")
        ax.barh(labels, sub["auc_insample"], color="#4575b4", alpha=0.85)
        ax.set_xlim(0, 1)
        ax.set_xlabel("In-sample ROC AUC")
        ax.set_title("RANO non-PD: DVH volume vs PyRadiomics (t1gd GTV t0)")
        ax.axvline(0.5, color="k", ls="--", lw=0.6, alpha=0.5)
        for i, (_, row) in enumerate(sub.iterrows()):
            ax.text(row["auc_insample"] + 0.01, i, f"{row['auc_insample']:.2f}", va="center", fontsize=8)
        plt.tight_layout()
        if save_path is not None:
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path)
    return fig


def run_rano_prediction_suite(
    frame: pd.DataFrame,
    metrics_dir: Path,
    figure_dir: Path | None = None,
) -> Dict[str, pd.DataFrame]:
    """Run pooled, LOOCV, PyRadiomics analyses and save outputs."""
    metrics_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = figure_dir or FIGURES_DIR

    pooled = compare_pooled_rano_models(frame)
    loocv_40 = loocv_40gy_models(frame)
    sub_pooled = _with_scheme(rano_cohort(frame))
    sub_pooled["volume_x_scheme"] = sub_pooled["volume_cc"] * sub_pooled["scheme_60gy"]
    loocv_rows = [
        loocv_logistic_auc(sub_pooled, ["volume_cc", "age", "who_status", "scheme_60gy"]).assign(
            model="volume_age_ps_scheme"
        ),
        loocv_logistic_auc(
            sub_pooled,
            ["volume_cc", "age", "who_status", "scheme_60gy", "volume_x_scheme"],
        ).assign(model="volume_age_ps_scheme_interaction"),
    ]
    loocv_pooled = pd.concat(loocv_rows, ignore_index=True)

    pyro_cmp, pyro_top = compare_pyradiomics_vs_volume(frame)

    pooled.to_csv(metrics_dir / "rano_pooled_logistic_comparison.csv", index=False)
    loocv_40.to_csv(metrics_dir / "rano_loocv_40gy.csv", index=False)
    loocv_pooled.to_csv(metrics_dir / "rano_loocv_pooled.csv", index=False)
    pyro_cmp.to_csv(metrics_dir / "pyradiomics_vs_volume_comparison.csv", index=False)
    pyro_top.to_csv(metrics_dir / "pyradiomics_top_features_rano.csv", index=False)

    plot_pooled_rano_roc(frame, fig_dir / "08_pooled_rano_roc.png")
    if not pyro_cmp.empty:
        plot_pyradiomics_comparison(pyro_cmp, fig_dir / "08_pyradiomics_vs_volume_auc.png")

    return {
        "pooled": pooled,
        "loocv_40": loocv_40,
        "loocv_pooled": loocv_pooled,
        "pyro_comparison": pyro_cmp,
        "pyro_top": pyro_top,
    }


def main() -> None:
    from src.config import DATA_PROCESSED, REPORTS_DIR

    frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    out = run_rano_prediction_suite(frame, REPORTS_DIR / "metrics")
    print("Pooled models:\n", out["pooled"].to_string(index=False))
    print("\nLOOCV 40 Gy:\n", out["loocv_40"].to_string(index=False))
    print("\nLOOCV pooled:\n", out["loocv_pooled"].to_string(index=False))
    print("\nPyRadiomics comparison:\n", out["pyro_comparison"].to_string(index=False))


if __name__ == "__main__":
    main()
