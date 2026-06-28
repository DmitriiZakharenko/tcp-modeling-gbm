"""
Model comparison utilities: AIC/BIC table, ROC, calibration, Hosmer–Lemeshow.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score

from src.config import FIGURES_DIR, RANDOM_SEED
from src.models.base_model import TCPModel
from src.models.eud_tcp import EUDTCPModel
from src.models.logistic_tcp import LogisticTCPModel
from src.models.poisson_tcp import PoissonTCPModel
from src.models.probit_tcp import ProbitTCPModel

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


def hosmer_lemeshow(
    outcomes: np.ndarray,
    predicted: np.ndarray,
    n_bins: int = 10,
) -> Dict[str, float]:
    """
    Hosmer–Lemeshow goodness-of-fit test for binary outcomes.

    Parameters
    ----------
    outcomes : np.ndarray
        Observed binary outcomes.
    predicted : np.ndarray
        Predicted probabilities in (0, 1).
    n_bins : int, optional
        Number of risk groups (default 10).

    Returns
    -------
    dict
        Keys: ``chi2``, ``df``, ``p_value``, ``n_bins_used``.
    """
    from scipy.stats import chi2

    outcomes = np.asarray(outcomes, dtype=float)
    predicted = np.clip(np.asarray(predicted, dtype=float), 1e-9, 1 - 1e-9)

    order = np.argsort(predicted)
    outcomes = outcomes[order]
    predicted = predicted[order]
    splits = np.array_split(np.arange(len(predicted)), n_bins)

    chi2_stat = 0.0
    bins_used = 0
    for idx in splits:
        if idx.size == 0:
            continue
        obs = outcomes[idx]
        pred = predicted[idx]
        n_g = idx.size
        o1 = obs.sum()
        e1 = pred.sum()
        o0 = n_g - o1
        e0 = n_g - e1
        if e1 <= 0 or e0 <= 0:
            continue
        chi2_stat += (o1 - e1) ** 2 / e1 + (o0 - e0) ** 2 / e0
        bins_used += 1

    df = max(bins_used - 2, 1)
    p_value = float(1.0 - chi2.cdf(chi2_stat, df))
    return {"chi2": float(chi2_stat), "df": float(df), "p_value": p_value, "n_bins_used": float(bins_used)}


def _null_nll(outcomes: np.ndarray) -> float:
    p0 = float(outcomes.mean())
    return float(-np.sum(outcomes * np.log(p0) + (1.0 - outcomes) * np.log(1.0 - p0)))


def evaluate_fitted_model(
    name: str,
    model: TCPModel,
    doses: np.ndarray,
    outcomes: np.ndarray,
    k_params: int,
) -> Dict[str, Any]:
    """Compute comparison metrics for an already-fitted model."""
    preds = model.predict(doses)
    n = len(outcomes)
    nll = model.nll_
    nll_null = _null_nll(outcomes)
    hl = hosmer_lemeshow(outcomes, preds)
    return {
        "model": name,
        "k_params": k_params,
        "nll": nll,
        "aic": 2 * k_params + 2 * nll,
        "bic": k_params * np.log(n) + 2 * nll,
        "log_likelihood": -nll,
        "roc_auc": float(roc_auc_score(outcomes, preds)),
        "brier_score": float(brier_score_loss(outcomes, preds)),
        "hl_chi2": hl["chi2"],
        "hl_df": hl["df"],
        "hl_p_value": hl["p_value"],
        "predictions": preds,
    }


def fit_all_tcp_models(
    doses: np.ndarray,
    outcomes: np.ndarray,
    frame: Optional[pd.DataFrame] = None,
) -> List[Dict[str, Any]]:
    """
    Fit all four TCP models and return comparison metric dicts.

    Parameters
    ----------
    doses : np.ndarray
        EQD2 doses in Gy.
    outcomes : np.ndarray
        Binary outcomes.
    frame : pd.DataFrame, optional
        Required for EUD model (gEUD columns).

    Returns
    -------
    list of dict
        One entry per model with metrics and predictions.
    """
    results: List[Dict[str, Any]] = []

    poisson = PoissonTCPModel(d50_init=55.0, gamma50_init=1.5)
    poisson.fit(doses, outcomes)
    results.append(evaluate_fitted_model("poisson_tcp", poisson, doses, outcomes, k_params=2))

    logistic = LogisticTCPModel(d50_init=53.0, k_init=0.1)
    logistic.fit(doses, outcomes)
    results.append(evaluate_fitted_model("logistic_tcp", logistic, doses, outcomes, k_params=2))

    probit = ProbitTCPModel(d50_init=53.0, sigma_init=10.0)
    probit.fit(doses, outcomes)
    results.append(evaluate_fitted_model("probit_tcp", probit, doses, outcomes, k_params=2))

    if frame is not None:
        geud_cols = EUDTCPModel.geud_columns_from_frame(frame)
        eud = EUDTCPModel.fit_select_a(geud_cols, outcomes)
        eud_doses = geud_cols[float(eud.params_["geud_a"])]
        results.append(evaluate_fitted_model("eud_tcp", eud, eud_doses, outcomes, k_params=3))

    return results


def comparison_table(results: Sequence[Dict[str, Any]]) -> pd.DataFrame:
    """
    Build a summary comparison table from ``fit_all_tcp_models`` output.

    Parameters
    ----------
    results : sequence of dict
        Output of ``fit_all_tcp_models``.

    Returns
    -------
    pd.DataFrame
        One row per model with AIC, BIC, AUC, Brier, HL p-value.
    """
    rows = []
    for r in results:
        rows.append(
            {
                "model": r["model"],
                "k_params": r["k_params"],
                "log_likelihood": r["log_likelihood"],
                "aic": r["aic"],
                "bic": r["bic"],
                "roc_auc": r["roc_auc"],
                "brier_score": r["brier_score"],
                "hl_chi2": r["hl_chi2"],
                "hl_p_value": r["hl_p_value"],
            }
        )
    return pd.DataFrame(rows).sort_values("aic").reset_index(drop=True)


def plot_calibration(
    results: Sequence[Dict[str, Any]],
    outcomes: np.ndarray,
    save_path: Optional[Path] = None,
    n_bins: int = 10,
) -> plt.Figure:
    """
    Plot calibration curves for multiple fitted models.

    Parameters
    ----------
    results : sequence of dict
        Must include ``model`` and ``predictions`` keys.
    outcomes : np.ndarray
        Binary outcomes.
    save_path : pathlib.Path, optional
        If set, save figure at 300 dpi.
    n_bins : int, optional
        Number of calibration bins.

    Returns
    -------
    matplotlib.figure.Figure
    """
    outcomes = np.asarray(outcomes, dtype=float)
    colors = ["#2c7bb6", "#d7191c", "#fdae61", "#1a9641"]

    with plt.rc_context(DEFAULT_RC_PARAMS):
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        ax.plot([0, 1], [0, 1], "k--", lw=1, label="Perfect calibration")

        for i, r in enumerate(results):
            preds = np.clip(r["predictions"], 1e-9, 1 - 1e-9)
            bins = np.linspace(0.0, 1.0, n_bins + 1)
            bin_centers = []
            obs_rates = []
            for b in range(n_bins):
                mask = (preds >= bins[b]) & (preds < bins[b + 1] if b < n_bins - 1 else preds <= bins[b + 1])
                if mask.sum() == 0:
                    continue
                bin_centers.append(preds[mask].mean())
                obs_rates.append(outcomes[mask].mean())
            ax.plot(
                bin_centers,
                obs_rates,
                "o-",
                lw=1.8,
                ms=5,
                color=colors[i % len(colors)],
                label=r["model"],
            )

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Predicted TCP")
        ax.set_ylabel("Observed outcome rate")
        ax.set_title("Calibration plot (OS ≥ median proxy)")
        ax.legend(fontsize=9, frameon=False)
        ax.grid(True, alpha=0.25, lw=0.6)
        plt.tight_layout()

        if save_path is not None:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path)

    return fig


def run_model_comparison(
    doses: np.ndarray,
    outcomes: np.ndarray,
    frame: Optional[pd.DataFrame] = None,
    figure_path: Optional[Path] = None,
) -> Tuple[pd.DataFrame, plt.Figure]:
    """
    Fit all models, build comparison table, and plot calibration.

    Parameters
    ----------
    doses : np.ndarray
        EQD2 doses.
    outcomes : np.ndarray
        Binary outcomes.
    frame : pd.DataFrame, optional
        Modeling table for EUD model.
    figure_path : pathlib.Path, optional
        Calibration plot output path.

    Returns
    -------
    tuple
        ``(comparison_table, figure)``
    """
    results = fit_all_tcp_models(doses, outcomes, frame=frame)
    table = comparison_table(results)
    fig = plot_calibration(
        results,
        outcomes,
        save_path=figure_path or FIGURES_DIR / "03_model_calibration.png",
    )
    return table, fig


def main() -> None:
    """Run four-model comparison on modeling cohort."""
    from src.config import DATA_PROCESSED

    frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    doses = frame["eqd2_gy"].to_numpy()
    median_os = frame["survival_weeks"].median()
    outcomes = (frame["survival_weeks"] >= median_os).astype(float).to_numpy()

    table, _ = run_model_comparison(doses, outcomes, frame=frame)
    print(f"Cohort: n={len(frame)}, median OS split at {median_os:.0f} wk")
    print(f"Random seed (project default): {RANDOM_SEED}")
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
