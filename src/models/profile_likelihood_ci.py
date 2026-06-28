"""
Profile-likelihood confidence intervals for Poisson TCP parameters.

Complements bootstrap CIs (Part IV): fixes one parameter at the MLE and scans
the other, using the chi-square(1) threshold for 95% intervals.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from src.models.poisson_tcp import PoissonTCPModel, poisson_tcp

_CHI2_95_1DF = 3.841459  # 95% CI, 1 degree of freedom (NLL shift = 0.5 * this)


def _nll_at(
    doses: np.ndarray,
    outcomes: np.ndarray,
    d50: float,
    gamma50: float,
) -> float:
    p = np.clip(poisson_tcp(doses, d50, gamma50), 1e-9, 1 - 1e-9)
    return float(-np.sum(outcomes * np.log(p) + (1 - outcomes) * np.log(1 - p)))


def profile_likelihood_ci_poisson(
    doses: np.ndarray,
    outcomes: np.ndarray,
    alpha: float = 0.05,
    n_grid: int = 120,
    d50_bounds: Tuple[float, float] = (30.0, 80.0),
    gamma_bounds: Tuple[float, float] = (0.5, 8.0),
) -> pd.DataFrame:
    """
    Compute profile-likelihood 95% CIs for D50 and gamma50 (Poisson TCP).

    Parameters
    ----------
    doses, outcomes : arrays
        Training data for MLE.
    alpha : float
        Two-sided significance level (default 0.05 -> 95% CI).
    n_grid : int
        Grid points for bracketing each profile curve.
    d50_bounds, gamma_bounds : tuple
        Search ranges when profiling the complementary parameter.

    Returns
    -------
    pd.DataFrame
        Columns: parameter, estimate, ci_lower, ci_upper, method.
    """
    doses = np.asarray(doses, dtype=float)
    outcomes = np.asarray(outcomes, dtype=float)
    chi2_cut = -np.log(alpha) * 2  # two-sided; 0.05 -> 3.841

    base = PoissonTCPModel(d50_init=55.0, gamma50_init=2.0)
    base.fit(doses, outcomes)
    d50_mle = base.params_["D50_gy"]
    g50_mle = base.params_["gamma50"]
    nll_min = base.nll_
    threshold = nll_min + 0.5 * chi2_cut

    def profile_d50(d50: float) -> float:
        def obj(g: float) -> float:
            return _nll_at(doses, outcomes, d50, g)

        res = minimize_scalar(
            obj,
            bounds=gamma_bounds,
            method="bounded",
        )
        return float(res.fun)

    def profile_gamma(gamma: float) -> float:
        def obj(d: float) -> float:
            return _nll_at(doses, outcomes, d, gamma)

        res = minimize_scalar(
            obj,
            bounds=d50_bounds,
            method="bounded",
        )
        return float(res.fun)

    d50_grid = np.linspace(d50_bounds[0], d50_bounds[1], n_grid)
    g50_grid = np.linspace(gamma_bounds[0], gamma_bounds[1], n_grid)

    d50_prof = np.array([profile_d50(d) for d in d50_grid])
    g50_prof = np.array([profile_gamma(g) for g in g50_grid])

    d50_in = d50_grid[d50_prof <= threshold]
    g50_in = g50_grid[g50_prof <= threshold]

    rows: List[Dict[str, object]] = [
        {
            "parameter": "D50_gy",
            "estimate": d50_mle,
            "ci_lower": float(d50_in.min()) if d50_in.size else np.nan,
            "ci_upper": float(d50_in.max()) if d50_in.size else np.nan,
            "method": "profile_likelihood",
        },
        {
            "parameter": "gamma50",
            "estimate": g50_mle,
            "ci_lower": float(g50_in.min()) if g50_in.size else np.nan,
            "ci_upper": float(g50_in.max()) if g50_in.size else np.nan,
            "method": "profile_likelihood",
        },
    ]
    return pd.DataFrame(rows)


def main() -> None:
    from src.config import DATA_PROCESSED

    frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    doses = frame["eqd2_gy"].to_numpy()
    median_os = frame["survival_weeks"].median()
    outcomes = (frame["survival_weeks"] >= median_os).astype(float).to_numpy()
    table = profile_likelihood_ci_poisson(doses, outcomes)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
