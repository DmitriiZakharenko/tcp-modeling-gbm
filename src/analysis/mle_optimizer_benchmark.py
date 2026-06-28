"""
Compare scipy.optimize methods for Poisson TCP MLE (Part IV documentation).

Default production code uses L-BFGS-B (bounded, fast, stable on this smooth
likelihood). This module benchmarks alternatives and records whether they
recover the same D50/gamma50.
"""

from __future__ import annotations

import time
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize

from src.models.poisson_tcp import PoissonTCPModel, poisson_tcp

_EPS = 1e-9


def _neg_log_lik(params: np.ndarray, doses: np.ndarray, outcomes: np.ndarray) -> float:
    d50, gamma = float(params[0]), float(params[1])
    if d50 <= 0 or gamma <= 0:
        return 1e12
    p = np.clip(poisson_tcp(doses, d50, gamma), _EPS, 1 - _EPS)
    return float(-np.sum(outcomes * np.log(p) + (1 - outcomes) * np.log(1 - p)))


def _fit_with_method(
    doses: np.ndarray,
    outcomes: np.ndarray,
    method: str,
    x0: np.ndarray,
    bounds: Tuple[Tuple[float, float], Tuple[float, float]],
) -> Dict[str, Any]:
    t0 = time.perf_counter()
    if method == "differential_evolution":
        result = differential_evolution(
            _neg_log_lik,
            bounds=bounds,
            args=(doses, outcomes),
            seed=42,
            maxiter=300,
            polish=True,
            tol=1e-6,
        )
    else:
        kw: Dict[str, Any] = {"method": method, "args": (doses, outcomes)}
        if method in ("L-BFGS-B", "TNC", "SLSQP"):
            kw["bounds"] = bounds
        result = minimize(_neg_log_lik, x0, **kw)
    elapsed = time.perf_counter() - t0
    return {
        "method": method,
        "success": bool(result.success),
        "D50_gy": float(result.x[0]),
        "gamma50": float(result.x[1]),
        "nll": float(result.fun),
        "seconds": elapsed,
        "message": str(getattr(result, "message", "")),
    }


def run_mle_optimizer_benchmark(
    doses: np.ndarray,
    outcomes: np.ndarray,
    metrics_dir: Any = None,
) -> pd.DataFrame:
    """
    Benchmark MLE optimizers on Poisson TCP for the same dose/outcome vectors.

    Returns
    -------
    pd.DataFrame
        One row per optimizer; reference row from production ``PoissonTCPModel``.
    """
    doses = np.asarray(doses, dtype=float)
    outcomes = np.asarray(outcomes, dtype=float)
    bounds = ((30.0, 80.0), (0.5, 8.0))
    x0 = np.array([55.0, 2.0])

    ref = PoissonTCPModel(d50_init=55.0, gamma50_init=2.0)
    ref.fit(doses, outcomes)
    ref_d50 = ref.params_["D50_gy"]
    ref_g50 = ref.params_["gamma50"]

    methods = [
        "L-BFGS-B",
        "TNC",
        "SLSQP",
        "Powell",
        "Nelder-Mead",
        "differential_evolution",
    ]
    rows: List[Dict[str, Any]] = []
    for method in methods:
        try:
            row = _fit_with_method(doses, outcomes, method, x0, bounds)
        except Exception as exc:  # pragma: no cover
            row = {
                "method": method,
                "success": False,
                "D50_gy": np.nan,
                "gamma50": np.nan,
                "nll": np.nan,
                "seconds": np.nan,
                "message": str(exc),
            }
        row["delta_D50_vs_lbfgs"] = (
            abs(row["D50_gy"] - ref_d50) if np.isfinite(row["D50_gy"]) else np.nan
        )
        row["delta_gamma_vs_lbfgs"] = (
            abs(row["gamma50"] - ref_g50) if np.isfinite(row["gamma50"]) else np.nan
        )
        rows.append(row)

    df = pd.DataFrame(rows)
    if metrics_dir is not None:
        from pathlib import Path

        out = Path(metrics_dir) / "mle_optimizer_benchmark.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, index=False)
    return df


def main() -> None:
    from src.config import DATA_PROCESSED, REPORTS_DIR

    frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    doses = frame["eqd2_gy"].to_numpy()
    outcomes = (frame["survival_weeks"] >= frame["survival_weeks"].median()).astype(float).to_numpy()
    df = run_mle_optimizer_benchmark(doses, outcomes, REPORTS_DIR / "metrics")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
