"""
Bootstrap confidence intervals for fitted TCP models.

Resamples patients with replacement (default n=1000, seed=42) and
refits any ``TCPModel`` subclass to obtain percentile CIs.
"""

from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Type

import numpy as np
import pandas as pd

from src.config import RANDOM_SEED
from src.models.base_model import TCPModel


def bootstrap_tcp_params(
    model_cls: Type[TCPModel],
    doses: np.ndarray,
    outcomes: np.ndarray,
    n_bootstrap: int = 1000,
    seed: int = RANDOM_SEED,
    model_factory: Optional[Callable[[], TCPModel]] = None,
    fit_kwargs: Optional[Dict[str, Any]] = None,
    percentiles: Tuple[float, float] = (2.5, 97.5),
) -> pd.DataFrame:
    """
    Bootstrap percentile confidence intervals for TCP model parameters.

    Parameters
    ----------
    model_cls : type
        ``TCPModel`` subclass to refit on each bootstrap sample.
    doses : np.ndarray, shape (n_patients,)
        Dose or gEUD values in Gy.
    outcomes : np.ndarray, shape (n_patients,)
        Binary outcomes.
    n_bootstrap : int, optional
        Number of bootstrap resamples.
    seed : int, optional
        Random seed for reproducibility.
    model_factory : callable, optional
        Returns a fresh model instance; defaults to ``model_cls()``.
    fit_kwargs : dict, optional
        Extra keyword arguments passed to ``fit`` (e.g. ``geud_a``).
    percentiles : tuple of float, optional
        Lower and upper percentile for CI (default 95%).

    Returns
    -------
    pd.DataFrame
        Columns: parameter, estimate, ci_lower, ci_upper, bootstrap_std.

    Raises
    ------
    ValueError
        If inputs are empty or lengths mismatch.
    RuntimeError
        If the point-estimate fit fails.
    """
    doses = np.asarray(doses, dtype=float)
    outcomes = np.asarray(outcomes, dtype=float)
    fit_kwargs = fit_kwargs or {}

    if doses.size == 0:
        raise ValueError("doses must not be empty")
    if outcomes.size != doses.size:
        raise ValueError("doses and outcomes must have the same length")

    factory = model_factory or model_cls
    point_model = factory()
    point_model.fit(doses, outcomes, **fit_kwargs)
    param_names = list(point_model.params_.keys())

    rng = np.random.default_rng(seed)
    n = doses.size
    boot_values: Dict[str, List[float]] = {name: [] for name in param_names}
    failures = 0

    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        sample_doses = doses[idx]
        sample_outcomes = outcomes[idx]
        if np.unique(sample_outcomes).size < 2:
            failures += 1
            continue
        try:
            boot_model = factory()
            boot_model.fit(sample_doses, sample_outcomes, **fit_kwargs)
        except (RuntimeError, ValueError):
            failures += 1
            continue
        for name in param_names:
            boot_values[name].append(float(boot_model.params_[name]))

    if not boot_values[param_names[0]]:
        raise RuntimeError("All bootstrap fits failed")

    rows = []
    for name in param_names:
        samples = np.array(boot_values[name], dtype=float)
        rows.append(
            {
                "parameter": name,
                "estimate": float(point_model.params_[name]),
                "ci_lower": float(np.percentile(samples, percentiles[0])),
                "ci_upper": float(np.percentile(samples, percentiles[1])),
                "bootstrap_std": float(samples.std(ddof=1)),
                "n_bootstrap_success": int(samples.size),
                "n_bootstrap_failures": failures,
            }
        )

    return pd.DataFrame(rows)


def main() -> None:
    """Bootstrap CI for Poisson TCP on EQD2 (median-OS proxy)."""
    from src.config import DATA_PROCESSED
    from src.models.poisson_tcp import PoissonTCPModel

    frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    doses = frame["eqd2_gy"].to_numpy()
    median_os = frame["survival_weeks"].median()
    outcomes = (frame["survival_weeks"] >= median_os).astype(float).to_numpy()

    ci = bootstrap_tcp_params(
        PoissonTCPModel,
        doses,
        outcomes,
        n_bootstrap=1000,
        model_factory=lambda: PoissonTCPModel(d50_init=55.0, gamma50_init=1.5),
    )

    print(f"Cohort: n={len(frame)}, median OS split at {median_os:.0f} wk")
    print(f"Bootstrap: n=1000, seed={RANDOM_SEED}")
    print(ci.to_string(index=False))


if __name__ == "__main__":
    main()
