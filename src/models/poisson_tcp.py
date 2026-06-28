"""
Poisson TCP dose-response model (Webb–Nahum parameterisation).

TCP(D) = exp( −ln(2) · exp( γ50 · (1 − D/D50) ) )
"""

from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from src.models.base_model import TCPModel

LN2 = np.log(2.0)
_EPS = 1e-9


def poisson_tcp(doses: np.ndarray, d50: float, gamma50: float) -> np.ndarray:
    """
    Compute Poisson TCP for an array of doses.

    Parameters
    ----------
    doses : np.ndarray
        Dose values in Gy (EQD2-corrected).
    d50 : float
        Dose yielding TCP = 0.5 (Gy).
    gamma50 : float
        Normalised dose–response gradient at D50.

    Returns
    -------
    np.ndarray
        Predicted TCP in (0, 1).
    """
    doses = np.asarray(doses, dtype=float)
    exponent = gamma50 * (1.0 - doses / d50)
    return np.exp(-LN2 * np.exp(exponent))


class PoissonTCPModel(TCPModel):
    """
    Poisson TCP model fitted by maximum likelihood on binary outcomes.

    Parameters
    ----------
    d50_init : float, optional
        Initial guess for D50 (Gy).
    gamma50_init : float, optional
        Initial guess for γ50.
    bounds : tuple of tuple, optional
        ``((d50_min, d50_max), (gamma50_min, gamma50_max))`` for MLE.
    """

    PARAM_NAMES: Tuple[str, ...] = ("D50_gy", "gamma50")

    def __init__(
        self,
        d50_init: float = 50.0,
        gamma50_init: float = 2.0,
        bounds: Optional[Tuple[Tuple[float, float], Tuple[float, float]]] = None,
    ) -> None:
        super().__init__()
        self.d50_init = d50_init
        self.gamma50_init = gamma50_init
        self.bounds = bounds or ((1.0, 100.0), (0.01, 20.0))

    def predict(self, doses: np.ndarray) -> np.ndarray:
        """
        Predict TCP for given doses using fitted parameters.

        Parameters
        ----------
        doses : np.ndarray
            Dose values in Gy.

        Returns
        -------
        np.ndarray
            Predicted TCP values.

        Raises
        ------
        RuntimeError
            If the model has not been fitted.
        """
        self._check_fitted()
        return poisson_tcp(doses, self.params_["D50_gy"], self.params_["gamma50"])

    def log_likelihood(
        self,
        params: np.ndarray,
        doses: np.ndarray,
        outcomes: np.ndarray,
    ) -> float:
        """
        Negative log-likelihood for Bernoulli outcomes.

        Parameters
        ----------
        params : np.ndarray
            ``[D50, gamma50]``.
        doses : np.ndarray
            Dose values in Gy.
        outcomes : np.ndarray
            Binary outcomes (1 = success, 0 = failure).

        Returns
        -------
        float
            Negative log-likelihood (minimisation target).
        """
        d50, gamma50 = float(params[0]), float(params[1])
        if d50 <= 0.0 or gamma50 <= 0.0:
            return np.inf

        tcp = poisson_tcp(doses, d50, gamma50)
        tcp = np.clip(tcp, _EPS, 1.0 - _EPS)
        outcomes = np.asarray(outcomes, dtype=float)
        nll = -np.sum(outcomes * np.log(tcp) + (1.0 - outcomes) * np.log(1.0 - tcp))
        return float(nll)

    def fit(self, doses: np.ndarray, outcomes: np.ndarray) -> "PoissonTCPModel":
        """
        Fit D50 and γ50 by maximum likelihood.

        Parameters
        ----------
        doses : np.ndarray, shape (n_patients,)
            Dose values in Gy (EQD2-corrected).
        outcomes : np.ndarray, shape (n_patients,)
            Binary outcomes.

        Returns
        -------
        self

        Raises
        ------
        ValueError
            If inputs are empty or outcomes are not binary.
        """
        doses = np.asarray(doses, dtype=float)
        outcomes = np.asarray(outcomes, dtype=float)

        if doses.size == 0:
            raise ValueError("doses must not be empty")
        if outcomes.size != doses.size:
            raise ValueError("doses and outcomes must have the same length")
        if not np.all(np.isin(outcomes, (0.0, 1.0))):
            raise ValueError("outcomes must be binary (0 or 1)")

        x0 = np.array([self.d50_init, self.gamma50_init], dtype=float)
        result = minimize(
            self.log_likelihood,
            x0,
            args=(doses, outcomes),
            method="L-BFGS-B",
            bounds=self.bounds,
        )

        if not result.success:
            raise RuntimeError(f"Poisson TCP MLE failed: {result.message}")

        self.params_ = {"D50_gy": float(result.x[0]), "gamma50": float(result.x[1])}
        self.fitted_ = True
        self.nll_ = float(result.fun)
        self.n_obs_ = int(doses.size)
        return self

    def summary(self) -> pd.DataFrame:
        """
        Return fitted parameter estimates.

        Returns
        -------
        pd.DataFrame
            Columns: parameter, estimate.

        Raises
        ------
        RuntimeError
            If called before fit().
        """
        self._check_fitted()
        return pd.DataFrame(
            {"parameter": list(self.PARAM_NAMES), "estimate": [self.params_[k] for k in self.PARAM_NAMES]}
        )


def main() -> None:
    """Fit Poisson TCP on modeling cohort EQD2 with median-OS binary proxy (sanity check)."""
    from src.config import DATA_PROCESSED, RANDOM_SEED

    frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    doses = frame["eqd2_gy"].to_numpy()
    median_os = frame["survival_weeks"].median()
    outcomes = (frame["survival_weeks"] >= median_os).astype(float).to_numpy()

    model = PoissonTCPModel(d50_init=55.0, gamma50_init=1.5)
    model.fit(doses, outcomes)

    print(f"Cohort: n={len(frame)}, median OS split at {median_os:.0f} wk")
    print(model.summary().to_string(index=False))
    print(f"Negative log-likelihood: {model.nll_:.2f}")
    print(f"Random seed (project default): {RANDOM_SEED}")


if __name__ == "__main__":
    main()
