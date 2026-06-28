"""
Logistic TCP dose-response model.

TCP(D) = 1 / (1 + exp(−k · (D − D50)))
"""

from typing import Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from src.models.base_model import TCPModel

_EPS = 1e-9


def logistic_tcp(doses: np.ndarray, d50: float, k: float) -> np.ndarray:
    """
    Compute logistic TCP for an array of doses.

    Parameters
    ----------
    doses : np.ndarray
        Dose values in Gy (EQD2-corrected).
    d50 : float
        Dose yielding TCP = 0.5 (Gy).
    k : float
        Logistic steepness (Gy⁻¹).

    Returns
    -------
    np.ndarray
        Predicted TCP in (0, 1).
    """
    doses = np.asarray(doses, dtype=float)
    return 1.0 / (1.0 + np.exp(-k * (doses - d50)))


class LogisticTCPModel(TCPModel):
    """
    Logistic TCP model fitted by maximum likelihood on binary outcomes.

    Parameters
    ----------
    d50_init : float, optional
        Initial guess for D50 (Gy).
    k_init : float, optional
        Initial guess for steepness k (Gy⁻¹).
    bounds : tuple of tuple, optional
        ``((d50_min, d50_max), (k_min, k_max))`` for MLE.
    """

    PARAM_NAMES: Tuple[str, ...] = ("D50_gy", "k")

    def __init__(
        self,
        d50_init: float = 50.0,
        k_init: float = 0.1,
        bounds: Optional[Tuple[Tuple[float, float], Tuple[float, float]]] = None,
    ) -> None:
        super().__init__()
        self.d50_init = d50_init
        self.k_init = k_init
        self.bounds = bounds or ((1.0, 100.0), (0.001, 5.0))

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
        return logistic_tcp(doses, self.params_["D50_gy"], self.params_["k"])

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
            ``[D50, k]``.
        doses : np.ndarray
            Dose values in Gy.
        outcomes : np.ndarray
            Binary outcomes (1 = success, 0 = failure).

        Returns
        -------
        float
            Negative log-likelihood (minimisation target).
        """
        d50, k = float(params[0]), float(params[1])
        if d50 <= 0.0 or k <= 0.0:
            return np.inf

        tcp = logistic_tcp(doses, d50, k)
        tcp = np.clip(tcp, _EPS, 1.0 - _EPS)
        outcomes = np.asarray(outcomes, dtype=float)
        nll = -np.sum(outcomes * np.log(tcp) + (1.0 - outcomes) * np.log(1.0 - tcp))
        return float(nll)

    def fit(self, doses: np.ndarray, outcomes: np.ndarray) -> "LogisticTCPModel":
        """
        Fit D50 and k by maximum likelihood.

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
        RuntimeError
            If optimisation fails.
        """
        doses = np.asarray(doses, dtype=float)
        outcomes = np.asarray(outcomes, dtype=float)

        if doses.size == 0:
            raise ValueError("doses must not be empty")
        if outcomes.size != doses.size:
            raise ValueError("doses and outcomes must have the same length")
        if not np.all(np.isin(outcomes, (0.0, 1.0))):
            raise ValueError("outcomes must be binary (0 or 1)")

        x0 = np.array([self.d50_init, self.k_init], dtype=float)
        result = minimize(
            self.log_likelihood,
            x0,
            args=(doses, outcomes),
            method="L-BFGS-B",
            bounds=self.bounds,
        )

        if not result.success:
            raise RuntimeError(f"Logistic TCP MLE failed: {result.message}")

        self.params_ = {"D50_gy": float(result.x[0]), "k": float(result.x[1])}
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
    """Fit Logistic TCP on modeling cohort EQD2 with median-OS binary proxy (sanity check)."""
    from sklearn.metrics import roc_auc_score

    from src.config import DATA_PROCESSED, RANDOM_SEED

    frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    doses = frame["eqd2_gy"].to_numpy()
    median_os = frame["survival_weeks"].median()
    outcomes = (frame["survival_weeks"] >= median_os).astype(float).to_numpy()

    model = LogisticTCPModel(d50_init=53.0, k_init=0.1)
    model.fit(doses, outcomes)
    preds = model.predict(doses)

    print(f"Cohort: n={len(frame)}, median OS split at {median_os:.0f} wk")
    print(model.summary().to_string(index=False))
    print(f"Negative log-likelihood: {model.nll_:.2f}")
    print(f"ROC AUC: {roc_auc_score(outcomes, preds):.4f}")
    print(f"Random seed (project default): {RANDOM_SEED}")


if __name__ == "__main__":
    main()
