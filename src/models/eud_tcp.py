"""
gEUD-based TCP model with selectable volume-effect exponent ``a``.

Uses the Poisson TCP form applied to gEUD instead of uniform dose::

    TCP(gEUD) = exp( −ln(2) · exp( γ50 · (1 − gEUD/D50) ) )

When multiple gEUD exponents are available (a ∈ {−10, 1, 10}), ``fit_select_a``
picks the value minimising negative log-likelihood.
"""

from typing import Dict, Mapping, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from src.models.base_model import TCPModel
from src.models.poisson_tcp import LN2, poisson_tcp

_EPS = 1e-9


class EUDTCPModel(TCPModel):
    """
    gEUD-based TCP model fitted by maximum likelihood on binary outcomes.

    Parameters
    ----------
    d50_init : float, optional
        Initial guess for D50 (Gy).
    gamma50_init : float, optional
        Initial guess for γ50.
    geud_a : float, optional
        gEUD exponent used for the current fit.
    bounds : tuple of tuple, optional
        ``((d50_min, d50_max), (gamma50_min, gamma50_max))`` for MLE.
    """

    PARAM_NAMES: Tuple[str, ...] = ("D50_gy", "gamma50", "geud_a")

    def __init__(
        self,
        d50_init: float = 50.0,
        gamma50_init: float = 2.0,
        geud_a: float = 10.0,
        bounds: Optional[Tuple[Tuple[float, float], Tuple[float, float]]] = None,
    ) -> None:
        super().__init__()
        self.d50_init = d50_init
        self.gamma50_init = gamma50_init
        self.geud_a = geud_a
        self.bounds = bounds or ((1.0, 100.0), (0.01, 20.0))

    def predict(self, doses: np.ndarray) -> np.ndarray:
        """
        Predict TCP for given gEUD values using fitted parameters.

        Parameters
        ----------
        doses : np.ndarray
            gEUD values in Gy.

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
            gEUD values in Gy.
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

    def fit(
        self,
        doses: np.ndarray,
        outcomes: np.ndarray,
        geud_a: Optional[float] = None,
    ) -> "EUDTCPModel":
        """
        Fit D50 and γ50 by maximum likelihood for a fixed gEUD exponent.

        Parameters
        ----------
        doses : np.ndarray, shape (n_patients,)
            gEUD values in Gy for the chosen exponent ``a``.
        outcomes : np.ndarray, shape (n_patients,)
            Binary outcomes.
        geud_a : float, optional
            gEUD exponent; defaults to ``self.geud_a``.

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

        a_value = float(self.geud_a if geud_a is None else geud_a)
        x0 = np.array([self.d50_init, self.gamma50_init], dtype=float)
        result = minimize(
            self.log_likelihood,
            x0,
            args=(doses, outcomes),
            method="L-BFGS-B",
            bounds=self.bounds,
        )

        if not result.success:
            raise RuntimeError(f"EUD TCP MLE failed: {result.message}")

        self.params_ = {
            "D50_gy": float(result.x[0]),
            "gamma50": float(result.x[1]),
            "geud_a": a_value,
        }
        self.geud_a = a_value
        self.fitted_ = True
        self.nll_ = float(result.fun)
        self.n_obs_ = int(doses.size)
        return self

    @classmethod
    def fit_select_a(
        cls,
        geud_by_a: Mapping[float, np.ndarray],
        outcomes: np.ndarray,
        **kwargs,
    ) -> "EUDTCPModel":
        """
        Fit the model for each candidate ``a`` and return the best by NLL.

        Parameters
        ----------
        geud_by_a : mapping
            Maps gEUD exponent ``a`` to per-patient gEUD arrays (Gy).
        outcomes : np.ndarray
            Binary outcomes.
        **kwargs
            Forwarded to ``EUDTCPModel`` constructor.

        Returns
        -------
        EUDTCPModel
            Fitted model with lowest negative log-likelihood.

        Raises
        ------
        ValueError
            If ``geud_by_a`` is empty.
        """
        if not geud_by_a:
            raise ValueError("geud_by_a must not be empty")

        best: Optional[EUDTCPModel] = None
        for a_value, doses in geud_by_a.items():
            candidate = cls(geud_a=float(a_value), **kwargs)
            candidate.fit(np.asarray(doses, dtype=float), outcomes, geud_a=float(a_value))
            if best is None or candidate.nll_ < best.nll_:
                best = candidate
        assert best is not None
        return best

    @staticmethod
    def geud_columns_from_frame(frame: pd.DataFrame) -> Dict[float, np.ndarray]:
        """
        Extract precomputed gEUD columns from a modeling table.

        Parameters
        ----------
        frame : pd.DataFrame
            Must contain ``gEUD_am10_gy``, ``gEUD_a1_gy``, ``gEUD_a10_gy``.

        Returns
        -------
        dict
            Maps ``a`` to gEUD array.
        """
        mapping = {
            -10.0: "gEUD_am10_gy",
            1.0: "gEUD_a1_gy",
            10.0: "gEUD_a10_gy",
        }
        return {a: frame[col].to_numpy(dtype=float) for a, col in mapping.items() if col in frame.columns}

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
    """Fit EUD TCP selecting best a on modeling cohort (median-OS proxy)."""
    from sklearn.metrics import roc_auc_score

    from src.config import DATA_PROCESSED, RANDOM_SEED

    frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
    median_os = frame["survival_weeks"].median()
    outcomes = (frame["survival_weeks"] >= median_os).astype(float).to_numpy()

    geud_cols = EUDTCPModel.geud_columns_from_frame(frame)
    model = EUDTCPModel.fit_select_a(geud_cols, outcomes)
    doses = geud_cols[model.params_["geud_a"]]
    preds = model.predict(doses)

    print(f"Cohort: n={len(frame)}, median OS split at {median_os:.0f} wk")
    print(f"Candidate a values: {sorted(geud_cols)}")
    print(model.summary().to_string(index=False))
    print(f"Negative log-likelihood: {model.nll_:.2f}")
    print(f"ROC AUC: {roc_auc_score(outcomes, preds):.4f}")
    print(f"Random seed (project default): {RANDOM_SEED}")


if __name__ == "__main__":
    main()
