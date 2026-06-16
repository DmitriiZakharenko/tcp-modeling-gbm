"""
Abstract base class for all TCP dose-response models.

Every model must implement fit, predict, and log_likelihood.
"""

from abc import ABC, abstractmethod
from typing import Dict

import numpy as np
import pandas as pd


class TCPModel(ABC):
    """
    Abstract base class for Tumor Control Probability dose-response models.

    Subclasses must implement:
        fit(doses, outcomes)
        predict(doses)
        log_likelihood(params, doses, outcomes)
        summary()
    """

    def __init__(self) -> None:
        self.params_: Dict[str, float] = {}
        self.fitted_: bool = False

    @abstractmethod
    def fit(self, doses: np.ndarray, outcomes: np.ndarray) -> "TCPModel":
        """
        Fit model parameters by maximum likelihood estimation.

        Parameters
        ----------
        doses : np.ndarray, shape (n_patients,)
            Dose values in Gy (EQD2-corrected).
        outcomes : np.ndarray, shape (n_patients,)
            Binary outcomes (1 = tumor controlled, 0 = failure).

        Returns
        -------
        self
        """

    @abstractmethod
    def predict(self, doses: np.ndarray) -> np.ndarray:
        """
        Predict TCP for a given array of doses.

        Parameters
        ----------
        doses : np.ndarray, shape (n,)
            Dose values in Gy.

        Returns
        -------
        np.ndarray, shape (n,)
            Predicted TCP values in [0, 1].
        """

    @abstractmethod
    def log_likelihood(self, params: np.ndarray, doses: np.ndarray, outcomes: np.ndarray) -> float:
        """
        Compute negative log-likelihood for use with scipy.optimize.minimize.

        Parameters
        ----------
        params : np.ndarray
            Model parameters in the order expected by the subclass.
        doses : np.ndarray, shape (n_patients,)
            Dose values in Gy.
        outcomes : np.ndarray, shape (n_patients,)
            Binary outcomes.

        Returns
        -------
        float
            Negative log-likelihood (minimization target).
        """

    @abstractmethod
    def summary(self) -> pd.DataFrame:
        """
        Return a DataFrame summarising fitted parameters.

        Returns
        -------
        pd.DataFrame
            Columns: parameter, estimate. One row per parameter.

        Raises
        ------
        RuntimeError
            If called before fit().
        """

    def _check_fitted(self) -> None:
        """Raise RuntimeError if model has not been fitted."""
        if not self.fitted_:
            raise RuntimeError(f"{self.__class__.__name__} must be fitted before calling this method.")
