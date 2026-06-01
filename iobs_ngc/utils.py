"""User-facing utilities for XRD pattern fitting with iobs_ngc."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, Optional, Any

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

from .fitter import IOBSFitter, FitResult
from .parameter_bounds import DEFAULT_PARAM_BOUNDS

_DEFAULT_PARAMS_FILE = Path(__file__).parent / "default_params.json"


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def open_xy(path: str, file_name: str) -> np.ndarray:
    with open(os.path.join(path, file_name), 'r') as f:
        data = f.readlines()
        data_clean = [line.replace('\t', ' ').split(' ') for line in data]
        data_array = np.array(data_clean).astype(float)
        sorted_data = data_array[data_array[:, 0].argsort()]
    return sorted_data


def twoThetaToS(x: np.ndarray, wavelength: float = 1.5418) -> np.ndarray:
    return 2 / wavelength * np.sin(x / 2 * np.pi / 180)


def interp_spectra(path: str, file_name: str):
    data = open_xy(path, file_name)
    x = twoThetaToS(data[:, 0])
    test_s = np.linspace(x.min(), x.max(), 1000)
    f_interp = interp1d(x, data[:, 1])
    y_interp = f_interp(test_s)
    return test_s, y_interp


# ---------------------------------------------------------------------------
# FitPattern
# ---------------------------------------------------------------------------

class FitPattern:
    """Load an XRD pattern, fit it, and inspect or plot the result.

    Parameters
    ----------
    xy_file : str or Path
        Path to the 2theta-format ``.xy`` data file.
    params_json : str or Path, optional
        Path to a JSON file with IOBSParameters values.  Any keys present
        override the package defaults; missing keys are filled from the
        bundled ``default_params.json``.
    ls_kwargs : dict, optional
        Extra keyword arguments forwarded to ``scipy.optimize.least_squares``
        (``max_nfev``, ``ftol``, ``xtol``, ``gtol``, ``diff_step``).
    """

    def __init__(
        self,
        xy_file: str | Path,
        params_json: Optional[str | Path] = None,
        ls_kwargs: Optional[Dict[str, Any]] = None,
    ):
        self._xy_file = Path(xy_file)
        self._ls_kwargs = ls_kwargs

        self._template_params, self._initial_params = self._build_params(
            params_json
        )

        self.result: Optional[FitResult] = None
        self._s: Optional[np.ndarray] = None
        self._y: Optional[np.ndarray] = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def fit(self) -> FitResult:
        """Load the data and run the multistage fit. Returns the FitResult."""
        self._s, self._y = interp_spectra(
            str(self._xy_file.parent), self._xy_file.name
        )
        fitter = IOBSFitter(
            template_params=self._template_params,
            ls_kwargs=self._ls_kwargs,
        )
        self.result = fitter.fit(self._s, self._y, self._initial_params)
        return self.result

    @property
    def fitted_params(self) -> Dict[str, float]:
        """Final fitted parameters in physical units. Requires fit() first."""
        if self.result is None:
            raise RuntimeError("Call fit() before accessing fitted_params.")
        return self.result.parameters

    def plot(
        self,
        save_path: Optional[str | Path] = None,
        show: bool = True,
    ) -> plt.Figure:
        """Plot the experimental data, fit, and residuals.

        Parameters
        ----------
        save_path : str or Path, optional
            If given, save the figure to this path (150 dpi).
        show : bool
            If True, call plt.show().
        """
        if self.result is None:
            raise RuntimeError("Call fit() before plotting.")

        fig, (ax_fit, ax_res) = plt.subplots(
            2, 1, figsize=(8, 7), sharex=True,
            gridspec_kw={"height_ratios": [3, 1]},
        )

        ax_fit.plot(self._s, self._y, lw=1, label="Experiment",
                    color="steelblue")
        if self.result.fitted_values is not None:
            ax_fit.plot(
                self._s, self.result.fitted_values, lw=1.5,
                label="Fit", color="tomato", linestyle="--",
            )
        ax_fit.set_ylabel("Intensity (a.u.)")
        ax_fit.legend()
        status = "converged" if self.result.success else "not converged"
        ax_fit.set_title(
            f"{self._xy_file.name}  —  "
            f"R² = {self.result.r_squared:.4f}  ({status})"
        )

        if self.result.residuals is not None:
            ax_res.plot(self._s, self.result.residuals, lw=0.8, color="gray")
        ax_res.axhline(0, color="black", lw=0.8, linestyle="--")
        ax_res.set_xlabel("s (Å⁻¹)")
        ax_res.set_ylabel("Residual")

        fig.tight_layout()

        if save_path is not None:
            fig.savefig(save_path, dpi=150)
        if show:
            plt.show()
            plt.close(fig)
            return None

        return fig

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_params(self, params_json):
        with open(_DEFAULT_PARAMS_FILE) as fh:
            template = json.load(fh)

        if params_json is not None:
            with open(params_json) as fh:
                user = json.load(fh)
            template.update(user)

        initial = {
            k: (lo + hi) / 2.0
            for k, (lo, hi) in DEFAULT_PARAM_BOUNDS.items()
        }
        return template, initial
