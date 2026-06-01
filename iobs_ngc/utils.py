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

# Fixed value of ``nu`` used when neither an explicit argument nor a template
# parameter is supplied (matches bundled default_params.json).
_DEFAULT_NU = 4.0


# ---------------------------------------------------------------------------
# Microstructure report
# ---------------------------------------------------------------------------

# Refined parameters: name -> (group, description).  Descriptions are taken
# from example/param_descriptors.md ("Relationship to microstructure
# parameters").  Order here defines the order shown in the report.
_REFINED_PARAM_INFO = {
    'a3':     ('microstructure', 'Average layer distance (Å)'),
    'da3':    ('microstructure', 'Avg − minimal layer distance (Å)'),
    'sig3':   ('microstructure', 'Disorder of the stacks (std dev of a3)'),
    'eta':    ('microstructure', 'Homogeneity of the stacks'),
    'mu':     ('microstructure', 'Shape factor for stack-height distribution'),
    'beta':   ('microstructure', 'Shape factor for stack-height distribution'),
    'alpha':  ('microstructure', 'Shape factor for layer-size distribution'),
    'sig1':   ('microstructure', 'Disorder of the layers (stress and strain)'),
    'lcc':    ('microstructure', 'Average C–C bond length (Å)'),
    'q':      ('microstructure', 'Preferred orientation'),
    'u3':     ('microstructure', 'Thermal motion'),
    'k':      ('normalization', 'Normalization (intensity scale) constant'),
    'const1': ('normalization', 'Constant background shift'),
    'const2': ('normalization', 'Linear background shift'),
    'g':      ('normalization', 'Exponential damping factor'),
}

# Calculated parameters: (name, formula, description, fn(refined, nu) -> float).
# ``fn`` may raise ZeroDivisionError; callers render those as "—".
_CALCULATED_PARAMS = [
    ('a3min', 'a3 − da3',        'Minimal layer distance (Å)',
     lambda p, nu: p['a3'] - p['da3']),
    ('N',     '(mu+1)/beta',     'Avg. graphene layers per stack',
     lambda p, nu: (p['mu'] + 1.0) / p['beta']),
    ('Nm',    'mu/beta',         'Avg. graphene layers per stack',
     lambda p, nu: p['mu'] / p['beta']),
    ('Lc',    'N·a3',            'Average stack height (Å)',
     lambda p, nu: (p['mu'] + 1.0) / p['beta'] * p['a3']),
    ('kapc',  '1/mu',            'Polydispersity of stack height',
     lambda p, nu: 1.0 / p['mu']),
    ('eps3',  'da3/a3min',       'Stack disorder from local strains',
     lambda p, nu: p['da3'] / (p['a3'] - p['da3'])),
    ('La',    '(nu+1)/alpha',    'Average graphene layer size (Å)',
     lambda p, nu: (nu + 1.0) / p['alpha']),
    ('lm',    'nu/alpha',        'Average chord length (Å)',
     lambda p, nu: nu / p['alpha']),
    ('kapa',  '1/nu',            'Polydispersity of chord length',
     lambda p, nu: 1.0 / nu),
    ('kapr',  '3π²(1/nu+1)/32−1', 'Polydispersity of layer radius',
     lambda p, nu: 3.0 * np.pi ** 2 * (1.0 / nu + 1.0) / 32.0 - 1.0),
]


def _fmt_value(value: float) -> str:
    """Format a numeric value for the report, falling back to a dash."""
    if value is None or not np.isfinite(value):
        return "—"
    return f"{value:.4f}"


def microstructure_report(
    result: FitResult,
    template_params: Optional[Dict[str, Any]] = None,
    nu: Optional[float] = None,
) -> str:
    """Build a human-readable microstructure report from a fit result.

    The report lists every refined parameter and the microstructure quantities
    calculated from them (each with a short description), preceded by a minimal
    fit-quality header (R², RSS, nfev).  It is the physical-interpretation
    counterpart to :meth:`FitResult.summary`, which reports timing/convergence.

    Parameters
    ----------
    result : FitResult
        Result returned by :class:`~iobs_ngc.IOBSFitter`.
    template_params : dict, optional
        Full parameter template.  Only used to source ``nu`` (a fixed, non-fitted
        parameter) for the calculated quantities.
    nu : float, optional
        Explicit value of ``nu``.  Resolution order: this argument →
        ``template_params['nu']`` → default (``4``).

    Returns
    -------
    str
        Formatted report.  The caller decides whether to print it.
    """
    if nu is None:
        if template_params is not None and 'nu' in template_params:
            nu = float(template_params['nu'])
        else:
            nu = _DEFAULT_NU

    params = result.parameters

    lines = [
        "─── Microstructure Report " + "─" * 31,
        f"  R²    : {result.r_squared:.8f}",
        f"  RSS   : {result.rss:.6e}",
        f"  nfev  : {result.nfev}",
    ]

    for group, heading in (
        ('microstructure', "Refined — microstructure"),
        ('normalization', "Refined — normalization / background"),
    ):
        rows = [
            (name, params[name], desc)
            for name, (grp, desc) in _REFINED_PARAM_INFO.items()
            if grp == group and name in params
        ]
        if not rows:
            continue
        lines.append("")
        lines.append(f"  {heading}")
        lines.append(f"    {'name':<8} {'value':>10}   description")
        for name, value, desc in rows:
            lines.append(f"    {name:<8} {_fmt_value(value):>10}   {desc}")

    lines.append("")
    lines.append("  Calculated — microstructure")
    lines.append(
        f"    {'name':<8} {'value':>10}   {'formula':<16} description"
    )
    for name, formula, desc, fn in _CALCULATED_PARAMS:
        try:
            value = float(fn(params, nu))
        except (ZeroDivisionError, KeyError, ValueError):
            value = float('nan')
        lines.append(
            f"    {name:<8} {_fmt_value(value):>10}   {formula:<16} {desc}"
        )

    lines.append("─" * 57)
    return "\n".join(lines)


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

    def report(self) -> str:
        """Return the microstructure report for the fit. Requires fit() first.

        See :func:`microstructure_report` for the report contents.
        """
        if self.result is None:
            raise RuntimeError("Call fit() before report().")
        return microstructure_report(self.result, self._template_params)

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
