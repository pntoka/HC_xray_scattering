"""Multi-stage curve fitter for iobs_ngc."""

from __future__ import annotations

import time
import warnings
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any

from scipy.optimize import least_squares

from .parameter_bounds import DEFAULT_PARAM_BOUNDS, DEFAULT_STAGE_DEFINITIONS
from .iobs_ngc import IOBSParameters as _IOBSParameters, IOBSCalculator as _IOBSCalculator


@dataclass
class StageResult:
    name: str
    parameters: Dict[str, float]
    success: bool
    rss: float
    nfev: int
    convergence_info: str
    retry_count: int = 0
    wall_time_s: float = 0.0          # total wall time for this stage


@dataclass
class FitResult:
    """Result returned by IOBSFitter."""
    parameters: Dict[str, float]          # final fitted parameters (physical units)
    parameters_scaled: Dict[str, float]   # final fitted parameters (0-1 scaled)
    r_squared: float
    rss: float
    nfev: int
    success: bool
    stage_results: List[StageResult] = field(default_factory=list)
    fitted_values: Optional[np.ndarray] = None
    residuals: Optional[np.ndarray] = None
    wall_time_s: float = 0.0              # total wall time for the full fit

    def summary(self) -> str:
        """Return a human-readable timing and quality summary string."""
        lines = [
            "─── IOBSFitter Result " + "─" * 35,
            f"  Status       : {'success' if self.success else 'FAILED'}",
            f"  R²           : {self.r_squared:.8f}",
            f"  RSS          : {self.rss:.6e}",
            f"  nfev (total) : {self.nfev}",
            f"  Wall time    : {self.wall_time_s:.3f} s",
            "",
            f"  {'Stage':<22} {'nfev':>6}  {'time(s)':>8}  {'rss':>12}",
            "  " + "─" * 55,
        ]
        for sr in self.stage_results:
            status = "✓" if sr.success else "✗"
            lines.append(
                f"  {status} {sr.name:<20} {sr.nfev:>6}  {sr.wall_time_s:>8.3f}  "
                f"{sr.rss:>12.4e}"
            )
        lines.append("─" * 57)
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.summary()


def _scale(value: float, lo: float, hi: float) -> float:
    return (value - lo) / (hi - lo)


def _unscale(scaled: float, lo: float, hi: float) -> float:
    return scaled * (hi - lo) + lo


class IOBSFitter:
    """
    Multi-stage least-squares fitter for the iobs_ngc calculator.

    Parameters
    ----------
    template_params : dict
        Full IOBSParameters dict used as a fixed template.  Only the
        parameters listed in ``param_bounds`` are varied during fitting.
    param_bounds : dict, optional
        {name: (min, max)} defining physical-space bounds and the
        min-max scaling applied internally.  Defaults to
        ``DEFAULT_PARAM_BOUNDS``.
    stage_definitions : list of dict, optional
        Sequential fitting stages.  Each entry must have keys
        ``'name'`` and ``'param_names'``.  Defaults to
        ``DEFAULT_STAGE_DEFINITIONS``.
    max_retries : int
        Number of times a stage is retried with perturbed p0 on failure.
    ls_kwargs : dict, optional
        Extra keyword arguments forwarded to ``scipy.optimize.least_squares``
        (e.g. ``max_nfev``, ``ftol``, ``xtol``, ``gtol``, ``diff_step``).
    """

    def __init__(
        self,
        template_params: Dict[str, Any],
        param_bounds: Optional[Dict[str, Tuple[float, float]]] = None,
        stage_definitions: Optional[List[Dict]] = None,
        max_retries: int = 3,
        ls_kwargs: Optional[Dict] = None,
    ):
        self.template_params = dict(template_params)
        self.param_bounds = param_bounds or DEFAULT_PARAM_BOUNDS
        self.stage_definitions = stage_definitions or DEFAULT_STAGE_DEFINITIONS
        self.max_retries = max_retries
        self.ls_kwargs = ls_kwargs or {
            'max_nfev': 500,
            'ftol': 1e-12,
            'xtol': 1e-12,
            'gtol': 1e-12,
            'diff_step': 1e-6,
        }
        # Cache a single stateless calculator instance (avoids per-evaluation overhead)
        self._calc = _IOBSCalculator()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def fit(
        self,
        s: np.ndarray,
        y: np.ndarray,
        initial_params: Dict[str, float],
    ) -> FitResult:
        """
        Fit the iobs model to experimental data.

        Parameters
        ----------
        s : array-like
            Scattering vector values.
        y : array-like
            Observed intensity values.
        initial_params : dict
            Starting values in **physical units** for the parameters
            listed in ``param_bounds``.

        Returns
        -------
        FitResult
        """
        s = np.asarray(s, dtype=float)
        y = np.asarray(y, dtype=float)

        # Build scaled [0, 1] starting point for every fittable parameter.
        # Priority: initial_params → template_params → midpoint (0.5).
        current_scaled: Dict[str, float] = {}
        for name, (lo, hi) in self.param_bounds.items():
            if name in initial_params:
                raw = float(initial_params[name])
            elif name in self.template_params:
                raw = float(self.template_params[name])
                warnings.warn(
                    f"Parameter '{name}' not in initial_params; "
                    f"using template value {raw}."
                )
            else:
                current_scaled[name] = 0.5
                warnings.warn(
                    f"Parameter '{name}' not found anywhere; using mid-range default."
                )
                continue
            current_scaled[name] = float(np.clip(_scale(raw, lo, hi), 0.0, 1.0))

        # ------------------------------------------------------------------
        # Pre-filter: identify s-values where the model is ill-conditioned
        # at the starting parameters and exclude them from the optimisation.
        # The element-level NaN guard in _residuals acts as a safety net for
        # any additional bad points encountered during the stage sweeps.
        # ------------------------------------------------------------------
        y_init = self._evaluate(s, current_scaled)
        valid_mask = np.isfinite(y_init)
        n_bad = int((~valid_mask).sum())
        if n_bad > 0:
            warnings.warn(
                f"{n_bad}/{len(s)} s-values produce NaN/Inf at the starting "
                f"parameters and will be excluded from optimisation."
            )
        s_fit = s[valid_mask]
        y_fit_data = y[valid_mask]

        stage_results: List[StageResult] = []
        total_nfev = 0
        t_fit_start = time.perf_counter()

        for stage_def in self.stage_definitions:
            stage_name = stage_def['name']
            # Filter to params that are actually in current_scaled
            names_to_fit = [n for n in stage_def['param_names'] if n in current_scaled]
            if not names_to_fit:
                warnings.warn(f"Stage '{stage_name}' has no valid parameters; skipping.")
                continue

            result, nfev = self._fit_stage(
                s_fit, y_fit_data, current_scaled, names_to_fit, stage_name
            )
            total_nfev += nfev

            if result.success:
                current_scaled.update(result.parameters)

            stage_results.append(result)

        # Build final physical-space parameters
        final_physical = {
            name: _unscale(val, *self.param_bounds[name])
            for name, val in current_scaled.items()
        }

        # Final evaluation on the full s range; compute quality metrics on
        # valid (finite) points only.
        y_full = self._evaluate(s, current_scaled)
        valid_final = np.isfinite(y_full)
        y_obs_valid = y[valid_final]
        y_pred_valid = y_full[valid_final]
        residuals_valid = y_obs_valid - y_pred_valid
        rss = float(np.sum(residuals_valid ** 2))
        ss_tot = float(np.sum((y_obs_valid - y_obs_valid.mean()) ** 2))
        r2 = 1.0 - rss / ss_tot if ss_tot > 0 else 0.0

        # Full-length arrays for plotting (NaN where model is ill-conditioned)
        residuals_full = np.full(len(s), np.nan)
        residuals_full[valid_final] = residuals_valid

        return FitResult(
            parameters=final_physical,
            parameters_scaled=dict(current_scaled),
            r_squared=r2,
            rss=rss,
            nfev=total_nfev,
            success=bool(stage_results) and all(sr.success for sr in stage_results),
            stage_results=stage_results,
            fitted_values=y_full,
            residuals=residuals_full,
            wall_time_s=time.perf_counter() - t_fit_start,
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _evaluate(self, s: np.ndarray, scaled_params: Dict[str, float]) -> np.ndarray:
        """Run the C++ calculator with the given scaled parameters.

        Returns the raw result array.  Individual NaN/Inf values indicate
        s-values where the model is numerically ill-conditioned; callers are
        responsible for handling them (e.g. by zeroing residuals at those
        points or pre-filtering the data).
        """
        params_dict = dict(self.template_params)
        for name, scaled_val in scaled_params.items():
            lo, hi = self.param_bounds[name]
            params_dict[name] = _unscale(scaled_val, lo, hi)

        iobs_params = _IOBSParameters(params_dict)
        return np.asarray(self._calc.calculate(iobs_params, s), dtype=float)

    def _residuals(
        self,
        stage_vals: np.ndarray,
        s: np.ndarray,
        y: np.ndarray,
        current_scaled: Dict[str, float],
        names_to_fit: List[str],
    ) -> np.ndarray:
        """Residual function passed to least_squares.

        NaN/Inf values returned by the model at specific s-points are replaced
        with zero residual so that those points are excluded from the
        optimisation rather than poisoning the entire cost function.
        """
        trial = dict(current_scaled)
        for name, val in zip(names_to_fit, stage_vals):
            trial[name] = float(val)

        try:
            y_pred = self._evaluate(s, trial)
            residuals = y_pred - y
            bad = ~np.isfinite(residuals)
            if bad.any():
                residuals[bad] = 0.0
            return residuals
        except Exception:
            return np.ones_like(y) * 1e6

    def _fit_stage(
        self,
        s: np.ndarray,
        y: np.ndarray,
        current_scaled: Dict[str, float],
        names_to_fit: List[str],
        stage_name: str,
    ) -> Tuple[StageResult, int]:
        """Run one fitting stage with retry logic."""
        p0 = np.array([current_scaled[n] for n in names_to_fit])
        # All scaled params live in [0, 1]
        lb = np.zeros(len(names_to_fit))
        ub = np.ones(len(names_to_fit))

        total_nfev = 0
        last_error = ""

        for attempt in range(self.max_retries):
            try:
                t_stage = time.perf_counter()
                res = least_squares(
                    fun=self._residuals,
                    x0=np.clip(p0, lb + 1e-9, ub - 1e-9),
                    bounds=(lb, ub),
                    method='trf',
                    args=(s, y, current_scaled, names_to_fit),
                    **self.ls_kwargs,
                )
                stage_wall = time.perf_counter() - t_stage
                total_nfev += res.nfev

                fitted = {
                    name: float(val)
                    for name, val in zip(names_to_fit, res.x)
                }
                rss = float(2.0 * res.cost)

                return StageResult(
                    name=stage_name,
                    parameters=fitted,
                    success=res.success,
                    rss=rss,
                    nfev=total_nfev,
                    convergence_info=res.message,
                    retry_count=attempt,
                    wall_time_s=stage_wall,
                ), total_nfev

            except Exception as exc:
                last_error = str(exc)
                warnings.warn(
                    f"Stage '{stage_name}' attempt {attempt + 1} failed: {exc}"
                )
                # Perturb p0 slightly for retry
                p0 = p0 * (1.0 + 0.02 * np.random.randn(len(p0)))
                p0 = np.clip(p0, lb + 1e-9, ub - 1e-9)

        # All retries failed
        return StageResult(
            name=stage_name,
            parameters={n: current_scaled[n] for n in names_to_fit},
            success=False,
            rss=np.inf,
            nfev=total_nfev,
            convergence_info=f"Failed after {self.max_retries} attempts: {last_error}",
            retry_count=self.max_retries,
        ), total_nfev