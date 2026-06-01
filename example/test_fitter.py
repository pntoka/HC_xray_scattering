"""Tests for IOBSFitter using the bundled example data.

Run as a pytest suite:
    pytest example/test_fitter.py -v

Or directly as a script:
    python example/test_fitter.py
"""

import json
import os

import numpy as np
import pytest

EXAMPLE_DIR = os.path.dirname(os.path.abspath(__file__))
FIT_DATA_FILE = os.path.join(EXAMPLE_DIR, "example_fit.xy")
PARAMS_FILE = os.path.join(EXAMPLE_DIR, "example_params.json")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_data():
    """Load the example scattering data (s, intensity)."""
    data = np.loadtxt(FIT_DATA_FILE)
    return data[:, 0], data[:, 1]


def _load_params():
    """Load the example parameter dict (the 'true' params that generated the data)."""
    with open(PARAMS_FILE) as fh:
        return json.load(fh)


def _extract_initial(params, param_bounds):
    """Pull only the fittable keys out of the full params dict."""
    return {k: float(params[k]) for k in param_bounds if k in params}


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_fitter_result_structure():
    """FitResult and StageResult objects must have the expected types/shapes."""
    from iobs_ngc import IOBSFitter, DEFAULT_PARAM_BOUNDS, FitResult, StageResult

    s, y = _load_data()
    params = _load_params()
    initial = _extract_initial(params, DEFAULT_PARAM_BOUNDS)

    fitter = IOBSFitter(
        template_params=params,
        # Tight nfev cap – we only care about structure here, not quality
        ls_kwargs={
            "max_nfev": 250,
            "ftol": 1e-15,
            "xtol": 1e-15,
            "gtol": 1e-15,
            "diff_step": 1e-8,
        },
    )
    result = fitter.fit(s, y, initial)

    # Top-level types
    assert isinstance(result, FitResult)
    assert isinstance(result.parameters, dict)
    assert isinstance(result.parameters_scaled, dict)
    assert isinstance(result.r_squared, float)
    assert isinstance(result.rss, float)
    assert isinstance(result.nfev, int)
    assert isinstance(result.success, bool)
    assert isinstance(result.stage_results, list)
    assert len(result.stage_results) > 0, "Expected at least one stage result"

    # Stage result types
    for sr in result.stage_results:
        assert isinstance(sr, StageResult)
        assert isinstance(sr.name, str)
        assert isinstance(sr.parameters, dict)
        assert isinstance(sr.success, bool)
        assert isinstance(sr.rss, float)
        assert isinstance(sr.nfev, int)
        assert isinstance(sr.convergence_info, str)

    # Array shapes
    assert result.fitted_values is not None
    assert result.residuals is not None
    assert result.fitted_values.shape == s.shape
    assert result.residuals.shape == y.shape

    # Scaled params are in [0, 1]
    for k, v in result.parameters_scaled.items():
        assert 0.0 <= v <= 1.0, f"Scaled param {k}={v} is outside [0, 1]"

    # Physical params are within their declared bounds
    for k, v in result.parameters.items():
        lo, hi = DEFAULT_PARAM_BOUNDS[k]
        assert lo - 1e-8 <= v <= hi + 1e-8, (
            f"Physical param {k}={v:.6f} is outside bounds [{lo}, {hi}]"
        )

    # Timing fields are non-negative
    assert result.wall_time_s >= 0.0
    for sr in result.stage_results:
        assert sr.wall_time_s >= 0.0

    print(result.summary())


def test_fitter_known_params():
    """
    Starting from the parameters that generated the data should yield an
    essentially perfect fit (R² > 0.9999) with very few iterations.
    """
    from iobs_ngc import IOBSFitter, DEFAULT_PARAM_BOUNDS

    s, y = _load_data()
    params = _load_params()
    initial = _extract_initial(params, DEFAULT_PARAM_BOUNDS)

    fitter = IOBSFitter(
        template_params=params,
        ls_kwargs={
            "max_nfev": 200,
            "ftol": 1e-12,
            "xtol": 1e-12,
            "gtol": 1e-12,
            "diff_step": 1e-6,
        },
    )
    result = fitter.fit(s, y, initial)

    print(result.summary())

    assert result.r_squared > 0.9999, (
        f"Expected near-perfect fit from known params, got R²={result.r_squared:.6f}"
    )


def test_fitter_perturbed_params():
    """
    Starting from 5 % relatively-perturbed parameters should still converge
    to a good fit (R² > 0.99) using the multi-stage approach.
    """
    from iobs_ngc import IOBSFitter, DEFAULT_PARAM_BOUNDS

    s, y = _load_data()
    params = _load_params()

    rng = np.random.default_rng(42)
    initial = {}
    for k in DEFAULT_PARAM_BOUNDS:
        if k not in params:
            continue
        lo, hi = DEFAULT_PARAM_BOUNDS[k]
        val = float(params[k])
        # ±5 % relative noise, clipped strictly inside bounds
        noise = 0.05 * abs(val) * float(rng.uniform(-1.0, 1.0))
        initial[k] = float(np.clip(val + noise, lo + 1e-6, hi - 1e-6))

    fitter = IOBSFitter(template_params=params)
    result = fitter.fit(s, y, initial)

    print(result.summary())

    assert result.r_squared > 0.99, (
        f"Recovery from perturbed params failed: R²={result.r_squared:.6f}"
    )


# ---------------------------------------------------------------------------
# Standalone runner
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    tests = [
        test_fitter_result_structure,
        test_fitter_known_params,
        test_fitter_perturbed_params,
    ]
    failed = []
    for t in tests:
        print(f"\n{'=' * 55}")
        print(f"  {t.__name__}")
        print("=" * 55)
        try:
            t()
            print("  PASSED ✓")
        except AssertionError as exc:
            print(f"  FAILED ✗  {exc}")
            failed.append(t.__name__)
        except Exception as exc:
            print(f"  ERROR  ✗  {type(exc).__name__}: {exc}")
            failed.append(t.__name__)

    print(f"\n{'=' * 55}")
    if failed:
        print(f"  {len(failed)} test(s) FAILED: {', '.join(failed)}")
    else:
        print(f"  All {len(tests)} tests PASSED ✓")
