"""iobs_ngc – Python package for non-graphitic carbon scattering calculations."""

# Re-export the compiled extension symbols
from .iobs_ngc import (   # noqa: F401
    IOBSParameters,
    IOBSCalculator,
    Calculations,
    RadiationType,
    CSP,
)

from .fitter import IOBSFitter, FitResult, StageResult   # noqa: F401
from .parameter_bounds import (                          # noqa: F401
    DEFAULT_PARAM_BOUNDS,
    DEFAULT_STAGE_DEFINITIONS,
)

from .utils import FitPattern                            # noqa: F401