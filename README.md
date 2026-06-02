# HC_xray_scattering
Developing toolkit for refinement of non-graphitic carbon scattering data.

The refinement approach used in this repository is based on the method developed by W. Ruland and B. Smarsly ( J. Appl. Cryst., 35, 624-633, [doi:10.1107/S0021889802011007](https://doi.org/10.1107/S0021889802011007)) and the software implementation developed by O. Osswald (C - Journal of Carbon Research, 8(4), 78, [doi:10.3390/c8040078](https://doi.org/10.3390/c8040078))

## Overview

`iobs_ngc` models the X-ray scattering pattern of non-graphitic carbons following the
Ruland–Smarsly turbostratic disorder framework. The library provides:

- **`IOBSCalculator`** — low-level C++ calculator that evaluates the scattering model
  at a given set of s-values for a full parameter set.
- **`IOBSFitter`** — multi-stage least-squares fitter built on `scipy.optimize.least_squares`
  that refines 15 structural parameters sequentially before a final global refinement.
- **`FitPattern`** — high-level wrapper that handles data loading, parameter defaults, and
  result visualisation in a few lines of code.

### Parameters

The 15 fittable structural parameters and their physical meanings:

| Parameter | Description |
|-----------|-------------|
| `mu` | Mean number of layers |
| `beta` | Layer stacking disorder |
| `a3` | Mean interlayer spacing (Å) |
| `da3` | Interlayer spacing distribution width |
| `sig3` | Interlayer spacing standard deviation |
| `eta` | Lorentz/Gauss mixing for interlayer peak |
| `alpha` | In-plane crystallite size parameter |
| `sig1` | In-plane disorder parameter |
| `lcc` | In-plane C–C bond length (Å) |
| `q` | Probability of graphene-like stacking |
| `g` | Preferred orientation parameter |
| `u3` | Additional layer displacement parameter |
| `const1` | Background offset |
| `const2` | Background slope |
| `k` | Intensity scale factor |

For the full list of parameters and their detailed descriptions, see the
upstream [osswaldo/OctCarb](https://github.com/osswaldo/OctCarb) project.

#### Instrument-dependent parameters

Three additional parameters describe the measurement geometry rather than the
sample. They are **fixed during fitting (not refined)** and **must be supplied
by the user** to match their instrument — the bundled defaults are only
placeholders:

| Parameter | Description |
|-----------|-------------|
| `par_r` | Goniometer radius (cm; fixed by the experiment) |
| `par_delta` | Divergence angle (°; set by the user) |
| `par_l` | Irradiated sample length (cm; fixed during measurement) |

Provide them inline via the `params` dict argument:

```python
session = FitPattern("sample.xy", params={"par_r": 14, "par_delta": 4, "par_l": 7})
```

or in a `params_json` file. Any keys not given fall back to the bundled
defaults.

---

## Quick start

### Fitting an experimental pattern

The simplest workflow uses `FitPattern`. Your data file should be a two-column ASCII
`.xy` file with 2θ (degrees) in the first column and intensity in the second.

```python
from iobs_ngc import FitPattern

# Minimal — all parameters initialised to mid-range defaults
session = FitPattern("path/to/sample.xy")
session.fit()

# Inspect fitted parameters
print(session.fitted_params)

# Plot data, fit, and residuals
session.plot()                               # display interactively
session.plot(save_path="fit.png", show=False)  # save to file
```

If your data file already stores the scattering vector `s` (Å⁻¹) in the first
column instead of 2θ, pass `x_is_s=True` to skip the 2θ→s conversion:

```python
session = FitPattern("sample.xy", x_is_s=True)
```

Providing custom parameters overrides any matching keys; the rest are filled from
the bundled defaults. You can pass them inline as a dict or point to a JSON file:

```python
session = FitPattern("sample.xy", params={"par_r": 14, "par_delta": 4, "par_l": 7})
session = FitPattern("sample.xy", params_json="my_params.json")
```

### Getting a microstructure report

After fitting, `report()` returns a human-readable summary: a fit-quality
header (R², RSS, nfev), the refined structural and background parameters, and
the microstructure quantities calculated from them (e.g. crystallite sizes).

```python
session = FitPattern("sample.xy")
session.fit()
print(session.report())
```

For the lower-level workflow, build the report directly from a `FitResult`:

```python
from iobs_ngc import IOBSFitter, microstructure_report

result = IOBSFitter(template_params=template_params).fit(s, y, initial_params)
print(microstructure_report(result, template_params))
```

### Smoothing the data

Set `smooth=True` to apply a Savitzky-Golay filter to the (interpolated)
intensity before fitting — useful for suppressing high-frequency noise in
experimental patterns. The window length and polynomial order are configurable:

```python
session = FitPattern(
    "sample.xy",
    smooth=True,
    savgol_window=11,     # filter window length (odd, > polyorder)
    savgol_polyorder=3,   # polynomial order
)
session.fit()
print(session.report())
```

The same options are available on the lower-level `interp_spectra` helper
(`smooth`, `window_length`, `polyorder`).

Scipy solver tolerances can be tightened or relaxed via `ls_kwargs`:

```python
session = FitPattern(
    "sample.xy",
    ls_kwargs={
        "max_nfev": 500,
        "ftol": 1e-12,
        "xtol": 1e-12,
        "gtol": 1e-12,
        "diff_step": 1e-6,
    },
)
session.fit()
```

### Lower-level access

`IOBSFitter` and the data-loading helpers are available directly for custom workflows:

```python
import json
from iobs_ngc import IOBSFitter, DEFAULT_PARAM_BOUNDS
from iobs_ngc.utils import interp_spectra

# Load and convert data (2theta -> s, interpolated to 1000 points)
s, y = interp_spectra("data/", "sample.xy")

# Build initial params and template
with open("params.json") as fh:
    template_params = json.load(fh)
initial_params = {k: (lo + hi) / 2 for k, (lo, hi) in DEFAULT_PARAM_BOUNDS.items()}

fitter = IOBSFitter(template_params=template_params)
result = fitter.fit(s, y, initial_params)

print(result.summary())
print(result.parameters)   # final fitted values
```

---

## Installation

### Prerequisites

First, clone this repository:

```bash
git clone https://github.com/pntoka/HC_xray_scattering.git
cd HC_xray_scattering
```

All installation commands below should be run from within the cloned repository directory.

### Windows

#### Recommended: Install from Prebuilt Wheels

The easiest way to install on Windows is to use the prebuilt wheel files included in the repository:

```bash
pip install --find-links dist/ iobs_ngc
```

This will install `iobs_ngc` using the precompiled wheels from the `dist/` folder, avoiding the need to build from source.

#### Build from Source (Advanced)

If you need to build from source on Windows, use [vcpkg](https://github.com/microsoft/vcpkg) to install the required dependencies. This method is not recommended as may require additionally downloading Visual Studio and carrying out installation in the VS Developer Shell:

1. Install vcpkg:
   ```bash
   git clone https://github.com/Microsoft/vcpkg.git C:\vcpkg
   cd C:\vcpkg
   .\bootstrap-vcpkg.bat
   ```

2. Install dependencies (dynamic libraries):
   ```bash
   .\vcpkg.exe install gmp:x64-windows mpfr:x64-windows flint:x64-windows
   ```

3. Set environment variables and install:
   ```bash
   set VCPKG_ROOT=C:\vcpkg
   set CMAKE_TOOLCHAIN_FILE=C:\vcpkg\scripts\buildsystems\vcpkg.cmake
   pip install . --config-settings=cmake.define.USE_SYSTEM_LIBS=ON
   ```

### macOS

#### Recommended: Install from Prebuilt Wheels

Install using the prebuilt wheel file in the repository:

```bash
pip install --find-links dist/ iobs_ngc
```

#### Build from Source with System Libraries

To build from source using existing system libraries, first install the required dependencies via Homebrew:

```bash
brew install gmp mpfr flint boost libomp
```

Then install the package:

```bash
pip install . --config-settings=cmake.define.USE_SYSTEM_LIBS=ON --config-settings=cmake.define.FETCH_BOOST=OFF
```

### Ubuntu/Debian

For Ubuntu and Debian systems, you must first install GMP and MPFR system libraries before installing `iobs_ngc`. This is required because [FLINT version 3](https://flintlib.org/doc/) needs to be built from source during installation.

1. Install system dependencies:
   ```bash
   sudo apt-get update
   sudo apt-get install -y libgmp-dev libmpfr-dev build-essential
   ```

2. Install `iobs_ngc`:
   ```bash
   pip install .
   ```

The installation process will automatically build FLINT 3 from source using the system GMP and MPFR libraries.

If FLINT 3 is already installed then `iobs_ngc` can be installed using the command:
   ```bash
   pip install . --config-settings=cmake.define.USE_SYSTEM_LIBS=ON
   ```

### General Notes

- All installation methods require a C++17 compatible compiler
- Python 3.8 or higher is required
- For development installations, use `pip install -e .` instead of `pip install .` 
