# HC_xray_scattering
Developing toolkit for refinement of non-graphitic carbon scattering data

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
pip install . --find-links dist/
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
   .\vcpkg.exe install gmp:x64-windows mpfr:x64-windows flint:x64-windows arb:x64-windows
   ```

3. Set environment variables and install:
   ```bash
   set VCPKG_ROOT=C:\vcpkg
   set CMAKE_TOOLCHAIN_FILE=C:\vcpkg\scripts\buildsystems\vcpkg.cmake
   pip install . --config-settings=cmake.define.USE_SYSTEM_LIBS=ON
   ```

### macOS

#### Recommended: Install from Prebuilt Wheels

Install using the published wheel file:

```bash
pip install iobs_ngc
```

Or if installing from the repository:

```bash
pip install . --find-links dist/
```

#### Build from Source with System Libraries

To build from source using existing system libraries, first install the required dependencies via Homebrew:

```bash
brew install gmp mpfr flint boost
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
