# FindFLINT.cmake - Cross-platform finder supporting conda environments

# Build search paths list
set(_FLINT_SEARCH_PATHS "")

# 1. Check for conda environment (highest priority)
if(DEFINED ENV{CONDA_PREFIX})
    list(APPEND _FLINT_SEARCH_PATHS "$ENV{CONDA_PREFIX}")
    message(STATUS "FindFLINT: Using CONDA_PREFIX=$ENV{CONDA_PREFIX}")
endif()

# 2. Check for conda-build environment (used during conda build)
if(DEFINED ENV{PREFIX})
    list(APPEND _FLINT_SEARCH_PATHS "$ENV{PREFIX}")
    message(STATUS "FindFLINT: Using PREFIX=$ENV{PREFIX}")
endif()

# 3. Windows-specific conda-build paths
if(WIN32 AND DEFINED ENV{LIBRARY_PREFIX})
    list(APPEND _FLINT_SEARCH_PATHS "$ENV{LIBRARY_PREFIX}")
endif()

# 4. Platform-specific system paths (fallback)
if(WIN32)
    list(APPEND _FLINT_SEARCH_PATHS
        "C:/msys64/mingw64"
        "C:/vcpkg/installed/x64-windows"
    )
elseif(APPLE)
    list(APPEND _FLINT_SEARCH_PATHS
        "/opt/homebrew"
        "/usr/local"
    )
else()
    # Linux
    list(APPEND _FLINT_SEARCH_PATHS
        "/usr"
        "/usr/local"
    )
endif()

# Find FLINT header
find_path(FLINT_INCLUDE_DIR
    NAMES flint/flint.h
    PATHS ${_FLINT_SEARCH_PATHS}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
)

# Fallback to default paths if not found
if(NOT FLINT_INCLUDE_DIR)
    find_path(FLINT_INCLUDE_DIR NAMES flint/flint.h)
endif()

# Find FLINT library
find_library(FLINT_LIBRARY
    NAMES flint
    PATHS ${_FLINT_SEARCH_PATHS}
    PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu
    NO_DEFAULT_PATH
)

if(NOT FLINT_LIBRARY)
    find_library(FLINT_LIBRARY NAMES flint)
endif()

# Find GMP library
find_library(GMP_LIBRARY
    NAMES gmp
    PATHS ${_FLINT_SEARCH_PATHS}
    PATH_SUFFIXES lib lib64
    NO_DEFAULT_PATH
)

if(NOT GMP_LIBRARY)
    find_library(GMP_LIBRARY NAMES gmp)
endif()

# Find MPFR library
find_library(MPFR_LIBRARY
    NAMES mpfr
    PATHS ${_FLINT_SEARCH_PATHS}
    PATH_SUFFIXES lib lib64
    NO_DEFAULT_PATH
)

if(NOT MPFR_LIBRARY)
    find_library(MPFR_LIBRARY NAMES mpfr)
endif()

# Standard CMake package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLINT
    REQUIRED_VARS FLINT_LIBRARY FLINT_INCLUDE_DIR GMP_LIBRARY MPFR_LIBRARY
)

if(FLINT_FOUND)
    set(FLINT_LIBRARIES ${FLINT_LIBRARY} ${GMP_LIBRARY} ${MPFR_LIBRARY})
    set(FLINT_INCLUDE_DIRS ${FLINT_INCLUDE_DIR})
    
    message(STATUS "FLINT found:")
    message(STATUS "  Include: ${FLINT_INCLUDE_DIR}")
    message(STATUS "  Libraries: ${FLINT_LIBRARIES}")
endif()

mark_as_advanced(FLINT_INCLUDE_DIR FLINT_LIBRARY GMP_LIBRARY MPFR_LIBRARY)
