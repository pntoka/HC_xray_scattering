# FindFLINT.cmake
find_path(FLINT_INCLUDE_DIR
    NAMES flint/flint.h
    PATHS /usr/include /usr/local/include
)

find_library(FLINT_LIBRARY
    NAMES flint
    PATHS /usr/lib /usr/local/lib /usr/lib/x86_64-linux-gnu
)

# FLINT requires GMP and MPFR
find_library(GMP_LIBRARY NAMES gmp)
find_library(MPFR_LIBRARY NAMES mpfr)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLINT
    REQUIRED_VARS FLINT_LIBRARY FLINT_INCLUDE_DIR GMP_LIBRARY MPFR_LIBRARY
)

if(FLINT_FOUND)
    set(FLINT_LIBRARIES ${FLINT_LIBRARY} ${GMP_LIBRARY} ${MPFR_LIBRARY})
    set(FLINT_INCLUDE_DIRS ${FLINT_INCLUDE_DIR})
endif()

mark_as_advanced(FLINT_INCLUDE_DIR FLINT_LIBRARY GMP_LIBRARY MPFR_LIBRARY)
