# BuildDependencies.cmake - Builds FLINT and dependencies from source

include(ExternalProject)

set(DEPS_INSTALL_DIR ${CMAKE_BINARY_DIR}/deps)
set(DEPS_INCLUDE_DIR ${DEPS_INSTALL_DIR}/include)
set(DEPS_LIB_DIR ${DEPS_INSTALL_DIR}/lib)

# Determine parallel jobs
include(ProcessorCount)
ProcessorCount(NPROC)
if(NPROC EQUAL 0)
    set(NPROC 1)
endif()
message(STATUS "Building dependencies with ${NPROC} parallel jobs")

# ============================================================================
# GMP (required by FLINT)
# ============================================================================
message(STATUS "Will build GMP from source")
ExternalProject_Add(ep_gmp
    URL https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
    URL_HASH SHA256=a3c2b80201b89e68616f4ad30bc66aee4927c3ce50e33929ca819d5c43538898
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/downloads
    CONFIGURE_COMMAND <SOURCE_DIR>/configure
        --prefix=${DEPS_INSTALL_DIR}
        --enable-cxx
        --with-pic
        --disable-shared
        --enable-static
    BUILD_COMMAND make -j${NPROC}
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE FALSE
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

# ============================================================================
# MPFR (required by FLINT, depends on GMP)
# ============================================================================
message(STATUS "Will build MPFR from source")
ExternalProject_Add(ep_mpfr
    URL https://ftp.gnu.org/gnu/mpfr/mpfr-4.2.1.tar.xz
    URL_HASH SHA256=277807353a6726978996945af13e52829e3abd7a9a5b7fb2793894e18f1fcbb2
    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/downloads
    CONFIGURE_COMMAND <SOURCE_DIR>/configure
        --prefix=${DEPS_INSTALL_DIR}
        --with-gmp=${DEPS_INSTALL_DIR}
        --with-pic
        --disable-shared
        --enable-static
    BUILD_COMMAND make -j${NPROC}
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE FALSE
    DEPENDS ep_gmp
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

# ============================================================================
# FLINT (uses CMake, depends on GMP and MPFR)
# ============================================================================
message(STATUS "Will build FLINT from source")
ExternalProject_Add(ep_flint
    GIT_REPOSITORY https://github.com/flintlib/flint.git
    GIT_TAG v3.1.2
    GIT_SHALLOW TRUE
    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/downloads
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${DEPS_INSTALL_DIR}
        -DCMAKE_PREFIX_PATH=${DEPS_INSTALL_DIR}
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DCMAKE_BUILD_TYPE=Release
        -DBUILD_SHARED_LIBS=OFF
        -DBUILD_TESTING=OFF
        -DBUILD_DOCS=OFF
    BUILD_COMMAND ${CMAKE_COMMAND} --build . --parallel ${NPROC}
    INSTALL_COMMAND ${CMAKE_COMMAND} --install .
    DEPENDS ep_mpfr
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

# Create an interface target that depends on the external projects
add_library(flint_external INTERFACE)
add_dependencies(flint_external ep_flint)

target_include_directories(flint_external INTERFACE
    ${DEPS_INCLUDE_DIR}
)

# Link order matters: flint needs mpfr, mpfr needs gmp
target_link_libraries(flint_external INTERFACE
    ${DEPS_LIB_DIR}/libflint.a
    ${DEPS_LIB_DIR}/libmpfr.a
    ${DEPS_LIB_DIR}/libgmp.a
    ${DEPS_LIB_DIR}/libgmpxx.a  # Add C++ bindings if available
)

# Export for parent scope
set(FLINT_TARGET flint_external PARENT_SCOPE)
set(FLINT_INCLUDE_DIRS ${DEPS_INCLUDE_DIR} PARENT_SCOPE)