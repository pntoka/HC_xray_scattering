# BuildDependencies.cmake - Builds FLINT/ARB and dependencies from source

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

# Windows-specific: Use vcpkg for GMP/MPFR/ARB
if(WIN32)
    message(STATUS "Windows build: using vcpkg for GMP/MPFR/ARB")
    
    find_package(PkgConfig REQUIRED)
    
    pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)
    pkg_check_modules(gmpxx REQUIRED IMPORTED_TARGET gmpxx)
    pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
    
    # Find pthread
    find_library(PTHREAD_LIBRARY
        NAMES pthreadVC3
        PATHS 
            C:/vcpkg/installed/${VCPKG_TARGET_TRIPLET}/lib
            ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib
        NO_DEFAULT_PATH
        REQUIRED
    )
    
    # ARB often doesn't provide a .pc file, so find it manually
    find_library(ARB_LIBRARY
        NAMES arb
        PATHS 
            C:/vcpkg/installed/${VCPKG_TARGET_TRIPLET}/lib
            ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib
        NO_DEFAULT_PATH
        REQUIRED
    )
    
    find_path(ARB_INCLUDE_DIR
        NAMES arb.h
        PATHS 
            C:/vcpkg/installed/${VCPKG_TARGET_TRIPLET}/include
            ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include
        NO_DEFAULT_PATH
        REQUIRED
    )
    
    # Also need FLINT for arb
    find_library(FLINT_LIBRARY
        NAMES flint
        PATHS 
            C:/vcpkg/installed/${VCPKG_TARGET_TRIPLET}/lib
            ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib
        NO_DEFAULT_PATH
        REQUIRED
    )
    
    message(STATUS "Found ARB library: ${ARB_LIBRARY}")
    message(STATUS "Found ARB include: ${ARB_INCLUDE_DIR}")
    message(STATUS "Found FLINT library: ${FLINT_LIBRARY}")
    message(STATUS "Found pthread library: ${PTHREAD_LIBRARY}")
    
    # Create imported target for pthread
    add_library(pthread_imported STATIC IMPORTED)
    set_target_properties(pthread_imported PROPERTIES
        IMPORTED_LOCATION "${PTHREAD_LIBRARY}"
    )
    
    # Create imported target for arb
    add_library(arb_imported STATIC IMPORTED)
    set_target_properties(arb_imported PROPERTIES
        IMPORTED_LOCATION "${ARB_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${ARB_INCLUDE_DIR}"
    )
    
    # Create imported target for flint
    add_library(flint_imported STATIC IMPORTED)
    set_target_properties(flint_imported PROPERTIES
        IMPORTED_LOCATION "${FLINT_LIBRARY}"
    )
    
    # Bundle everything together
    add_library(flint_external INTERFACE)
    target_include_directories(flint_external INTERFACE ${ARB_INCLUDE_DIR})
    target_link_libraries(flint_external INTERFACE
        arb_imported
        flint_imported
        PkgConfig::mpfr
        PkgConfig::gmpxx
        PkgConfig::gmp
        pthread_imported
    )
    
    # Create dummy ep_flint target for compatibility
    add_custom_target(ep_flint)
    
    set(FLINT_TARGET flint_external)
    set(FLINT_INCLUDE_DIRS ${ARB_INCLUDE_DIR})
    return()
endif()

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
    BUILD_BYPRODUCTS 
        ${DEPS_LIB_DIR}/libgmp.a
        ${DEPS_LIB_DIR}/libgmpxx.a
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
    BUILD_BYPRODUCTS ${DEPS_LIB_DIR}/libmpfr.a
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

# ============================================================================
# FLINT (uses Autotools, depends on GMP and MPFR)
# ============================================================================
message(STATUS "Will build FLINT from source")

set(FLINT_ARCHIVE "${CMAKE_BINARY_DIR}/downloads/flint-3.1.2.tar.gz")
if(NOT EXISTS ${FLINT_ARCHIVE})
    message(STATUS "Downloading FLINT 3.1.2...")
    file(DOWNLOAD
        https://github.com/flintlib/flint/releases/download/v3.1.2/flint-3.1.2.tar.gz
        ${FLINT_ARCHIVE}
        EXPECTED_HASH SHA256=fdb3a431a37464834acff3bdc145f4fe8d0f951dd5327c4c6f93f4cbac5c2700
        SHOW_PROGRESS
        STATUS DOWNLOAD_STATUS
    )
    list(GET DOWNLOAD_STATUS 0 STATUS_CODE)
    list(GET DOWNLOAD_STATUS 1 ERROR_MESSAGE)
    if(NOT STATUS_CODE EQUAL 0)
        message(FATAL_ERROR "Failed to download FLINT: ${ERROR_MESSAGE}")
    endif()
    message(STATUS "FLINT download complete")
endif()

ExternalProject_Add(ep_flint
    URL ${FLINT_ARCHIVE}
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    CONFIGURE_COMMAND <SOURCE_DIR>/configure
        --prefix=${DEPS_INSTALL_DIR}
        --with-gmp=${DEPS_INSTALL_DIR}
        --with-mpfr=${DEPS_INSTALL_DIR}
        --disable-shared
        --enable-static
        --with-pic
        CFLAGS=-fPIC
        CXXFLAGS=-fPIC
    BUILD_COMMAND make -j${NPROC}
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE TRUE
    DEPENDS ep_mpfr
    BUILD_BYPRODUCTS ${DEPS_LIB_DIR}/libflint.a
    LOG_CONFIGURE OFF
    LOG_BUILD OFF
    LOG_INSTALL OFF
)

# Create imported targets with proper DEPENDS
add_library(gmp_imported STATIC IMPORTED GLOBAL)
set_target_properties(gmp_imported PROPERTIES
    IMPORTED_LOCATION ${DEPS_LIB_DIR}/libgmp.a
)
add_dependencies(gmp_imported ep_gmp)

add_library(gmpxx_imported STATIC IMPORTED GLOBAL)
set_target_properties(gmpxx_imported PROPERTIES
    IMPORTED_LOCATION ${DEPS_LIB_DIR}/libgmpxx.a
)
add_dependencies(gmpxx_imported ep_gmp)

add_library(mpfr_imported STATIC IMPORTED GLOBAL)
set_target_properties(mpfr_imported PROPERTIES
    IMPORTED_LOCATION ${DEPS_LIB_DIR}/libmpfr.a
)
add_dependencies(mpfr_imported ep_mpfr)

add_library(flint_imported STATIC IMPORTED GLOBAL)
set_target_properties(flint_imported PROPERTIES
    IMPORTED_LOCATION ${DEPS_LIB_DIR}/libflint.a
)
add_dependencies(flint_imported ep_flint)

# Create an interface target that bundles everything
add_library(flint_external INTERFACE)
target_include_directories(flint_external INTERFACE ${DEPS_INCLUDE_DIR})
target_link_libraries(flint_external INTERFACE
    flint_imported
    mpfr_imported
    gmpxx_imported
    gmp_imported
)

set(FLINT_TARGET flint_external)
set(FLINT_INCLUDE_DIRS ${DEPS_INCLUDE_DIR})