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
    
    # Find DLLs for wheel repair
    find_file(ARB_DLL
        NAMES arb.dll
        PATHS 
            C:/vcpkg/installed/${VCPKG_TARGET_TRIPLET}/bin
            ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin
        NO_DEFAULT_PATH
    )
    find_file(FLINT_DLL
        NAMES flint.dll
        PATHS 
            C:/vcpkg/installed/${VCPKG_TARGET_TRIPLET}/bin
            ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin
        NO_DEFAULT_PATH
    )
    
    message(STATUS "Found ARB library: ${ARB_LIBRARY}")
    message(STATUS "Found ARB include: ${ARB_INCLUDE_DIR}")
    message(STATUS "Found FLINT library: ${FLINT_LIBRARY}")
    message(STATUS "Found pthread library: ${PTHREAD_LIBRARY}")
    if(ARB_DLL)
        message(STATUS "Found ARB DLL: ${ARB_DLL}")
    endif()
    if(FLINT_DLL)
        message(STATUS "Found FLINT DLL: ${FLINT_DLL}")
    endif()
    
    # Create imported target for pthread (could be static or dynamic)
    add_library(pthread_imported SHARED IMPORTED)
    set_target_properties(pthread_imported PROPERTIES
        IMPORTED_IMPLIB "${PTHREAD_LIBRARY}"
    )
    
    # Create imported target for arb
    add_library(arb_imported SHARED IMPORTED)
    set_target_properties(arb_imported PROPERTIES
        IMPORTED_IMPLIB "${ARB_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${ARB_INCLUDE_DIR}"
    )
    
    # Create imported target for flint
    add_library(flint_imported SHARED IMPORTED)
    set_target_properties(flint_imported PROPERTIES
        IMPORTED_IMPLIB "${FLINT_LIBRARY}"
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
# Linux/macOS: Use system GMP/MPFR, build FLINT from source
# ============================================================================
message(STATUS "Linux/macOS build: using system GMP/MPFR, building FLINT from source")

# Find system GMP and MPFR
find_library(GMP_LIBRARY NAMES gmp REQUIRED)
find_library(GMPXX_LIBRARY NAMES gmpxx REQUIRED)
find_library(MPFR_LIBRARY NAMES mpfr REQUIRED)
find_path(GMP_INCLUDE_DIR NAMES gmp.h REQUIRED)
find_path(MPFR_INCLUDE_DIR NAMES mpfr.h REQUIRED)

message(STATUS "Found system GMP: ${GMP_LIBRARY}")
message(STATUS "Found system GMP include: ${GMP_INCLUDE_DIR}")
message(STATUS "Found system MPFR: ${MPFR_LIBRARY}")
message(STATUS "Found system MPFR include: ${MPFR_INCLUDE_DIR}")

# ============================================================================
# FLINT (uses Autotools, depends on system GMP and MPFR)
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
        --disable-shared
        --enable-static
        --with-pic
        CFLAGS=-fPIC
        CXXFLAGS=-fPIC
    BUILD_COMMAND make -j${NPROC}
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE TRUE
    BUILD_BYPRODUCTS ${DEPS_LIB_DIR}/libflint.a
    LOG_CONFIGURE OFF
    LOG_BUILD OFF
    LOG_INSTALL OFF
)

# Create imported targets for system libraries
add_library(gmp_imported SHARED IMPORTED GLOBAL)
set_target_properties(gmp_imported PROPERTIES
    IMPORTED_LOCATION ${GMP_LIBRARY}
    INTERFACE_INCLUDE_DIRECTORIES ${GMP_INCLUDE_DIR}
)

add_library(gmpxx_imported SHARED IMPORTED GLOBAL)
set_target_properties(gmpxx_imported PROPERTIES
    IMPORTED_LOCATION ${GMPXX_LIBRARY}
    INTERFACE_INCLUDE_DIRECTORIES ${GMP_INCLUDE_DIR}
)

add_library(mpfr_imported SHARED IMPORTED GLOBAL)
set_target_properties(mpfr_imported PROPERTIES
    IMPORTED_LOCATION ${MPFR_LIBRARY}
    INTERFACE_INCLUDE_DIRECTORIES ${MPFR_INCLUDE_DIR}
)

# Create imported target for FLINT
add_library(flint_imported STATIC IMPORTED GLOBAL)
set_target_properties(flint_imported PROPERTIES
    IMPORTED_LOCATION ${DEPS_LIB_DIR}/libflint.a
    INTERFACE_INCLUDE_DIRECTORIES ${DEPS_INCLUDE_DIR}
)
add_dependencies(flint_imported ep_flint)

# Create an interface target that bundles everything
add_library(flint_external INTERFACE)
target_include_directories(flint_external INTERFACE 
    ${DEPS_INCLUDE_DIR}
    ${GMP_INCLUDE_DIR}
    ${MPFR_INCLUDE_DIR}
)
target_link_libraries(flint_external INTERFACE
    flint_imported
    mpfr_imported
    gmpxx_imported
    gmp_imported
)

set(FLINT_TARGET flint_external)
set(FLINT_INCLUDE_DIRS ${DEPS_INCLUDE_DIR})