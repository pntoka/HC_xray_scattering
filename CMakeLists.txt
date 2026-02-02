cmake_minimum_required(VERSION 3.15)
project(iobs_ngc LANGUAGES CXX C)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Debug: Print architecture information
message(STATUS "CMake Generator: ${CMAKE_GENERATOR}")
message(STATUS "CMake Generator Platform: ${CMAKE_GENERATOR_PLATFORM}")
message(STATUS "Pointer size: ${CMAKE_SIZEOF_VOID_P} bytes")
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(STATUS "Building for 64-bit (x64)")
else()
    message(STATUS "Building for 32-bit (Win32)")
endif()

# ============================================================================
# Options
# ============================================================================
option(USE_SYSTEM_LIBS "Use system-installed FLINT/GMP/MPFR instead of building from source" OFF)
option(FETCH_BOOST "Download Boost via FetchContent instead of using system Boost" ON)

# Add cmake modules path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# ============================================================================
# Conda environment support (only when using system libs)
# ============================================================================
if(USE_SYSTEM_LIBS)
    if(DEFINED ENV{CONDA_PREFIX})
        list(PREPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
        message(STATUS "Added CONDA_PREFIX to CMAKE_PREFIX_PATH: $ENV{CONDA_PREFIX}")
    endif()
endif()

# ============================================================================
# Python and pybind11
# ============================================================================
# Help CMake find Python in manylinux environments
if(DEFINED Python3_ROOT_DIR AND NOT Python3_INCLUDE_DIRS)
    set(Python3_FIND_STRATEGY LOCATION)
    set(Python3_FIND_REGISTRY NEVER)
endif()

find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)

# Try to find pybind11 via CMake config first, fall back to Python module
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    message(STATUS "pybind11 not found via CMake config, trying Python module...")
    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -m pybind11 --cmakedir
        OUTPUT_VARIABLE pybind11_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE pybind11_RESULT
    )
    if(pybind11_RESULT EQUAL 0)
        message(STATUS "Found pybind11 CMake dir: ${pybind11_DIR}")
        find_package(pybind11 CONFIG REQUIRED)
    else()
        message(FATAL_ERROR "pybind11 not found. Please install it via pip.")
    endif()
endif()

message(STATUS "Found pybind11: ${pybind11_VERSION}")
message(STATUS "Python3 executable: ${Python3_EXECUTABLE}")
message(STATUS "Python3 include dirs: ${Python3_INCLUDE_DIRS}")

# ============================================================================
# Boost
# ============================================================================
if(FETCH_BOOST)
    include(FetchContent)
    message(STATUS "Fetching Boost via FetchContent...")
    
    # Only include the Boost libraries we need (header-only for multiprecision)
    set(BOOST_INCLUDE_LIBRARIES multiprecision math)
    set(BOOST_ENABLE_CMAKE ON)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    
    FetchContent_Declare(
        Boost
        URL https://github.com/boostorg/boost/releases/download/boost-1.84.0/boost-1.84.0.tar.xz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        EXCLUDE_FROM_ALL
    )
    FetchContent_MakeAvailable(Boost)
    
    set(BOOST_TARGET Boost::multiprecision)
    message(STATUS "Boost will be fetched from source")
else()
    # Use system Boost (conda or system-installed)
    find_package(Boost REQUIRED)
    message(STATUS "Found system Boost: ${Boost_VERSION}")
    set(BOOST_TARGET Boost::headers)
endif()

# ============================================================================
# FLINT and dependencies (GMP, MPFR)
# ============================================================================
# Check if we should use pre-built FLINT from CI (with system GMP/MPFR)
if(DEFINED ENV{CI} AND EXISTS "/opt/flint-deps/lib/libflint.a")
    message(STATUS "Using pre-built FLINT from CI with system GMP/MPFR")
    
    set(FLINT_INSTALL_DIR "/opt/flint-deps")
    
    # Find system GMP and MPFR
    find_library(GMP_LIBRARY NAMES gmp REQUIRED)
    find_library(GMPXX_LIBRARY NAMES gmpxx REQUIRED)
    find_library(MPFR_LIBRARY NAMES mpfr REQUIRED)
    find_path(GMP_INCLUDE_DIR NAMES gmp.h REQUIRED)
    find_path(MPFR_INCLUDE_DIR NAMES mpfr.h REQUIRED)
    
    message(STATUS "Found system GMP: ${GMP_LIBRARY}")
    message(STATUS "Found system MPFR: ${MPFR_LIBRARY}")
    
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

    # Create imported target for pre-built FLINT
    add_library(flint_imported STATIC IMPORTED GLOBAL)
    set_target_properties(flint_imported PROPERTIES
        IMPORTED_LOCATION ${FLINT_INSTALL_DIR}/lib/libflint.a
        INTERFACE_INCLUDE_DIRECTORIES ${FLINT_INSTALL_DIR}/include
    )

    # Create an interface target that bundles everything
    add_library(flint_external INTERFACE)
    target_include_directories(flint_external INTERFACE 
        ${FLINT_INSTALL_DIR}/include
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
    set(USE_SYSTEM_LIBS OFF)  # Override since we're using our custom setup
elseif(USE_SYSTEM_LIBS)
    message(STATUS "Using system-installed FLINT")
    
    if(WIN32 AND DEFINED CMAKE_TOOLCHAIN_FILE AND CMAKE_TOOLCHAIN_FILE MATCHES "vcpkg")
        message(STATUS "Using vcpkg for Windows dependencies")
        
        find_package(PkgConfig REQUIRED)
        
        pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)
        pkg_check_modules(gmpxx REQUIRED IMPORTED_TARGET gmpxx)
        pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
        pkg_check_modules(flint REQUIRED IMPORTED_TARGET flint)
        
        add_library(flint_external INTERFACE)
        target_link_libraries(flint_external INTERFACE 
            PkgConfig::flint
            PkgConfig::mpfr
            PkgConfig::gmpxx
            PkgConfig::gmp
        )
        set(FLINT_TARGET flint_external)
    else()
        find_package(FLINT REQUIRED)
        
        add_library(flint_external INTERFACE)
        target_include_directories(flint_external INTERFACE ${FLINT_INCLUDE_DIRS})
        target_link_libraries(flint_external INTERFACE ${FLINT_LIBRARIES})
        set(FLINT_TARGET flint_external)
    endif()
else()
    message(STATUS "Building FLINT from source (GMP/MPFR from vcpkg on Windows)")
    include(BuildDependencies)
endif()

# ============================================================================
# Sources
# ============================================================================
set(SOURCES
    src/calculations.cpp
    src/iobs_calculator.cpp
    src/python_bindings.cpp
)

# ============================================================================
# Create the Python module
# ============================================================================
pybind11_add_module(iobs_ngc ${SOURCES})

# Add explicit dependency on FLINT build
if(NOT USE_SYSTEM_LIBS)
    add_dependencies(iobs_ngc ep_flint)
endif()

# Include directories
target_include_directories(iobs_ngc PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

# Add FLINT include directory when building from source
if(NOT USE_SYSTEM_LIBS)
    target_include_directories(iobs_ngc PRIVATE
        ${CMAKE_BINARY_DIR}/deps/include
    )
endif()

# Link libraries
target_link_libraries(iobs_ngc PRIVATE
    ${FLINT_TARGET}
    ${BOOST_TARGET}
)

# Platform-specific linking
if(WIN32)
    # Windows: link statically
    target_link_libraries(iobs_ngc PRIVATE ${FLINT_TARGET})
elseif(APPLE)
    target_link_options(iobs_ngc PRIVATE
        -Wl,-force_load,${DEPS_LIB_DIR}/libflint.a
        -Wl,-force_load,${DEPS_LIB_DIR}/libmpfr.a
        -Wl,-force_load,${DEPS_LIB_DIR}/libgmp.a
    )
elseif(UNIX AND NOT APPLE)
    target_link_options(iobs_ngc PRIVATE
        -Wl,--whole-archive
        ${DEPS_LIB_DIR}/libflint.a
        ${DEPS_LIB_DIR}/libmpfr.a
        ${DEPS_LIB_DIR}/libgmp.a
        -Wl,--no-whole-archive
    )
endif()

# Add math defines
target_compile_definitions(iobs_ngc PRIVATE
    _USE_MATH_DEFINES
)

# Platform-specific definitions
if(WIN32)
    target_compile_definitions(iobs_ngc PRIVATE USING_ARB)
    message(STATUS "Windows: Using ARB library")
else()
    target_compile_definitions(iobs_ngc PRIVATE USING_FLINT)
    message(STATUS "Unix: Using FLINT 3+ library")
endif()

# ============================================================================
# Compiler warnings
# ============================================================================
if(MSVC)
    target_compile_options(iobs_ngc PRIVATE /W4)
else()
    target_compile_options(iobs_ngc PRIVATE -Wall -Wextra -Wpedantic)
endif()

# ============================================================================
# Install
# ============================================================================
install(TARGETS iobs_ngc DESTINATION .)

# ============================================================================
# Summary
# ============================================================================
message(STATUS "")
message(STATUS "=== iobs_ngc Configuration Summary ===")
message(STATUS "  USE_SYSTEM_LIBS: ${USE_SYSTEM_LIBS}")
message(STATUS "  FETCH_BOOST:     ${FETCH_BOOST}")
message(STATUS "  Python:          ${Python_EXECUTABLE}")
message(STATUS "  pybind11:        ${pybind11_VERSION}")
message(STATUS "======================================")
message(STATUS "")