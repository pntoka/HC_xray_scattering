cmake_minimum_required(VERSION 3.15)
project(iobs_ngc LANGUAGES CXX C)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
message(STATUS "Found pybind11: ${pybind11_VERSION}")

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
if(USE_SYSTEM_LIBS)
    message(STATUS "Using system-installed FLINT")
    find_package(FLINT REQUIRED)
    
    # Create interface target for consistency
    add_library(flint_external INTERFACE)
    target_include_directories(flint_external INTERFACE ${FLINT_INCLUDE_DIRS})
    target_link_libraries(flint_external INTERFACE ${FLINT_LIBRARIES})
    set(FLINT_TARGET flint_external)
else()
    message(STATUS "Building FLINT and dependencies from source")
    include(BuildDependencies)
    # FLINT_TARGET is set by BuildDependencies.cmake
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

# Force complete static linking on macOS
if(APPLE)
    target_link_options(iobs_ngc PRIVATE
        -Wl,-force_load,${DEPS_LIB_DIR}/libflint.a
        -Wl,-force_load,${DEPS_LIB_DIR}/libmpfr.a
        -Wl,-force_load,${DEPS_LIB_DIR}/libgmp.a
    )
elseif(UNIX)
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