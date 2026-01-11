cmake_minimum_required(VERSION 3.15)
project(iobs_ngc LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Add conda environment to CMAKE_PREFIX_PATH if available
if(DEFINED ENV{CONDA_PREFIX})
    list(PREPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
    message(STATUS "Added CONDA_PREFIX to CMAKE_PREFIX_PATH: $ENV{CONDA_PREFIX}")
endif()

# Add cmake modules path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Find Python
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Find pybind11 (installed via conda)
find_package(pybind11 CONFIG REQUIRED)
message(STATUS "Found pybind11: ${pybind11_VERSION}")

# Find Boost
find_package(Boost REQUIRED)
message(STATUS "Found Boost: ${Boost_VERSION}")

# Find FLINT using our custom finder
find_package(FLINT REQUIRED)

# Sources
set(SOURCES
    src/bindings.cpp
    # Add other source files here:
    # src/other_file.cpp
)

# Create the Python module
pybind11_add_module(iobs_ngc ${SOURCES})

# Include directories
target_include_directories(iobs_ngc PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${FLINT_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# Link libraries
target_link_libraries(iobs_ngc PRIVATE
    ${FLINT_LIBRARIES}
)

# Compiler warnings (optional but recommended)
if(MSVC)
    target_compile_options(iobs_ngc PRIVATE /W4)
else()
    target_compile_options(iobs_ngc PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Install the module
install(TARGETS iobs_ngc DESTINATION .)