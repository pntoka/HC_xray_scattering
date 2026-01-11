cmake_minimum_required(VERSION 3.15)
project(iobs_ngc)

# Add cmake modules path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Find Python
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Find pybind11 (from pip install)
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
    OUTPUT_VARIABLE pybind11_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
find_package(pybind11 REQUIRED)

# Find FLINT using our custom finder
find_package(FLINT REQUIRED)

# Find Boost (header-only math library)
find_package(Boost REQUIRED)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Sources
set(SOURCES
    src/calculations.cpp
    src/iobs_calculator.cpp
    src/python_bindings.cpp
)

# Create Python module
pybind11_add_module(iobs_ngc ${SOURCES})

target_include_directories(iobs_ngc PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${FLINT_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

target_compile_definitions(iobs_ngc PRIVATE
    _USE_MATH_DEFINES
)

target_link_libraries(iobs_ngc PRIVATE
    ${FLINT_LIBRARIES}
    ${Boost_LIBRARIES}
)

# Install target
install(TARGETS iobs_ngc DESTINATION .)