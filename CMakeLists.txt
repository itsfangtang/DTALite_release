cmake_minimum_required(VERSION 3.14)
project(TAPLite VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/TAPLite/core)

# Build Type and Optimization
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "/O2 /GL")
set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od /RTC1")

# Specify output directories for binaries and libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/output/lib)  # Static libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/output/lib)  # Shared libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/output/bin)  # Executables
set(CMAKE_PYTHON_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/output/python) # Python bindings

# Find pybind11
find_package(pybind11 REQUIRED)

# Optional: OpenMP support
option(USE_OPENMP "Enable OpenMP support" ON)
if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_FOUND)
        message(STATUS "OpenMP found, enabling OpenMP support")
        set(OPENMP_FLAGS ${OpenMP_CXX_FLAGS})
        set(OPENMP_LIBRARIES OpenMP::OpenMP_CXX)
    else()
        message(WARNING "OpenMP requested but not found")
    endif()
endif()

# Source files
set(TAPLite_SOURCES
    TAPLite/bindings.cpp
    TAPLite/core/TAPLite.cpp
    TAPLite/core/TAPLite.h
)

# Define the Python extension module
pybind11_add_module(TAPLite ${TAPLite_SOURCES})

# Specify where the compiled Python module should go
set_target_properties(TAPLite PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_PYTHON_OUTPUT_DIRECTORY})

# Compiler flags and warning suppression
if(MSVC)
    target_compile_options(TAPLite PRIVATE /W3 /bigobj /wd4477 /wd4267 /wd4244 /wd4996 /wd4101)
else()
    target_compile_options(TAPLite PRIVATE -Wall -Wextra -Wpedantic)
endif()

# OpenMP support
if(USE_OPENMP AND OPENMP_FOUND)
    target_compile_options(TAPLite PRIVATE ${OPENMP_FLAGS})
    target_link_libraries(TAPLite PRIVATE ${OPENMP_LIBRARIES})
endif()

# Installation
install(
    TARGETS TAPLite
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.*/site-packages
)

# Project summary
message(STATUS "Project Configuration Summary:")
message(STATUS "  Project Name: ${PROJECT_NAME}")
message(STATUS "  Project Version: ${PROJECT_VERSION}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
if(USE_OPENMP)
    if(OPENMP_FOUND)
        message(STATUS "  OpenMP: Enabled")
    else()
        message(STATUS "  OpenMP: Disabled (not found)")
    endif()
else()
    message(STATUS "  OpenMP: Disabled (by user)")
endif()
