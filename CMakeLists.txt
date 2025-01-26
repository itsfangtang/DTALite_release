cmake_minimum_required (VERSION 3.8)

project ("TAPLite")

option(BUILD_EXE "BUILD Executable or not" ON)

if(WIN32)
    message("Building dll file on Windows")
elseif(UNIX)
    message("Building so or dylib file on ${CMAKE_SYSTEM_NAME}")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_OSX_ARCHITECTURES "arm64")

# note CMAKE_BUILD_TYPE WILL BE INGNORED for IDE generators, like Visual Studio and Xcode
# it only works for single-configuration generators, e.g., make and Ninja
# see https://cmake.org/cmake/help/v3.0/variable/CMAKE_BUILD_TYPE.html?highlight=cmake_build_types for details
set(CMAKE_BUILD_TYPE  "Release")
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/bin)


if (BUILD_EXE)
    add_executable(${PROJECT_NAME} TAPLite/core/TAPLite.cpp)
else()
    add_library(${PROJECT_NAME} SHARED TAPLite/core/TAPLite.cpp)
endif()


if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        set(CMAKE_CXX_FLAGS "-Xpreprocessor -fopenmp -fpic")
    else()
        set(CMAKE_CXX_FLAGS "-fopenmp -fpic")
    endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CMAKE_CXX_FLAGS "-fopenmp -fpic")
endif()

find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)