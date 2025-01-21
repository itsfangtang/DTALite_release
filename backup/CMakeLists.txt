cmake_minimum_required(VERSION 3.5.0)
project(TAPLite VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE Release)
add_executable(TAPLite TAPLite/TAPLite.cpp)

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