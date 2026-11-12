## @brief CMake main build script for stencil_rtm.

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(stencil_distribution LANGUAGES C CXX)

## Install into the repository tree.
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}" CACHE PATH
    "Install path prefix" FORCE)

message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")

## Extra CMake modules.
list(APPEND CMAKE_MODULE_PATH
    "${CMAKE_SOURCE_DIR}/cmake/modules"
)

message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

## Build options.
option(USE_CUDA "Enable CUDA support" OFF)
option(USE_HIP "Enable AMD HIP support" OFF)
option(PROJECT_DEBUG "Build project in debug mode" OFF)

## Default build type.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Build type" FORCE)
endif()

if(PROJECT_DEBUG)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Build type" FORCE)
endif()

## Dependencies.
find_package(OpenMP REQUIRED)

if(USE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)

    if(NOT CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 70)
    endif()
endif()

if(USE_HIP)
    enable_language(HIP)

    if(NOT CMAKE_HIP_ARCHITECTURES)
        set(CMAKE_HIP_ARCHITECTURES gfx942)
    endif()
endif()

## Build CPU/GPU library.
add_subdirectory(src)

## Build applications.
add_subdirectory(app)

## Optional tests.
# add_subdirectory(test)