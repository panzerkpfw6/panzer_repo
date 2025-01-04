## @brief CMake main build script for stencil_rtm.
##
cmake_minimum_required (VERSION 3.8 FATAL_ERROR)

project (stencil_distribution LANGUAGES C)

## overwrite the install prefix.
set (CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")

## add extra cmake modules dir.
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

## configurable options for how we want to build project.
#option(project_debug "Build project with the debug mode." OFF)
#option(project_debug "Enable debugging" ON)

## set the build type (DEFAULT is Release).
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)
if (project_debug)
    set (CMAKE_BUILD_TYPE Debug)
endif (project_debug)

## check whether CUDA is detected:
## - The CUDA library queered in usual directories
## - This behavior can be changed with environment setting:
##   CUDA_INCLUDE_DIRS  : contains the CUDA headers
##   CUDA_LIBRARIES_DIR : contains the CUDA runtime libraries
#find_package (CUDA REQUIRED)

## add CUDA arch options.
#include(${CMAKE_MODULE_PATH}/cuda_arch_options.cmake)
#if (project_debug)
#  set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -D__DEBUG")
#endif (project_debug)

### set the compilation flags.
#if ("${CMAKE_C_COMPILER_ID}"     MATCHES "Cray")
#    include(${CMAKE_MODULE_PATH}/cray_c.cmake)
#    #elseif ("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
#    #  include(${CMAKE_MODULE_PATH}/gnu_c.cmake)
#elseif ("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
#    include(${CMAKE_MODULE_PATH}/intel_c.cmake)
#elseif ("${CMAKE_C_COMPILER_ID}" MATCHES "PGI")
#    include(${CMAKE_MODULE_PATH}/pgi_c.cmake)
#elseif ("${CMAKE_C_COMPILER_ID}" MATCHES "MSVC")
#    include(${CMAKE_MODULE_PATH}/msvc_c.cmake)
#endif ("${CMAKE_C_COMPILER_ID}"  MATCHES "Cray")

### set the compilation flags.
add_compile_options(-qopenmp -march=core-avx2 -mtune=core-avx2)
#add_compile_options(-fopenmp -march=core-avx2 -mtune=core-avx2)
#add_compile_options(-fopenmp -g3 -O0 -march=core-avx2 -mtune=core-avx2) # for debug

## add OpenMP.
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS
            "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

## build the project C/C++ library.
add_subdirectory(src)

## build the project based applications.
add_subdirectory(app)

## build the project unit and functional tests.
#add_subdirectory(test)
