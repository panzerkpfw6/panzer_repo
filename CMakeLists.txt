##
## @copyright Copyright 2017- Issam Said. All rights reserved.
## This file is part of \b simwave.
##
## \b simwave is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## simwave is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with \b simwave.  If not, see <http://www.gnu.org/licenses/>.
##
## @file CMakeLists.txt
## @author Issam SAID
## @brief CMake main build script for simwave.
##
cmake_minimum_required (VERSION 3.8 FATAL_ERROR)

project (simwave_distribution LANGUAGES C)

## overwrite the install prefix.
set (CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")

## add extra cmake modules dir.
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

## configurable options for how we want to build simwave.
option(simwave_debug "Build simwave with the debug mode." OFF)

## set the build type (DEFAULT is Release).
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)
if (simwave_debug)
	set (CMAKE_BUILD_TYPE Debug)
endif (simwave_debug)

## check whether CUDA is detected:
## - The CUDA library queered in usual directories
## - This behavior can be changed with environment setting:
##   CUDA_INCLUDE_DIRS  : contains the CUDA headers
##   CUDA_LIBRARIES_DIR : contains the CUDA runtime libraries
#find_package (CUDA REQUIRED)

## set the compilation flags. 
if ("${CMAKE_C_COMPILER_ID}"     MATCHES "Cray")
  include(${CMAKE_MODULE_PATH}/cray_c.cmake)
#elseif ("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
#  include(${CMAKE_MODULE_PATH}/gnu_c.cmake)
elseif ("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  include(${CMAKE_MODULE_PATH}/intel_c.cmake)
elseif ("${CMAKE_C_COMPILER_ID}" MATCHES "PGI")
  include(${CMAKE_MODULE_PATH}/pgi_c.cmake)
elseif ("${CMAKE_C_COMPILER_ID}" MATCHES "MSVC")
 	include(${CMAKE_MODULE_PATH}/msvc_c.cmake)
endif ("${CMAKE_C_COMPILER_ID}"  MATCHES "Cray")


## add CUDA arch options.
#include(${CMAKE_MODULE_PATH}/cuda_arch_options.cmake)
#if (simwave_debug)
#  set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -D__DEBUG")
#endif (simwave_debug)

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

## build the simwave C/C++ library.
add_subdirectory(src)

## build the simwave based applications.
add_subdirectory(app)

## build the simwave unit and functional tests.
#add_subdirectory(test)
