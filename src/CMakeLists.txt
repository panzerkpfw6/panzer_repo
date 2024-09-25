##
## @copyright Copyright 2024- Pavel Plotnitskii. All rights reserved.
## This file is part of the \b stencil project.
##
## \b stencil is free software: you may redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## The stencil project is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with the \b stencil project. If not, see <http://www.gnu.org/licenses/>.
##
## @file src/CMakeLists.txt
## @author Pavel Plotnitskii
## @brief CMake build script for the stencil executable.
##
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(stencil LANGUAGES C)

## Check whether CUDA is detected:
## - The CUDA library queered in usual directories
## - This behavior can be changed with environment setting:
##   CUDA_INCLUDE_DIRS  : contains the CUDA headers
##   CUDA_LIBRARIES_DIR : contains the CUDA runtime libraries
#find_package (CUDA REQUIRED)

## include headers.
include_directories(${CMAKE_SOURCE_DIR}/include)

## compile C and CUDA files.
file(GLOB SRCS "*.c")

#message ("CMAKE_CUDA_FLAGS : " ${CMAKE_CUDA_FLAGS})

add_library(stencil ${SRCS})
#set_target_properties(stencil PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

## install the library.
install(TARGETS stencil ARCHIVE DESTINATION lib)
