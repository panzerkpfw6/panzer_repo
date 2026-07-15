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

project(stencil LANGUAGES C CXX CUDA)

find_package(CUDAToolkit REQUIRED)
#find_package (CUDA REQUIRED)
find_package(OpenMP REQUIRED)

option(USE_CUDA "Enable CUDA support" ON)

## include headers.
include_directories(${CMAKE_SOURCE_DIR}/include)

## compile C and CUDA files.
file(GLOB C_SRCS CONFIGURE_DEPENDS "*.c")
file(GLOB CUDA_SRCS CONFIGURE_DEPENDS "*.cu")

#message ("CMAKE_CUDA_FLAGS : " ${CMAKE_CUDA_FLAGS})

add_library(stencil
    ${C_SRCS}
    ${CUDA_SRCS}
)

if(USE_CUDA)
    target_compile_definitions(stencil PUBLIC USE_CUDA)

    set_target_properties(stencil PROPERTIES
        CUDA_ARCHITECTURES 70
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_STANDARD 14
        CUDA_STANDARD_REQUIRED ON
    )

    target_link_libraries(stencil PUBLIC CUDA::cudart)
endif()

target_link_libraries(stencil PUBLIC
    OpenMP::OpenMP_C
)

target_compile_options(stencil PRIVATE
    $<$<COMPILE_LANGUAGE:C>:-O3;-march=core-avx2;-mtune=core-avx2>
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-march=core-avx2;-mtune=core-avx2>
    $<$<COMPILE_LANGUAGE:CUDA>:-O3>
)

# add_compile_options(
#     $<$<COMPILE_LANGUAGE:C>:-qopenmp>
#     $<$<COMPILE_LANGUAGE:C>:-march=core-avx2>
#     $<$<COMPILE_LANGUAGE:C>:-mtune=core-avx2>

#     $<$<COMPILE_LANGUAGE:CXX>:-qopenmp>
#     $<$<COMPILE_LANGUAGE:CXX>:-march=core-avx2>
#     $<$<COMPILE_LANGUAGE:CXX>:-mtune=core-avx2>
# )

install(TARGETS stencil ARCHIVE DESTINATION lib)
