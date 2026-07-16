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
## @brief CMake build script for the stencil libraries.
##

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

find_package(OpenMP REQUIRED)

## CPU sources only.
file(GLOB C_SRCS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.c"
)

## Main CPU stencil library.
add_library(stencil STATIC
    ${C_SRCS}
)

target_include_directories(stencil
    PUBLIC
        "${CMAKE_SOURCE_DIR}/include"
)

target_compile_options(stencil
    PRIVATE
        $<$<COMPILE_LANGUAGE:C>:-O3>
        $<$<COMPILE_LANGUAGE:C>:-march=core-avx2>
        $<$<COMPILE_LANGUAGE:C>:-mtune=core-avx2>
)

target_link_libraries(stencil
    PUBLIC
        OpenMP::OpenMP_C
        m
)

install(
    TARGETS stencil
    ARCHIVE DESTINATION lib
)

## Separate CUDA library.
if(USE_CUDA)
    find_package(CUDAToolkit REQUIRED)

    add_library(stencil_cuda STATIC
        gpu_wave.cu
    )

    target_include_directories(stencil_cuda
        PUBLIC
            "${CMAKE_SOURCE_DIR}/include"
    )

    target_compile_definitions(stencil_cuda
        PRIVATE
            USE_CUDA
    )

    target_compile_options(stencil_cuda
        PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:-O3>
    )

    set_target_properties(stencil_cuda PROPERTIES
        CUDA_ARCHITECTURES 70
        CUDA_STANDARD 14
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION ON
        POSITION_INDEPENDENT_CODE ON
    )

    target_link_libraries(stencil_cuda
        PRIVATE
            CUDA::cudart
    )

    install(
        TARGETS stencil_cuda
        ARCHIVE DESTINATION lib
    )
endif()