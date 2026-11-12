##
## @file src/CMakeLists.txt
## @brief Build CPU and CUDA stencil libraries separately.
##

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

## CPU source files only.
file(GLOB C_SRCS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.c"
)

add_library(stencil STATIC
    ${C_SRCS}
)

target_include_directories(stencil
    PUBLIC
        "${CMAKE_SOURCE_DIR}/include"
)

target_compile_options(stencil
    PRIVATE
        $<$<COMPILE_LANGUAGE:C>:-O3;-march=core-avx2;-mtune=core-avx2>
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-march=core-avx2;-mtune=core-avx2>
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

# ------------------------------------------------------------
# AMD HIP implementation
# ------------------------------------------------------------
if(USE_HIP)
    add_library(stencil_hip STATIC
        gpu_wave_SB_HIP.hip
    )

    target_include_directories(stencil_hip
        PUBLIC
            "${CMAKE_SOURCE_DIR}/include"
    )

    target_compile_definitions(stencil_hip
        PRIVATE USE_HIP
    )

    set_target_properties(stencil_hip PROPERTIES
        HIP_STANDARD 17
        HIP_STANDARD_REQUIRED ON
        HIP_ARCHITECTURES gfx942
        POSITION_INDEPENDENT_CODE ON
    )

    install(
        TARGETS stencil_hip
        ARCHIVE DESTINATION lib
    )
endif()

## CUDA implementation in a separate library.
if(USE_CUDA)
    add_library(stencil_cuda STATIC
        gpu_wave.cu
        modeling_gpu.cu
        gpu_wave_tb.cu
        modeling_tb_gpu.cu
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
        CUDA_STANDARD 14
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
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