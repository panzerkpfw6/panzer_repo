##
## @file app/CMakeLists.txt
## @brief CMake build script for applications based on stencil_rtm.
##

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(stencil_apps LANGUAGES C CXX)

## Include project headers.
include_directories(${CMAKE_SOURCE_DIR}/include)

## Convenience target that builds all applications.
add_custom_target(stencil_apps)

# ============================================================
# gather: CPU-only application
# ============================================================

add_executable(gather
    gather.c
)

target_link_libraries(gather
    PRIVATE
        stencil
        OpenMP::OpenMP_C
        m
)

add_dependencies(stencil_apps gather)

install(
    TARGETS gather
    RUNTIME DESTINATION bin
)

# ============================================================
# modeling: CPU + optional CUDA application
# ============================================================

add_executable(modeling
    modeling.c
)

target_link_libraries(modeling
    PRIVATE
        stencil
        OpenMP::OpenMP_C
        m
)

if(USE_CUDA)
    target_compile_definitions(modeling
        PRIVATE
            USE_CUDA
    )

    target_link_libraries(modeling
        PRIVATE
            stencil_cuda
            CUDA::cudart
    )
endif()

add_dependencies(stencil_apps modeling)

install(
    TARGETS modeling
    RUNTIME DESTINATION bin
)

# ============================================================
# rtm: CPU + optional CUDA application
# ============================================================

add_executable(rtm
    rtm.c
)

target_link_libraries(rtm
    PRIVATE
        stencil
        OpenMP::OpenMP_C
        m
)

if(USE_CUDA)
    target_compile_definitions(rtm
        PRIVATE
            USE_CUDA
    )

    target_link_libraries(rtm
        PRIVATE
            stencil_cuda
            CUDA::cudart
    )
endif()

add_dependencies(stencil_apps rtm)

install(
    TARGETS rtm
    RUNTIME DESTINATION bin
)