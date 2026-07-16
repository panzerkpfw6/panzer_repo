##
## @file app/CMakeLists.txt
## @brief Build applications based on stencil.
##

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

include_directories("${CMAKE_SOURCE_DIR}/include")

add_custom_target(stencil_apps)

# ------------------------------------------------------------
# gather: CPU only
# ------------------------------------------------------------
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
install(TARGETS gather RUNTIME DESTINATION bin)

# ------------------------------------------------------------
# modeling: CPU plus optional CUDA
# ------------------------------------------------------------
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

    set_target_properties(modeling PROPERTIES
        LINKER_LANGUAGE C
        CUDA_RESOLVE_DEVICE_SYMBOLS OFF
    )
endif()

add_dependencies(stencil_apps modeling)
install(TARGETS modeling RUNTIME DESTINATION bin)

# ------------------------------------------------------------
# rtm: initially CPU only
# ------------------------------------------------------------
add_executable(rtm
    rtm.c
)

target_link_libraries(rtm
    PRIVATE
        stencil
        OpenMP::OpenMP_C
        m
)

# Add stencil_cuda here only when rtm actually uses GPU functions.
#
# if(USE_CUDA)
#     target_compile_definitions(rtm PRIVATE USE_CUDA)
#     target_link_libraries(rtm PRIVATE stencil_cuda CUDA::cudart)
# endif()

add_dependencies(stencil_apps rtm)
install(TARGETS rtm RUNTIME DESTINATION bin)