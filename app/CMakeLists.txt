## @file app/CMakeLists.txt
## @brief CMake build script for applications based on stencil_rtm.
##
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project (stencil_apps LANGUAGES C CXX)

## Include the stencil_rtm headers
include_directories(${CMAKE_SOURCE_DIR}/include)

### Check whether CUDA is detected:
### - The CUDA library queered in usual directories
### - This behavior can be changed with env setting:
###   CUDA_INCLUDE_DIRS  : contains the CUDA headers
###   CUDA_LIBRARIES_DIR : contains the CUDA runtime libraries
#find_package (CUDA REQUIRED)
#
### include headers.
#include_directories(${CUDA_INCLUDE_DIRS})
#
### add CUDA compile options.
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo ")

## wildcard source files.
file(GLOB C_SRCS "*.c")

## generate the executables.
add_custom_target(stencil_apps)

## set the compiler flags for each NV target
foreach(app_src ${C_SRCS})
  file(RELATIVE_PATH app_rel ${CMAKE_CURRENT_SOURCE_DIR} ${app_src})
  string(REPLACE ".c" "" app_exe ${app_rel})
  add_executable(${app_exe} ${app_src})
  target_link_libraries(${app_exe} stencil)
  target_link_libraries(${app_exe} m)
  target_link_libraries(${app_exe} stdc++)
  target_link_libraries(${app_exe} ${OpenMP_C_LIBRARIES})
  target_link_libraries(${app_exe} ${OpenMP_CXX_LIBRARIES})
  target_link_libraries(${app_exe} ${OpenMP_C_FLAGS})
  target_link_libraries(${app_exe} ${OpenMP_CXX_FLAGS})
  add_dependencies(stencil_apps ${app_exe})
  install(TARGETS ${app_exe} RUNTIME DESTINATION bin)
endforeach(app_src ${C_SRCS})
