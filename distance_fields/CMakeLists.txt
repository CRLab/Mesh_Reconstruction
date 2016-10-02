cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_BUILD_TYPE Debug)

project(dfields)

find_package(CUDA)
list(APPEND CUDA_NVCC_FLAGS "-Xcudafe --diag_suppress=boolean_controlling_expr_is_constant")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../grid)

cuda_add_library(dfields_lib SHARED dfields.h dfields.cu)

target_link_libraries (dfields_lib grid_lib)
cuda_add_executable (dfields main.cpp)

target_link_libraries (dfields dfields_lib)

