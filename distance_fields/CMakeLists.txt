cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_BUILD_TYPE Debug)

project(dfields)

find_package(CUDA)
list(APPEND CUDA_NVCC_FLAGS "-Xcudafe --diag_suppress=boolean_controlling_expr_is_constant")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../grid)

if(CUDA_FOUND)
cuda_add_library(dfields_lib SHARED dfields.h dfields.cu)
else()
add_library(dfields_lib SHARED dfields.h dfields.cpp)
endif(CUDA_FOUND)

target_link_libraries (dfields_lib grid_lib)
if(CUDA_FOUND)
cuda_add_executable (dfields main.cpp)
else()
add_executable(dfields main.cpp)
endif(CUDA_FOUND)

target_link_libraries (dfields dfields_lib)


install(TARGETS dfields_lib 
        RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)