cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_BUILD_TYPE Debug)

project(dfields)

include_directories(/home/adamjri/mesh_reconstruction/grid)

add_library(dfields_lib SHARED dfields.h dfields.cpp)

target_link_libraries (dfields_lib grid_lib)
add_executable (dfields main.cpp)

target_link_libraries (dfields dfields_lib)

