cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_BUILD_TYPE Debug)
project(grid)

find_package(PCL 1.2 REQUIRED)
#SET(CMAKE_CXX_FLAGS "-std=c++0x -Wno-deprecated")

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
add_library(grid_lib SHARED grid.h grid.cpp)

target_link_libraries (grid_lib ${PCL_LIBRARIES})
add_executable (grid main.cpp)

target_link_libraries (grid grid_lib)
