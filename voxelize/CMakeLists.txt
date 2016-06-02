cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_BUILD_TYPE Debug)
project(voxelize)

find_package(PCL 1.2 REQUIRED)
#SET(CMAKE_CXX_FLAGS "-std=c++0x -Wno-deprecated")

include_directories(${PCL_INCLUDE_DIRS} /home/adamjri/mesh_reconstruction/binvoxToPCL)
include_directories(${PCL_INCLUDE_DIRS} /home/adamjri/mesh_reconstruction/confidence_pcl)

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
add_library(voxelize_lib SHARED voxelize.h voxelize.cpp)

target_link_libraries (voxelize_lib  binvoxToPcl_lib assignConfidence_lib ${PCL_LIBRARIES})
add_executable (voxelize main.cpp)

target_link_libraries (voxelize voxelize_lib ${PCL_LIBRARIES})
