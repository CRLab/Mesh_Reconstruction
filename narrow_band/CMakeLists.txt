cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_BUILD_TYPE Debug)
project(narrowBand)

find_package(PCL 1.2 REQUIRED)
#SET(CMAKE_CXX_FLAGS "-std=c++0x -Wno-deprecated")

include_directories(${PCL_INCLUDE_DIRS} /home/adamjri/mesh_reconstruction/binvoxToPCL)
include_directories(${PCL_INCLUDE_DIRS} /home/adamjri/mesh_reconstruction/confidence_pcl)
include_directories(${PCL_INCLUDE_DIRS} /home/adamjri/mesh_reconstruction/voxelize)

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
add_library(narrowBand_lib SHARED grid.h grid.cpp narrowBand.h narrowBand.cpp getsqdist.h getsqdist.cpp)

target_link_libraries (narrowBand_lib  binvoxToPcl_lib assignConfidence_lib voxelize_lib ${PCL_LIBRARIES})
add_executable (narrowBand main.cpp)

target_link_libraries (narrowBand narrowBand_lib ${PCL_LIBRARIES})
