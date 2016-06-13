cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(CMAKE_BUILD_TYPE Debug)
project(optimize)

find_package(PCL 1.2 REQUIRED)
#SET(CMAKE_CXX_FLAGS "-std=c++0x -Wno-deprecated")

include_directories(${PCL_INCLUDE_DIRS} /home/adamjri/mesh_reconstruction/binvoxToPCL)
include_directories(${PCL_INCLUDE_DIRS} /home/adamjri/mesh_reconstruction/confidence_pcl)
include_directories(${PCL_INCLUDE_DIRS} /home/adamjri/mesh_reconstruction/voxelize)
include_directories(${PCL_INCLUDE_DIRS} /home/adamjri/mesh_reconstruction/narrow_band)


link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
add_library(optimize_lib SHARED quadprog.h quadprog.cpp primeqp.h primeqp.cpp)

target_link_libraries (optimize_lib narrowBand_lib  binvoxToPcl_lib assignConfidence_lib voxelize_lib ${PCL_LIBRARIES})
add_executable (optimize main.cpp)

target_link_libraries (optimize optimize_lib ${PCL_LIBRARIES})
