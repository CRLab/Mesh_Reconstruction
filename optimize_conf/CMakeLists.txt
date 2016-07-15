cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_BUILD_TYPE Debug)
project(optimize_conf)

find_package(PCL 1.2 REQUIRED)
#SET(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated")

include_directories(${PCL_INCLUDE_DIRS})
include_directories(/home/adamjri/mesh_reconstruction/binvoxToPCL)
include_directories(/home/adamjri/mesh_reconstruction/confidence_pcl)
include_directories(/home/adamjri/mesh_reconstruction/voxelize)
include_directories(/home/adamjri/mesh_reconstruction/grid)
include_directories(/home/adamjri/mesh_reconstruction/distance_fields)
include_directories(/home/adamjri/mesh_reconstruction/narrow_band)


link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
add_library(optimize_conf_lib SHARED qp_conf.h qp_conf.cpp primeqp_conf.h primeqp_conf.cpp)

target_link_libraries (optimize_conf_lib narrowBand_lib dfields_lib grid_lib voxelize_lib assignConfidence_lib binvoxToPcl_lib ${PCL_LIBRARIES})
add_executable (optimize_conf main.cpp)

target_link_libraries (optimize_conf optimize_conf_lib)
