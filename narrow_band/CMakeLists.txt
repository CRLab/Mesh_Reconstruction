cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_BUILD_TYPE Debug)

project(narrowBand)

if(NOT PCL_INCLUDE_DIRS)
    find_package(PCL 1.2 REQUIRED)
endif(NOT PCL_INCLUDE_DIRS)

include_directories(${PCL_INCLUDE_DIRS})
include_directories(/home/adamjri/mesh_reconstruction/binvoxToPCL)
include_directories(/home/adamjri/mesh_reconstruction/confidence_pcl)
include_directories(/home/adamjri/mesh_reconstruction/voxelize)
include_directories(/home/adamjri/mesh_reconstruction/grid)
include_directories(/home/adamjri/mesh_reconstruction/distance_fields)

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
add_library(narrowBand_lib SHARED narrowBand.h narrowBand.cpp)

target_link_libraries (narrowBand_lib dfields_lib grid_lib ${PCL_LIBRARIES})
add_executable (narrowBand main.cpp)

target_link_libraries (narrowBand narrowBand_lib voxelize_lib assignConfidence_lib binvoxToPcl_lib)
