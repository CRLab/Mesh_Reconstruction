cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_BUILD_TYPE Debug)

project(optimize)

if(NOT PCL_INCLUDE_DIRS)
    find_package(PCL 1.2 REQUIRED)
endif(NOT PCL_INCLUDE_DIRS)

find_package(Eigen3 REQUIRED)
include_directories(EIGEN3_INCLUDE_DIR)

include_directories(${PCL_INCLUDE_DIRS})
include_directories(/home/adamjri/mesh_reconstruction/binvoxToPCL)
include_directories(/home/adamjri/mesh_reconstruction/confidence_pcl)
include_directories(/home/adamjri/mesh_reconstruction/voxelize)
include_directories(/home/adamjri/mesh_reconstruction/narrow_band)
include_directories(/home/adamjri/mesh_reconstruction/grid)
include_directories(/home/adamjri/mesh_reconstruction/distance_fields)
include_directories(/home/adamjri/mesh_reconstruction/marching_cubes)

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
add_library(optimize_lib SHARED quadprog.h quadprog.cpp primeqp.h primeqp.cpp)

target_link_libraries (optimize_lib narrowBand_lib dfields_lib grid_lib ${PCL_LIBRARIES} voxelize_lib assignConfidence_lib binvoxToPcl_lib)
add_executable (optimize main.cpp)

target_link_libraries (optimize optimize_lib mcubes_lib voxelize_lib assignConfidence_lib binvoxToPcl_lib )
