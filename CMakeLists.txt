cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(mesh_reconstruction)

find_package(PCL 1.2 REQUIRED)

find_package(Eigen3 REQUIRED)
include_directories(EIGEN3_INCLUDE_DIR)

include_directories(${PCL_INCLUDE_DIRS})
include_directories(/home/adamjri/mesh_reconstruction/binvoxToPCL)
include_directories(/home/adamjri/mesh_reconstruction/confidence_pcl)
include_directories(/home/adamjri/mesh_reconstruction/voxelize)
include_directories(/home/adamjri/mesh_reconstruction/grid)
include_directories(/home/adamjri/mesh_reconstruction/distance_fields)
include_directories(/home/adamjri/mesh_reconstruction/narrow_band)
include_directories(/home/adamjri/mesh_reconstruction/optimize)
#include_directories(/home/adamjri/mesh_reconstruction/optimize_conf)
include_directories(/home/adamjri/mesh_reconstruction/marching_cubes)

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_subdirectory(binvoxToPCL)
add_subdirectory(confidence_pcl)
add_subdirectory(voxelize)
add_subdirectory(grid)
add_subdirectory(distance_fields)
add_subdirectory(narrow_band)
add_subdirectory(optimize)
#add_subdirectory(optimize_conf)
add_subdirectory(marching_cubes)

add_executable (mesh_reconstruction main.cpp)
target_link_libraries (mesh_reconstruction mcubes_lib optimize_lib narrowBand_lib grid_lib dfields_lib voxelize_lib assignConfidence_lib binvoxToPcl_lib ${PCL_LIBRARIES})
#target_link_libraries (mesh_reconstruction mcubes_lib optimize_conf_lib optimize_lib narrowBand_lib grid_lib dfields_lib voxelize_lib assignConfidence_lib binvoxToPcl_lib ${PCL_LIBRARIES})
