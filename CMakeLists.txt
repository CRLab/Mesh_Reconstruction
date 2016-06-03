cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(mesh_reconstruction)

find_package(PCL 1.2 REQUIRED)
#SET(CMAKE_CXX_FLAGS "-std=c++0x -Wno-deprecated")

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_subdirectory(binvoxToPCL)
add_subdirectory(confidence_pcl)
add_subdirectory(voxelize)
add_subdirectory(narrow_band)

add_executable (mesh_reconstruction main.cpp)
target_link_libraries (mesh_reconstruction narrowBand_lib voxelize_lib assignConfidence_lib binvoxToPcl_lib ${PCL_LIBRARIES})
