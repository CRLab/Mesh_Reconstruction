cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_BUILD_TYPE Debug)

project(voxelize)

if(NOT PCL_INCLUDE_DIRS)
    find_package(PCL 1.2 REQUIRED)
endif(NOT PCL_INCLUDE_DIRS)

include_directories(${PCL_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../binvoxToPCL)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../confidence_pcl)

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
add_library(voxelize_lib SHARED voxelize.h voxelize.cpp)

target_link_libraries (voxelize_lib ${PCL_LIBRARIES})
add_executable (voxelize main.cpp)

target_link_libraries (voxelize voxelize_lib binvoxToPcl_lib assignConfidence_lib)
