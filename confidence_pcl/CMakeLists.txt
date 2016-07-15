cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(assignConfidence)

if(NOT PCL_INCLUDE_DIRS)
    find_package(PCL 1.2 REQUIRED)
endif(NOT PCL_INCLUDE_DIRS)

include_directories(${PCL_INCLUDE_DIRS})
include_directories(/home/adamjri/mesh_reconstruction/binvoxToPCL)

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
add_library(assignConfidence_lib SHARED full_confidence.h full_confidence.cpp assign_confidence.h assign_confidence.cpp Confidencor.h ConstConf.h ConstConf.cpp GaussConf.h GaussConf.cpp)

target_link_libraries (assignConfidence_lib ${PCL_LIBRARIES})
add_executable (assignConfidence main.cpp)

target_link_libraries (assignConfidence assignConfidence_lib binvoxToPcl_lib)
