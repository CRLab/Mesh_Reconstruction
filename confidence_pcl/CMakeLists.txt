cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(assignConfidence)

find_package(PCL 1.2 REQUIRED)
SET(CMAKE_CXX_FLAGS "-std=c++0x")

include_directories(${PCL_INCLUDE_DIRS} /home/adamjri/mesh_reconstruction/binvoxToPCL)

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
add_library(assignConfidence_lib SHARED full_confidence.h full_confidence.cpp assign_confidence.h assign_confidence.cpp Confidencor.h ConstConf.h ConstConf.cpp)

target_link_libraries (assignConfidence_lib  binvoxToPcl_lib ${PCL_LIBRARIES})
add_executable (assignConfidence main.cpp)

target_link_libraries (assignConfidence assignConfidence_lib ${PCL_LIBRARIES})
