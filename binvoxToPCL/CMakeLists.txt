cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(convertBinvox)

if(NOT PCL_INCLUDE_DIRS)
    find_package(PCL 1.2 REQUIRED)
endif(NOT PCL_INCLUDE_DIRS)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
add_library(binvoxToPcl_lib SHARED binvoxToPcl.h binvoxToPcl.cpp)

target_link_libraries (binvoxToPcl_lib ${PCL_LIBRARIES})
add_executable (convertBinvox main.cpp)
target_link_libraries (convertBinvox binvoxToPcl_lib)
