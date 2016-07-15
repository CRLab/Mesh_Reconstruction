cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_BUILD_TYPE Debug)
project(mcubes)

#SET(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated")
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.54.0 REQUIRED)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

if(NOT PCL_INCLUDE_DIRS)
    find_package(PCL 1.2 REQUIRED)
endif(NOT PCL_INCLUDE_DIRS)

include_directories(${PCL_INCLUDE_DIRS})
include_directories(/home/adamjri/mesh_reconstruction/binvoxToPCL)
include_directories(/home/adamjri/mesh_reconstruction/confidence_pcl)
include_directories(/home/adamjri/mesh_reconstruction/voxelize)
include_directories(/home/adamjri/mesh_reconstruction/grid)
include_directories(/home/adamjri/mesh_reconstruction/distance_fields)
include_directories(/home/adamjri/mesh_reconstruction/narrow_band)
include_directories(/home/adamjri/mesh_reconstruction/optimize)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_library(mcubes_lib SHARED hasher.h hasher.cpp mcubes.h mcubes.cpp)

target_link_libraries (mcubes_lib grid_lib ${PCL_LIBRARIES})
add_executable (mcubes main.cpp)

target_link_libraries (mcubes mcubes_lib optimize_lib narrowBand_lib dfields_lib voxelize_lib assignConfidence_lib binvoxToPcl_lib)
