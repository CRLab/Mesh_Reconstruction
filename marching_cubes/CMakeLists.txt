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
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../binvoxToPCL)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../confidence_pcl)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../voxelize)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../grid)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../distance_fields)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../narrow_band)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../optimize)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_library(mcubes_lib SHARED hasher.h hasher.cpp mcubes.h mcubes.cpp)

target_link_libraries (mcubes_lib grid_lib ${PCL_LIBRARIES})
add_executable (mcubes_test main.cpp)

target_link_libraries (mcubes_test mcubes_lib optimize_lib narrowBand_lib dfields_lib voxelize_lib assignConfidence_lib binvoxToPcl_lib)

install(TARGETS mcubes_lib 
        RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)