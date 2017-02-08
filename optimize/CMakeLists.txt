cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_BUILD_TYPE Debug)

project(optimize)

if(NOT PCL_INCLUDE_DIRS)
    find_package(PCL 1.2 REQUIRED)
endif(NOT PCL_INCLUDE_DIRS)

find_package(CUDA)
find_package(Eigen3 REQUIRED)
include_directories(EIGEN3_INCLUDE_DIR)

include_directories(${PCL_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../binvoxToPCL)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../confidence_pcl)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../voxelize)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../narrow_band)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../grid)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../distance_fields)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../marching_cubes)

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
if(CUDA_FOUND)
CUDA_ADD_LIBRARY(optimize_lib SHARED quadprog.h quadprog.cu primeqp.h primeqp.cpp)
else()
add_library(optimize_lib SHARED quadprog.h quadprog.cpp primeqp.h primeqp.cpp)
endif(CUDA_FOUND)

target_link_libraries (optimize_lib narrowBand_lib dfields_lib grid_lib ${PCL_LIBRARIES} voxelize_lib assignConfidence_lib binvoxToPcl_lib)
if(CUDA_FOUND)
cuda_add_executable (optimize main.cpp)
else()
add_executable(optimize main.cpp)
endif(CUDA_FOUND)

target_link_libraries (optimize optimize_lib mcubes_lib voxelize_lib assignConfidence_lib binvoxToPcl_lib )

install(TARGETS optimize_lib 
        RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)
