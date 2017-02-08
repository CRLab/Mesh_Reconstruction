cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(CMAKE_BUILD_TYPE Debug)

project(mesh_reconstruction)

find_package(PCL 1.2 REQUIRED)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(Eigen3 REQUIRED)
find_package(CUDA)
include_directories(EIGEN3_INCLUDE_DIR)

#find_package(Gperftools COMPONENTS profiler)

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

include_directories(${PCL_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/binvoxToPCL)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/confidence_pcl)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/voxelize)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/grid)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/distance_fields)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/narrow_band)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/optimize)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/optimize_conf)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/marching_cubes)

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

if(CUDA_FOUND)
  cuda_add_executable (mesh_reconstruction mesh_reconstruction.cpp)
  cuda_add_executable (m_cubes m_cubes.cpp)
else()
  add_executable (mesh_reconstruction mesh_reconstruction.cpp)
  add_executable (m_cubes m_cubes.cpp)
endif(CUDA_FOUND)

target_link_libraries (m_cubes mcubes_lib grid_lib dfields_lib voxelize_lib assignConfidence_lib ${PCL_LIBRARIES})
target_link_libraries (mesh_reconstruction mcubes_lib optimize_lib narrowBand_lib grid_lib dfields_lib voxelize_lib assignConfidence_lib binvoxToPcl_lib ${PCL_LIBRARIES} ${Boost_LIBRARIES})
#target_link_libraries (mesh_reconstruction mcubes_lib optimize_conf_lib optimize_lib narrowBand_lib grid_lib dfields_lib voxelize_lib assignConfidence_lib binvoxToPcl_lib ${PCL_LIBRARIES})

install(TARGETS m_cubes mesh_reconstruction 
        RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)