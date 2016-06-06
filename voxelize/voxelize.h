#ifndef VOXELIZE_H
#define VOXELIZE_H

#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/surface/mls.h>
#include <boost/thread/thread.hpp>
#include <pcl/common/common_headers.h>
#include <pcl/features/normal_3d.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include <pcl/conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <pcl/PCLPointCloud2.h>
#include <pcl/octree/octree.h>

#include <iostream>

using namespace std;

typedef boost::shared_ptr< pcl::VoxelGrid<pcl::InterestPoint> > VoxelGridPtr;

struct voxelized_data{
    float resolution;
    pcl::PointCloud<pcl::InterestPoint>::Ptr input_cloud;
    pcl::PointCloud<pcl::InterestPoint>::Ptr filtered_cloud;
    VoxelGridPtr grid_data;
    pcl::octree::OctreePointCloudSearch<pcl::InterestPoint>::Ptr octree;
};
typedef boost::shared_ptr< voxelized_data> voxelized_dataPtr;

//create voxel grid data
VoxelGridPtr voxelize(pcl::PointCloud<pcl::InterestPoint>::Ptr input, pcl::PointCloud<pcl::InterestPoint>::Ptr output, float resolution);

//primes octree with cloud data for neighbor searches
pcl::octree::OctreePointCloudSearch<pcl::InterestPoint>::Ptr primeOctree(pcl::PointCloud<pcl::InterestPoint>::Ptr cloud, float resolution);

//returns list of indexes for points in the input cloud i.e. cloud->points[index]
vector<int> voxelNeighbors(pcl::octree::OctreePointCloudSearch<pcl::InterestPoint>::Ptr octree, pcl::InterestPoint searchPoint);

//go through voxelized data and modify strength of centroid to be average of voxel neighbors
//input:: voxel grid, filtered point cloud, original cloud, primed octree
void modifyStrengths(VoxelGridPtr grid, pcl::PointCloud<pcl::InterestPoint>::Ptr filtered_cloud,
                     pcl::PointCloud<pcl::InterestPoint>::Ptr cloud, pcl::octree::OctreePointCloudSearch<pcl::InterestPoint>::Ptr octree);


//function for creating voxelized point cloud
voxelized_dataPtr voxelizeData(pcl::PointCloud<pcl::InterestPoint>::Ptr cloud, float resolution);

//show point clouds in visualizer
void visualizeData(voxelized_data* data);

#endif
