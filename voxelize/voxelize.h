#ifndef VOXELIZE_H
#define VOXELIZE_H

#include <pcl/kdtree/kdtree_flann.h>
#include <boost/thread/thread.hpp>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/octree/octree.h>

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

float getResolution(pcl::PointCloud<pcl::InterestPoint>::Ptr cloud);

//create voxel grid data
VoxelGridPtr voxelize(pcl::PointCloud<pcl::InterestPoint>::Ptr input, pcl::PointCloud<pcl::InterestPoint>::Ptr output, float resolution);

//primes octree with cloud data for neighbor searches
pcl::octree::OctreePointCloudSearch<pcl::InterestPoint>::Ptr primeOctree(pcl::PointCloud<pcl::InterestPoint>::Ptr cloud, float resolution);

//returns list of indexes for points in the input cloud i.e. cloud->points[index]
vector<int> voxelNeighbors(pcl::octree::OctreePointCloudSearch<pcl::InterestPoint>::Ptr octree, pcl::InterestPoint searchPoint);

//go through voxelized data and modify strength of centroid to be average of voxel neighbors
//input:: voxel grid, filtered point cloud, original cloud, primed octree
void modifyStrengths(pcl::PointCloud<pcl::InterestPoint>::Ptr filtered_cloud,
                     pcl::PointCloud<pcl::InterestPoint>::Ptr cloud, pcl::octree::OctreePointCloudSearch<pcl::InterestPoint>::Ptr octree);


//function for creating voxelized point cloud
voxelized_dataPtr voxelizeData(pcl::PointCloud<pcl::InterestPoint>::Ptr cloud);

//show point clouds in visualizer
void visualizeData(voxelized_dataPtr data);

#endif
