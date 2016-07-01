#ifndef GETSQDIST_H
#define GETSQDIST_H

#include "grid.h"

using namespace std;

typedef boost::shared_ptr< pcl::VoxelGrid<pcl::InterestPoint> > VoxelGridPtr;
typedef boost::shared_ptr<grid> gridPtr;

//create grid with voxel values=confidence values of point cloud
gridPtr createGrid(pcl::PointCloud<pcl::InterestPoint>::Ptr grid_cloud, VoxelGridPtr vox, int res_factor);

//create binary volume grid from confidence grid
gridPtr getBinaryVolume(gridPtr grid_cloud);

//copy grid
gridPtr copyGrid(gridPtr in);

//add 2 grids voxel by voxel
gridPtr addGrids(gridPtr in1, gridPtr in2);

//extract perimeter of binary volume
gridPtr fastPerim(gridPtr volume_grid);

//implementation of Saito and Toriwaki distance field
gridPtr getsqdist(gridPtr volume_grid);

//get square root of grid
gridPtr getsqrt(gridPtr g);

#endif
