#ifndef GETSQDIST_H
#define GETSQDIST_H

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/filters/voxel_grid.h>

#include <iostream>

using namespace std;

typedef boost::shared_ptr< pcl::VoxelGrid<pcl::InterestPoint> > VoxelGridPtr;

struct grid{
    float*** voxels;
    Eigen::Vector3i t_;
    Eigen::Vector3i dims;
    ~grid() {
        for(int i=0; i<dims[0]; i++){
            for(int j=0; j<dims[1]; j++){
                delete [] voxels[i][j];
            }
            delete [] voxels[i];
        }
        delete [] voxels;
    }
};
typedef boost::shared_ptr<grid> gridPtr;

float*** allocGrid(Eigen::Vector3i dims);

//create grid with voxel values=confidence values of point cloud
gridPtr createGrid(pcl::PointCloud<pcl::InterestPoint>::Ptr grid_cloud, VoxelGridPtr vox);

//create binary volume grid from confidence grid
gridPtr getBinaryVolume(gridPtr grid_cloud);

//copy grid
gridPtr copyGrid(gridPtr in);

//get negative of grid
gridPtr getNegGrid(gridPtr in);

//add 2 grids voxel by voxel
gridPtr addGrids(gridPtr in1, gridPtr in2);

//extract perimeter of binary volume
gridPtr fastPerim(gridPtr volume_grid);

//implementation of Saito and Toriwaki distance field
gridPtr getsqdist(gridPtr volume_grid);

//get square root of grid
gridPtr getsqrt(gridPtr g);

#endif
