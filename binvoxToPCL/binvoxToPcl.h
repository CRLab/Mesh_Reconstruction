#ifndef BINVOXTOPCL_H
#define BINVOXTOPCL_H

#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree_flann.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

typedef unsigned char byte;

struct binvox{
    int version;
    int depth, height, width;
    int size;
    byte *voxels;
    float tx, ty, tz;
    float scale;

};

int get_index(const binvox &vox, int x, int y, int z);

int read_binvox(binvox* vox, string filespec);

pcl::PointCloud<pcl::PointXYZ>::Ptr binvoxToPCL(string filespec, int res_factor=1);

//function returns density of point cloud i.e. #points/volume of bounding box
float getDensity(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);

int getResolutionFactor(string filespec, pcl::PointCloud<pcl::PointXYZ>::Ptr otherCloud);

//automatically adjust point cloud resolutions based on density ratio:
pcl::PointCloud<pcl::PointXYZ>::Ptr binvoxToPCL_autores(string filespec, pcl::PointCloud<pcl::PointXYZ>::Ptr otherCloud);

#endif
