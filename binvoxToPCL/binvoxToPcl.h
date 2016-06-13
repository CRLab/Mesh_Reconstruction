#ifndef BINVOXTOPCL_H
#define BINVOXTOPCL_H

#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/surface/mls.h>

#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

typedef unsigned char byte;



static int version;
static int depth, height, width;
static int size;
static byte *voxels = 0;
static float tx, ty, tz;
static float scale;

struct binvox{
    int version;
    int depth, height, width;
    int size;
    byte *voxels;
    float tx, ty, tz;
    float scale;

};

int get_index(binvox vox, int x, int y, int z);

//function returns density of point cloud i.e. #points/volume of bounding box
float getDensity(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud){
    //get min/max points of bounding box
    float min_x=10000;
    float max_x=-10000;
    float min_y=10000;
    float max_y=-10000;
    float min_z=10000;
    float max_z=-10000;
    for(int i=0; i<cloud->points.size(); i++){
        if(cloud->points[i].x<min_x) min_x=cloud->points[i].x;
        if(cloud->points[i].x>max_x) max_x=cloud->points[i].x;
        if(cloud->points[i].y<min_y) min_y=cloud->points[i].y;
        if(cloud->points[i].y>max_y) max_y=cloud->points[i].y;
        if(cloud->points[i].z<min_z) min_z=cloud->points[i].z;
        if(cloud->points[i].z>max_z) max_z=cloud->points[i].z;
    }
    pcl::PointXYZ dims;
    dims.x=max_x-min_x; dims.y=max_y-min_y; dims.z = max_z-min_z;

    float volume=dims.x*dims.y*dims.z;
    return ((float)cloud->points.size())/volume;
}

int read_binvox(binvox* vox, string filespec);

pcl::PointCloud<pcl::PointXYZ>::Ptr binvoxToPCL(string filespec, int res_factor=1);

//automatically adjust point cloud resolutions based on density ratio:
pcl::PointCloud<pcl::PointXYZ>::Ptr binvoxToPCL_autores(string filespec, pcl::PointCloud<pcl::PointXYZ>::Ptr otherCloud){
    float otherDensity = getDensity(otherCloud);
    pcl::PointCloud<pcl::PointXYZ>::Ptr initCloud = binvoxToPCL(filespec);
    float density = getDensity(initCloud);
    if(density>otherDensity) return initCloud;
    else{
        int res_factor = otherDensity/density;
        return binvoxToPCL(filespec, res_factor);
    }
}

#endif
