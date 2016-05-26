#ifndef BINVOXTOPCL_H
#define BINVOXTOPCL_H

#endif // BINVOXTOPCL_H
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

int read_binvox(binvox* vox, string filespec);

pcl::PointCloud<pcl::PointXYZ>::Ptr binvoxToPCL(string filespec);
