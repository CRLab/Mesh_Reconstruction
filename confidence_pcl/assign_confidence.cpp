#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/surface/mls.h>

#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "full_confidence.h"
#include "Confidencor.h"

using namespace std;

typedef unsigned char byte;

pcl::PointCloud<pcl::InterestPoint>::Ptr assign_confidence(pcl::PointCloud<pcl::InterestPoint>::Ptr seed, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, Confidencor* f){
    for(int i=0; i<cloud->points.size(); i++){
        pcl::InterestPoint pnt;
        pnt.x=cloud->points[i].x; pnt.y=cloud->points[i].y; pnt.z=cloud->points[i].z;
        pnt.strength=-1;
        seed->push_back(pnt);
    }
    f->conf_assigner(seed);

    return seed;
}
