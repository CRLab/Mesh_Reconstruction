#ifndef FULL_CONFIDENCE_H
#define FULL_CONFIDENCE_H

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

pcl::PointCloud<pcl::InterestPoint>::Ptr full_confidence(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);

#endif
