#ifndef FULL_CONFIDENCE_H
#define FULL_CONFIDENCE_H

#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>

#include <stdlib.h>

using namespace std;

typedef unsigned char byte;

pcl::PointCloud<pcl::InterestPoint>::Ptr full_confidence(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);

#endif
