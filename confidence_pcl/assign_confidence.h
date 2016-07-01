#ifndef ASSIGN_CONFIDENCE_H
#define ASSIGN_CONFIDENCE_H

#include <pcl/point_types.h>

#include <stdlib.h>

#include "full_confidence.h"
#include "Confidencor.h"

using namespace std;

typedef unsigned char byte;

pcl::PointCloud<pcl::InterestPoint>::Ptr assign_confidence(pcl::PointCloud<pcl::InterestPoint>::Ptr seed, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, Confidencor* f);

#endif
