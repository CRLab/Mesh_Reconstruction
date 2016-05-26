#ifndef ASSIGN_CONFIDENCE_H
#define ASSIGN_CONFIDENCE_H

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

pcl::PointCloud<pcl::InterestPoint>::Ptr assign_confidence(pcl::PointCloud<pcl::InterestPoint>::Ptr seed, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, Confidencor* f);

#endif
