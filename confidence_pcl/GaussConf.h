#ifndef GAUSSCONF_H
#define GAUSSCONF_H

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

class GaussConf: public Confidencor{
public:
    GaussConf();
    virtual void conf_assigner(pcl::PointCloud<pcl::InterestPoint>::Ptr pc);
};

#endif
