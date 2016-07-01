#ifndef GAUSSCONF_H
#define GAUSSCONF_H

#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>

#include <stdlib.h>

#include "Confidencor.h"

using namespace std;

typedef unsigned char byte;

class GaussConf: public Confidencor{
private:
    float variance;
public:
    GaussConf(float var);
    virtual void conf_assigner(pcl::PointCloud<pcl::InterestPoint>::Ptr pc);
};

#endif
