#ifndef CONSTCONF_H
#define CONSTCONF_H

#include <pcl/point_types.h>

#include <stdlib.h>

#include "full_confidence.h"
#include "Confidencor.h"

using namespace std;

typedef unsigned char byte;

class ConstConf: public Confidencor{
private:
    float confidence;
public:
    ConstConf(float conf);
    virtual void conf_assigner(pcl::PointCloud<pcl::InterestPoint>::Ptr pc);
};

#endif
