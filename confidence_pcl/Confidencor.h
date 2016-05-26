#ifndef CONFIDENCOR_H
#define CONFIDENCOR_H

#include <stdlib.h>

#include "full_confidence.h"

using namespace std;

typedef unsigned char byte;

class Confidencor
{
public:
    Confidencor(){}
    virtual ~Confidencor(){}
    virtual void conf_assigner(pcl::PointCloud<pcl::InterestPoint>::Ptr pc) = 0;

};

#endif
