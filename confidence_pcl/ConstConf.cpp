#include <pcl/point_types.h>

#include <stdlib.h>

#include "full_confidence.h"
#include "Confidencor.h"
#include "ConstConf.h"

using namespace std;

typedef unsigned char byte;

ConstConf::ConstConf(float conf):confidence(conf){}

void ConstConf::conf_assigner(pcl::PointCloud<pcl::InterestPoint>::Ptr pc){
    for(int i=0; i<pc->points.size(); i++){
        if(pc->points[i].strength<0){
            pc->points[i].strength = confidence;
        }
    }
}
