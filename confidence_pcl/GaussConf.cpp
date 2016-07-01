#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>

#include <stdlib.h>

#include "Confidencor.h"
#include "GaussConf.h"

using namespace std;

typedef unsigned char byte;

GaussConf::GaussConf(float var):variance(var){}

void GaussConf::conf_assigner(pcl::PointCloud<pcl::InterestPoint>::Ptr pc){
    //create a copy of the observed point cloud
    pcl::PointCloud<pcl::InterestPoint>::Ptr observed (new pcl::PointCloud<pcl::InterestPoint>());
    for(int i=0; i<pc->points.size(); i++){
        if(pc->points[i].strength==1.0){
            observed->push_back(pc->points[i]);
        }
    }

    //create kd tree
    pcl::KdTreeFLANN<pcl::InterestPoint> kdtree;
    kdtree.setInputCloud(observed);

    for(int i=0; i<pc->points.size(); i++){
        if(pc->points[i].strength<0.0){
            //perform nearest neighbor search
            vector<int> pointIdxNKNSearch(1);
            vector<float> pointNKNSquaredDistance(1);
            kdtree.nearestKSearch(pc->points[i], 1, pointIdxNKNSearch, pointNKNSquaredDistance);
            //get distance
            float dist = (float) sqrt((double)pointNKNSquaredDistance[0]);
            //compute new confidence from gaussian
            float exp = -dist*dist/(2*variance); //<-------------------------------------------------need a way to compute reasonable variance
            float gauss = pow(2.71828f,exp);
            pc->points[i].strength = gauss;
        }
    }
}
