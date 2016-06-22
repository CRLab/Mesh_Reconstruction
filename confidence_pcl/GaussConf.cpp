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
#include "GaussConf.h"

using namespace std;

typedef unsigned char byte;

GaussConf::GaussConf(){}

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
            //get confidence and distance
            float conf = observed->points[pointIdxNKNSearch[0]].strength;
            float dist = (float) sqrt((double)pointNKNSquaredDistance[0]);
            //compute new confidence from gaussian
            float var = 2.0;   //<-------------------------------------------------need a way to compute reasonable variance
            float exp = -dist*dist/(var*var);
            float gauss = (1.0f/((float)sqrt(var*2.0)*3.1415927f))*pow(2.71828f,exp);
            pc->points[i].strength = gauss*conf;
        }
    }
}
