#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>

#include <stdlib.h>

#include "full_confidence.h"

pcl::PointCloud<pcl::InterestPoint>::Ptr full_confidence(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud){
    pcl::PointCloud<pcl::InterestPoint>::Ptr pc (new pcl::PointCloud<pcl::InterestPoint> ());
    for(int i=0; i<cloud->points.size(); i++){
        pcl::PointXYZ pnt = cloud->points[i];
        pcl::InterestPoint p;
        p.x=pnt.x; p.y=pnt.y; p.z=pnt.z;
        p.strength=1.0;
        //p.strength=0.0;
        pc->push_back(p);
    }

    return pc;
}
