#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/surface/mls.h>
#include <boost/thread/thread.hpp>
#include <pcl/common/common_headers.h>
#include <pcl/features/normal_3d.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include <pcl/conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <pcl/octree/octree.h>

#include <iostream>

#include "voxelize.h"

using namespace std;

VoxelGridPtr voxelize(pcl::PointCloud<pcl::InterestPoint>::Ptr input, pcl::PointCloud<pcl::InterestPoint>::Ptr output, float resolution){

    cout << "PointCloud before filtering: " << input->points.size()
           << " data points"<<endl;

    //create voxelized data
    VoxelGridPtr sor (new pcl::VoxelGrid<pcl::InterestPoint>());

    sor->setInputCloud(input);
    sor->setLeafSize (resolution, resolution, resolution);
    sor->setSaveLeafLayout(true);
    sor->setDownsampleAllData(true);
    sor->filter(*output);

    cout << "PointCloud after filtering: " << output->points.size()<<" data points"<<endl;

    return sor;
}

//primes octree with cloud data for neighbor searches
pcl::octree::OctreePointCloudSearch<pcl::InterestPoint>::Ptr primeOctree(pcl::PointCloud<pcl::InterestPoint>::Ptr cloud, float resolution){
    pcl::octree::OctreePointCloudSearch<pcl::InterestPoint>::Ptr octree (new pcl::octree::OctreePointCloudSearch<pcl::InterestPoint>(resolution));
    octree->setInputCloud (cloud);
    octree->addPointsFromInputCloud ();
    return octree;
}

//returns list of indexes for points in the input cloud i.e. cloud->points[index]
vector<int> voxelNeighbors(pcl::octree::OctreePointCloudSearch<pcl::InterestPoint>::Ptr octree, pcl::InterestPoint searchPoint){
    vector<int> pointIdxVec;
    octree->voxelSearch(searchPoint, pointIdxVec);
    return pointIdxVec;

}

//go through voxelized data and modify strength of centroid to be average of voxel neighbors
//input:: voxel grid, filtered point cloud, original cloud, primed octree
void modifyStrengths(pcl::PointCloud<pcl::InterestPoint>::Ptr filtered_cloud,
                     pcl::PointCloud<pcl::InterestPoint>::Ptr cloud, pcl::octree::OctreePointCloudSearch<pcl::InterestPoint>::Ptr octree){
    for(int i=0; i<filtered_cloud->points.size(); i++){
        //get neighbors in voxel
        vector<int> pointIdxVec = voxelNeighbors(octree, filtered_cloud->points[i]);
        //take average of strengths
        float avg=0.0;
        for(int ind=0; ind<pointIdxVec.size(); ind++){
            avg+=cloud->points[pointIdxVec[ind]].strength;
        }
        avg=avg/(float)pointIdxVec.size();
        //replace centroid strength with average
        filtered_cloud->points[i].strength=avg;
    }//end of loop through voxel grid
}

voxelized_dataPtr voxelizeData(pcl::PointCloud<pcl::InterestPoint>::Ptr cloud, float resolution){
    //create dummy object
    pcl::PointCloud<pcl::InterestPoint>::Ptr dummy (new pcl::PointCloud<pcl::InterestPoint>());

    voxelized_dataPtr data (new voxelized_data());
    data->resolution=resolution;
    data->input_cloud=cloud;
    data->filtered_cloud = dummy;
    data->grid_data = voxelize(data->input_cloud, data->filtered_cloud, data->resolution);
    data->octree = primeOctree(data->input_cloud, data->resolution);
    modifyStrengths(data->filtered_cloud, data->input_cloud, data->octree);
    return data;
}

void visualizeData(voxelized_dataPtr data){
    //create rgb point cloud for visualization
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr confPclRGB (new pcl::PointCloud<pcl::PointXYZRGB>());
    for(int i=0; i<data->input_cloud->points.size(); i++){
        pcl::PointXYZRGB pnt;
        pnt.x=data->input_cloud->points[i].x;pnt.y=data->input_cloud->points[i].y;pnt.z=data->input_cloud->points[i].z;
        pnt.r=0; pnt.b=0;
        pnt.g=100+155*data->input_cloud->points[i].strength;
        confPclRGB->push_back(pnt);
    }
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr filtered_cloudRGB (new pcl::PointCloud<pcl::PointXYZRGB>());
    Eigen::Vector3i dims_min = data->grid_data->getMinBoxCoordinates();
    Eigen::Vector3i dims_max = data->grid_data->getMaxBoxCoordinates();
    for(int i=dims_min[0]; i<=dims_max[0]; i++){
        for(int j=dims_min[1]; j<=dims_max[1]; j++){
            for(int k=dims_min[2]; k<=dims_max[2]; k++){
                Eigen::Vector3i pnt;
                pnt[0]=i;pnt[1]=j;pnt[2]=k;
                int index = data->grid_data->getCentroidIndexAt(pnt);
                if(index!=-1){
                    pcl::PointXYZRGB p;
                    p.x=data->filtered_cloud->points[index].x;p.y=data->filtered_cloud->points[index].y;p.z=data->filtered_cloud->points[index].z;
                    p.g=0; p.b=0;
                    p.r=100+155*data->filtered_cloud->points[i].strength;
                    filtered_cloudRGB->push_back(p);
                }
            }
        }
    }

    //display in visualizor
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->setBackgroundColor (0, 0, 0);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(confPclRGB);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb2(filtered_cloudRGB);
    viewer->addPointCloud<pcl::PointXYZRGB> (confPclRGB, rgb, "cloud");
    viewer->addPointCloud<pcl::PointXYZRGB> (filtered_cloudRGB, rgb2, "filtered cloud");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "cloud");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "filtered cloud");
    viewer->addCoordinateSystem (1.0);
    viewer->initCameraParameters ();
    while (!viewer->wasStopped ())
    {
        viewer->spinOnce (100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }

}

