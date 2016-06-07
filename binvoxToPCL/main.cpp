#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/surface/mls.h>
#include <boost/thread/thread.hpp>
#include <pcl/common/common_headers.h>
#include <pcl/features/normal_3d.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include "binvoxToPcl.h"

using namespace std;

typedef unsigned char byte;

int main(int argc, char **argv){
    if(argc != 2){
        cout <<"Usage: convertBinvox <binvox filename>" << endl << endl;
        exit(1);
    }

    binvox vox;
    if(!read_binvox(&vox, argv[1])){
        cout << "Error reading [" << argv[1] << "]" << endl << endl;
        exit(1);
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = binvoxToPCL(argv[1]);

    stringstream ss;
    ss << "cloud.pcd";
    pcl::PCDWriter writer;
    writer.write<pcl::PointXYZ> (ss.str (), *cloud, false);

    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
     viewer->setBackgroundColor (0, 0, 0);
     pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color(cloud, 0, 255, 0);
     viewer->addPointCloud<pcl::PointXYZ> (cloud, single_color, "sample cloud");
     viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "sample cloud");
     viewer->addCoordinateSystem (1.0);
     viewer->initCameraParameters ();
     while (!viewer->wasStopped ())
     {
         viewer->spinOnce (100);
         boost::this_thread::sleep (boost::posix_time::microseconds (100000));
     }
}
