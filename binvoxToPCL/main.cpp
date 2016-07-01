#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include "binvoxToPcl.h"

using namespace std;

typedef unsigned char byte;

int main(int argc, char **argv){
    //convert binvox to pcl
    if(argc != 3){
        cout <<"Usage: Optimize Mesh with Confidences <binvox completed filename> <pcd observed filename>" << endl << endl;
        exit(1);
    }

    //read in data from pcd file
    pcl::PointCloud<pcl::PointXYZ>::Ptr observeCloud (new pcl::PointCloud<pcl::PointXYZ>());

    if (pcl::io::loadPCDFile<pcl::PointXYZ> (argv[2], *observeCloud) == -1) //* load the file
    {
    PCL_ERROR ("Couldn't read pcd file \n");
    exit(1);
    }

    //pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = binvoxToPCL(argv[1]);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = binvoxToPCL_autores(argv[1], observeCloud);

    stringstream ss;
    ss << "cloud.pcd";
    pcl::PCDWriter writer;
    writer.write<pcl::PointXYZ> (ss.str (), *cloud, false);

    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->setBackgroundColor (0, 0, 0);
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color1(observeCloud, 255, 0, 0);
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color2(cloud, 0, 255, 0);
    viewer->addPointCloud<pcl::PointXYZ> (observeCloud, single_color1, "sample cloud1");
    viewer->addPointCloud<pcl::PointXYZ> (cloud, single_color2, "sample cloud2");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "sample cloud1");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "sample cloud2");
    viewer->addCoordinateSystem (1.0);
    viewer->initCameraParameters ();
    while (!viewer->wasStopped ())
    {
     viewer->spinOnce (100);
     boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }
}
