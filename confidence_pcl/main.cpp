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
#include "assign_confidence.h"

//include different confidencors to test
#include "ConstConf.h"

using namespace std;

typedef unsigned char byte;

int main(int argc, char **argv){
    //convert binvox to pcl
    if(argc != 3){
        cout <<"Usage: assignConfidence <binvox filename> <pcd filename>" << endl << endl;
        exit(1);
    }

    //read in data from pcd file
    pcl::PointCloud<pcl::PointXYZ>::Ptr observeCloud (new pcl::PointCloud<pcl::PointXYZ>());

    if (pcl::io::loadPCDFile<pcl::PointXYZ> (argv[2], *observeCloud) == -1) //* load the file
    {
    PCL_ERROR ("Couldn't read pcd file \n");
    exit(1);
    }

    binvox vox;
    if(!read_binvox(&vox, argv[1])){
        cout << "Error reading [" << argv[1] << "]" << endl << endl;
        exit(1);
    }
    pcl::PointCloud<pcl::PointXYZ>::Ptr predictCloud = binvoxToPCL(argv[1]);






    //combine into pcl_conf with confidences
    Confidencor *confidence_assigner = new ConstConf(0); //<--- change confidencor function here

    //assign full confidence to observeCloud
    pcl::PointCloud<pcl::InterestPoint>::Ptr confPCL=full_confidence(observeCloud);

    //asign confidence to everything
    assign_confidence(confPCL, predictCloud, confidence_assigner);
    //cout<<"confPCL size: "<<confPCL->points.size()<<endl<<endl;


    //create rgb point cloud for visualization
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr confPclRGB (new pcl::PointCloud<pcl::PointXYZRGB>());
    for(int i=0; i<confPCL->points.size(); i++){
        pcl::PointXYZRGB pnt;
        pnt.x=confPCL->points[i].x;pnt.y=confPCL->points[i].y;pnt.z=confPCL->points[i].z;
        pnt.r=0; pnt.b=0;
        pnt.g=100+155*confPCL->points[i].strength;
        confPclRGB->push_back(pnt);
    }

    //display in visualizor
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->setBackgroundColor (0, 0, 0);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(confPclRGB);
    viewer->addPointCloud<pcl::PointXYZRGB> (confPclRGB, rgb, "sample cloud");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "sample cloud");
    viewer->addCoordinateSystem (1.0);
    viewer->initCameraParameters ();
    while (!viewer->wasStopped ())
    {
        viewer->spinOnce (100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }

    return 1;
}
