#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include "binvoxToPcl.h"
#include "assign_confidence.h"

#include "ConstConf.h"

#include "voxelize.h"

#include "narrowBand.h"


using namespace std;

typedef unsigned char byte;

typedef boost::shared_ptr<grid> gridPtr;

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

    int res = getResolutionFactor(argv[1], observeCloud);
    pcl::PointCloud<pcl::PointXYZ>::Ptr predictCloud = binvoxToPCL(argv[1], res);

    //combine into pcl_conf with confidences
    Confidencor *confidence_assigner = new ConstConf(0); //<--- change confidencor function here

    //assign full confidence to observeCloud
    pcl::PointCloud<pcl::InterestPoint>::Ptr confPCL=full_confidence(observeCloud);

    //asign confidence to everything
    assign_confidence(confPCL, predictCloud, confidence_assigner);

    //voxelize the data
    voxelized_dataPtr data = voxelizeData(confPCL); //<--test different resolutions

    //create grids
    gridPtr grid_cloud = createGrid(data->filtered_cloud, data->grid_data, res);
    gridPtr volume = getBinaryVolume(grid_cloud);
    gridPtr margin = getsqrt(getsqdist(volume));

    //create bands
    bandsPtr bnds = createBands(margin, 1.0);




    //create point clouds from bands for visualization
    pcl::PointCloud<pcl::InterestPoint>::Ptr pcl_band (new pcl::PointCloud<pcl::InterestPoint>());
    for(int i=0; i<bnds->band->dims[0]; i++){
        for(int j=0; j<bnds->band->dims[1]; j++){
            for(int k=0; k<bnds->band->dims[2]; k++){
                if((*(bnds->band))[i][j][k]==1.0){
                    pcl::InterestPoint pnt;
                    pnt.x=(float)i; pnt.y=(float)j; pnt.z=(float)k;
                    pnt.strength=0;
                    pcl_band->push_back(pnt);
                }
            }
        }
    }
    pcl::PointCloud<pcl::InterestPoint>::Ptr pcl_tightband (new pcl::PointCloud<pcl::InterestPoint>());
    for(int i=0; i<bnds->tight_band->dims[0]; i++){
        for(int j=0; j<bnds->tight_band->dims[1]; j++){
            for(int k=0; k<bnds->tight_band->dims[2]; k++){
                if((*(bnds->tight_band))[i][j][k]==1.0){
                    pcl::InterestPoint pnt;
                    pnt.x=(float)i+0.01; pnt.y=(float)j+0.01; pnt.z=(float)k+0.01;
                    pnt.strength=0;
                    pcl_tightband->push_back(pnt);
                }
            }
        }
    }

    //visualize
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->setBackgroundColor (0, 0, 0);
    //pcl::visualization::PointCloudColorHandlerCustom<pcl::InterestPoint> single_color(data->filtered_cloud, 0, 255, 0);
    pcl::visualization::PointCloudColorHandlerCustom<pcl::InterestPoint> single_color_tightband(pcl_tightband, 0, 0, 255);
    pcl::visualization::PointCloudColorHandlerCustom<pcl::InterestPoint> single_color_band(pcl_band, 255, 0, 0);
    //viewer->addPointCloud<pcl::InterestPoint> (data->filtered_cloud, single_color, "cloud");
    viewer->addPointCloud<pcl::InterestPoint> (pcl_tightband, single_color_tightband, "tight band");
    viewer->addPointCloud<pcl::InterestPoint> (pcl_band, single_color_band, "band");
    //viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "cloud");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "tight band");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "band");
    viewer->addCoordinateSystem (1.0);
    viewer->initCameraParameters ();
    while (!viewer->wasStopped ())
    {
        viewer->spinOnce (100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }


    return 1;
}
