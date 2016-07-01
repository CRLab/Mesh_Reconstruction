
#include "voxelize.h"

#include "binvoxToPcl.h"
#include "assign_confidence.h"

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

    //voxelize the data
    voxelized_dataPtr data = voxelizeData(confPCL); //<--test different resolutions

    //visualize
    visualizeData(data);

    return 1;
}
