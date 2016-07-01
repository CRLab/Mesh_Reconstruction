
#include "binvoxToPcl.h"
#include "assign_confidence.h"

//include confidencor functions
#include "ConstConf.h"
#include "GaussConf.h"

#include "voxelize.h"

#include "getsqdist.h"
#include "narrowBand.h"

#include "qp_conf.h"

using namespace std;
using namespace conf;

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


    int res = getResolutionFactor(argv[1], observeCloud);
    pcl::PointCloud<pcl::PointXYZ>::Ptr predictCloud = binvoxToPCL(argv[1], res);

    //combine into pcl_conf with confidences
    Confidencor *confidence_assigner = new GaussConf(2.0); //<--- change confidencor function here

    //assign full confidence to observeCloud
    pcl::PointCloud<pcl::InterestPoint>::Ptr confPCL=full_confidence(observeCloud);

    //asign confidence to everything
    assign_confidence(confPCL, predictCloud, confidence_assigner);

    //voxelize the data
    voxelized_dataPtr data = voxelizeData(confPCL); //<--test different resolutions

    //create grids
    gridPtr grid_cloud = createGrid(data->filtered_cloud, data->grid_data, res);
    gridPtr volume = getBinaryVolume(grid_cloud);

    //get imbedding function
    gridPtr F = optimize(grid_cloud, volume);
    //gridPtr F = volume;

    //write to file
    ofstream myfile;
    myfile.open("../embed_funcs/test_conf.txt");
    for(int i=0; i<F->dims[0]; i++){
        for(int j=0; j<F->dims[1]; j++){
            for(int k=0; k<F->dims[2]; k++){
                myfile<<i<<","<<j<<","<<k<<","<<(*F)[i][j][k]<<endl;
            }
        }
    }
    myfile.close();


    //visualize
    visualizeGrid(F);

    return 1;
}
