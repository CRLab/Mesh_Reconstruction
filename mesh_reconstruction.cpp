
#include "binvoxToPcl.h"
#include "assign_confidence.h"

#include "ConstConf.h"

#include "voxelize.h"

#include "narrowBand.h"

#include "quadprog.h"

#include "mcubes.h"

#include <iostream>

using namespace std;

int main(int argc, char **argv){
    if(argc!=5 || (strcmp(argv[4], "0")!=0 && strcmp(argv[4], "1")!=0)){
        cout<<"Usage: <binvox file> <pcd file> <output file> <0 or 1 to toggle feature handling>"<<endl;
        return 1;
    }
    bool USING_FEATURES;
    if(strcmp(argv[4], "1")==0){
        cout<<"Using Feature Detection"<<endl;
        USING_FEATURES=true;
    }
    else{
        cout<<"Not Using Features"<<endl;
        USING_FEATURES=false;
    }

    //thresholds for feature detection
    //threshold values must be in (0,1)
    const float FEATURE_THRESHOLD = 0.75; //<----------------increasing raises sensitivity
    const float CORNER_THRESHOLD = 0.8; //<-----------------decreasing raises sensitivity
    //**************************************************************************************

    //perform algorithm

    /* Get full path to pcd, binvox, output files */
    string pcd_path (argv[2]);
    string binvox_path (argv[1]);
    string output_path (argv[3]);

    /* Read in data from pcd and binvox files */
    pcl::PointCloud<pcl::PointXYZ>::Ptr observeCloud (new pcl::PointCloud<pcl::PointXYZ>());
    if (pcl::io::loadPCDFile<pcl::PointXYZ> (pcd_path.c_str(), *observeCloud) == -1){ //* load the file
        PCL_ERROR ("Couldn't read pcd file \n");
        exit(1);
    }
    int res = getResolutionFactor(binvox_path.c_str(), observeCloud);
    pcl::PointCloud<pcl::PointXYZ>::Ptr predictCloud = binvoxToPCL(binvox_path.c_str(), res);


    /* Combine into pcl_conf with confidences */
    Confidencor *confidence_assigner = new ConstConf(1); //<--- change confidencor function here
    //assign full confidence to observeCloud
    pcl::PointCloud<pcl::InterestPoint>::Ptr confPCL=full_confidence(observeCloud);
    //asign confidence to everything
    assign_confidence(confPCL, predictCloud, confidence_assigner);


    /* Voxelize the data */
    voxelized_dataPtr data = voxelizeData(confPCL);
    //create grids
    gridPtr grid_cloud = createGrid(data->filtered_cloud, data->grid_data, res);
    gridPtr volume = getBinaryVolume(grid_cloud);


    /* Perform feature detection */
    //detect features
    vector<int> surface = getSurface(volume);
    gridPtr surfaceMap = getIndexMap(volume, surface);
    vector<Eigen::Vector3f> normals = getSurfaceNormals(volume, surface);
    gridPtr featureMap = getFeatureMap(volume, surfaceMap, normals, FEATURE_THRESHOLD);

    /* Perform smoothing */
    gridPtr F = optimize(volume, featureMap, USING_FEATURES);

    /* Extract mesh and write to file */
    mcubes(F, surfaceMap, normals, 0.0, FEATURE_THRESHOLD, CORNER_THRESHOLD, output_path.c_str(), USING_FEATURES);

    return 1;
}
