
#include "assign_confidence.h"

#include "ConstConf.h"

#include "voxelize.h"

#include "dfields.h"

#include "mcubes.h"

#include <iostream>
#include <string>

using namespace std;


int main(int argc, char **argv){
    //************************************************************
    //feature values not used
    const float FEATURE_THRESHOLD = 0.75;
    const float CORNER_THRESHOLD = 0.8;
    //************************************************************

    if((argc!=3) and (argc != 4)){
        cerr<<"Usage: <pcd filename> <output filename (.ply)> [optional]leaf_size"<<endl;
        return 1;
    }

    cout<<"Running Standard Marching Cubes"<<endl;

    string PCD_PATH (argv[1]);
    string OUT_PATH (argv[2]);
    bool resolution_in_args = false;
    float resolution;
    if (argc == 4)
      {
	resolution_in_args = true;
      }
    
    
    /* Read in data from pcd and binvox files */
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>());
    if (pcl::io::loadPCDFile<pcl::PointXYZ> (PCD_PATH.c_str(), *cloud) == -1){ //* load the file
        PCL_ERROR ("Couldn't read pcd file \n");
        exit(1);
    }
    //assign full confidence to observeCloud
    pcl::PointCloud<pcl::InterestPoint>::Ptr confPCL=full_confidence(cloud);

    /* Voxelize the data */
    //create dummy object
    pcl::PointCloud<pcl::InterestPoint>::Ptr dummy (new pcl::PointCloud<pcl::InterestPoint>());
    if (resolution_in_args)
      {
	resolution = atof(argv[3]);
    }else{
      resolution = getResolution(confPCL);
    }
    cout<<"Leaf size is "<<resolution<<endl;
    voxelized_dataPtr data (new voxelized_data());
    data->resolution=2*resolution;
    data->input_cloud=confPCL;
    data->filtered_cloud = dummy;
    data->grid_data = voxelize(data->input_cloud, data->filtered_cloud, data->resolution);
    data->octree = primeOctree(data->input_cloud, data->resolution);
    modifyStrengths(data->filtered_cloud, data->input_cloud, data->octree);

    //create grids
    gridPtr grid_cloud (new grid(data->filtered_cloud, data->grid_data));
    gridPtr volume = getBinaryVolume(grid_cloud);

    /* Perform feature detection */
    //detect features
    vector<int> surface = getSurface(volume);
    gridPtr surfaceMap = getIndexMap(volume, surface);
    vector<Eigen::Vector3f> normals = getSurfaceNormals(volume, surface);

    /* Extract mesh and write to file */
    mcubes(volume, surfaceMap, normals, 0.5, FEATURE_THRESHOLD, CORNER_THRESHOLD, OUT_PATH.c_str(), false);
}
