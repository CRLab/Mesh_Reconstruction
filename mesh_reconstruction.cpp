
#include "binvoxToPcl.h"
#include "assign_confidence.h"

#include "ConstConf.h"

#include "voxelize.h"

#include "narrowBand.h"

#include "quadprog.h"

#include "mcubes.h"

#include <iostream>

#include <exception>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv){

    //thresholds for feature detection
    //threshold values must be in (0,1)
    float FEATURE_THRESHOLD = 0.75; //<----------------increasing raises sensitivity
    float CORNER_THRESHOLD = 0.8; //<-----------------decreasing raises sensitivity
    bool USING_FEATURES = false; //<------------------determines if feature detection is used
    bool USING_CUDA = false;     //<------------------determines if using GPU based algorithm
    //**************************************************************************************

    try {
        po::options_description desc("Allowed options");
        desc.add_options()
                ("help", "produce help message")
                ("feature-detection", po::bool_switch(&USING_FEATURES), "Toggle feature handling")
                ("cuda", po::bool_switch(&USING_CUDA), "Toggle CUDA option")
                ("feature-threshold", po::value<float>(), "Increasing raises feature sensitivity. Default: 0.75")
                ("corner-threshold", po::value<float>(), "Decreasing raises feature sensitivity. Default: 0.8")
                ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if(argc < 4) {
            cout<<"Usage: mesh_reconstruction <binvox file> <pcd file> <output file> <options>"<<endl;
            cout << desc << endl;
            return 1;
        }

        if (vm.count("help")) {
            cout<<"Usage: mesh_reconstruction <binvox file> <pcd file> <output file> <options>"<<endl;
            cout << desc << "\n";
            return 0;
        }

        if(USING_FEATURES) {
            cout << "Using feature detection" << endl;
            USING_FEATURES = true;
        } else {
            cout << "Not using feature detection" << endl;
            USING_FEATURES = false;
        }

        if(USING_CUDA) {
            cout << "Using CUDA" << endl;
            USING_CUDA = true;
        } else {
            cout << "Not using CUDA" << endl;
            USING_CUDA = false;
        }
        if(vm.count("feature-threshold")) {
            FEATURE_THRESHOLD = vm["feature-threshold"].as<float>();
        }
        cout << "Feature threshold is " << FEATURE_THRESHOLD  << endl;

        if(vm.count("corner-threshold")) {
            CORNER_THRESHOLD = vm["corner-threshold"].as<float>();
        }
        cout << "Corner threshold is " << CORNER_THRESHOLD << endl;
    }
    catch(std::exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }

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
    gridPtr F = optimize(volume, featureMap, USING_FEATURES, USING_CUDA);

    /* Extract mesh and write to file */
    mcubes(F, surfaceMap, normals, 0.0, FEATURE_THRESHOLD, CORNER_THRESHOLD, output_path.c_str(), USING_FEATURES);

    return 1;
}
