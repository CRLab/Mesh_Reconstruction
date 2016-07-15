
#include "voxelize.h"

#include "binvoxToPcl.h"
#include "assign_confidence.h"

#include "ConstConf.h"

#include "narrowBand.h"

#include "quadprog.h"

#include "mcubes.h"

#include <iostream>
#include <string>
#include <errno.h>
#include <dirent.h>
#include <sys/stat.h>

using namespace std;

//these directory handling functions only work on Linux
//*******************************************************************
bool hasDirectory(string dir){
    struct stat sb;
    return (stat(dir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode));
}

vector<string> getdirs (string dir){
    vector<string> files;
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error opening directory" << dir << endl;
        exit(1);
    }

    while ((dirp = readdir(dp)) != NULL) {
        files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    vector<string> dirs;
    for(int i=0; i<files.size(); i++){
        string d (dir+files[i]+"/");
        if(hasDirectory(d)){
            if(files[i].compare(".")!=0 && files[i].compare("..")!=0){
                dirs.push_back(files[i]+"/");
            }
        }
    }
    return dirs;
}

void createDirectory(string dir){
    if(!hasDirectory(dir)){
        const int dir_err = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err)
        {
            cerr<<"Error creating directory!"<<endl;
            exit(1);
        }
    }
}
//*****************************************************************



int main(int argc, char **argv){

    //************************************************************
    //variables for allowing processing of all data at once
    const bool CAPTURE_ALL_SHAPES = false;
    bool CAPTURE_ALL_SHAPE_DATA = false;

    if(CAPTURE_ALL_SHAPES) CAPTURE_ALL_SHAPE_DATA=true;

    //variable for toggling feature handling
    const bool USING_FEATURES = false;

    //thresholds for feature detection
    const float FEATURE_THRESHOLD = 0.85;
    const float CORNER_THRESHOLD = 0.7;

    //Highest level directories for data input and output
    string DATA_DIR ("/home/adamjri/aligned_completion_results/");
    string MESH_DIR ("/home/adamjri/completed_meshes/");

    //Lowest level filenames for pcd and binvox files
    string PCD_FILE ("c_partial.pcd");
    string BINVOX_FILE ("c_completed.binvox");
    //************************************************************

    //names of directories for shapes and shape capture angles
    vector<string> SHAPE_NAMES;
    vector<vector<string> > SHAPE_CAPTURE_ANGLES;

    /* Get shape directories */
    if(CAPTURE_ALL_SHAPES){
        /* Automatically get all shape directories */
        //get list of directories of DATA_DIR
        SHAPE_NAMES = getdirs(DATA_DIR);
    }
    else{
        /* Prompt user for shape directory */
        string shape_name="";
        cout<<"Here are the shape directories:"<<endl;
        vector<string> sub_dirs = getdirs(DATA_DIR);
        for(int i=0; i<sub_dirs.size(); i++){
            cout<<sub_dirs[i]<<endl;
        }
        while(shape_name.length()==0){
            cout<<"Please enter the name of the shape (with '/' at end): "<<endl;
            getline(cin, shape_name);
            if(shape_name.length()==0)cout<<"Shape not entered"<<endl;
        }
        if(!hasDirectory(DATA_DIR+shape_name)){
            cout<<"Could not find directory: "<<shape_name<<endl;
            exit(1);
        }
        SHAPE_NAMES.push_back(shape_name);
    }

    /* Get capture angle directories */
    for(int i=0; i<SHAPE_NAMES.size(); i++){
        if(CAPTURE_ALL_SHAPE_DATA){
            /* Automatically get all shape capture angle directories */
            //get list of directories for capture angles of shape
            SHAPE_CAPTURE_ANGLES.push_back(getdirs(DATA_DIR+SHAPE_NAMES[i]));
        }
        else{
            /* Prompt user for capture angle */
            //no capture angle may be selected
            string angle_name="";
            cout<<"Here are the angle directories:"<<endl;
            vector<string> sub_dirs = getdirs(DATA_DIR+SHAPE_NAMES[i]);
            for(int j=0; j<sub_dirs.size(); j++){
                cout<<sub_dirs[j]<<endl;
            }
            cout<<"Please enter the capture angle directory (with '/' at end).You may enter nothing: "<<endl;
            getline(cin, angle_name);
            if(angle_name.length()==0){
                //remove shape from shape list
                SHAPE_NAMES.erase(SHAPE_NAMES.begin()+i);
                i--;
                continue;
            }
            if(!hasDirectory(DATA_DIR+SHAPE_NAMES[i]+angle_name)){
                cout<<"Could not find directory: "<<angle_name<<endl;
                exit(1);
            }
            vector<string> names;
            names.push_back(angle_name);
            SHAPE_CAPTURE_ANGLES.push_back(names);
        }
    }

    //**************************************************************************************
    //perform algorithm
    for(int shape=0; shape<SHAPE_NAMES.size(); shape++){
        for(int angle=0; angle<SHAPE_CAPTURE_ANGLES[shape].size(); angle++){

            /* Get full path to pcd and binvox files */
            string pcd_path (DATA_DIR+SHAPE_NAMES[shape]+SHAPE_CAPTURE_ANGLES[shape][angle]+PCD_FILE);
            string binvox_path (DATA_DIR+SHAPE_NAMES[shape]+SHAPE_CAPTURE_ANGLES[shape][angle]+BINVOX_FILE);


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
            voxelized_dataPtr data = voxelizeData(confPCL); //<--test different resolutions
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

            /* Create directory for output file */
            createDirectory(MESH_DIR+SHAPE_NAMES[shape]);
            string out_file (SHAPE_CAPTURE_ANGLES[shape][angle]);
            out_file = out_file.substr(0, out_file.length()-1);
            out_file+= "_mesh.ply";
            string out_path (MESH_DIR+SHAPE_NAMES[shape]+out_file);

            /* Extract mesh and write to file */
            mcubes(F, surfaceMap, normals, 0.0, FEATURE_THRESHOLD, CORNER_THRESHOLD, out_path.c_str(), USING_FEATURES);
        }
    }

    return 1;
}
