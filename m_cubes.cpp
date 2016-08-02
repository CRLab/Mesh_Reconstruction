
#include "assign_confidence.h"

#include "ConstConf.h"

#include "voxelize.h"

#include "dfields.h"

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
    vector<string> list;
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error opening directory" << dir << endl;
        exit(1);
    }

    while ((dirp = readdir(dp)) != NULL) {
        list.push_back(string(dirp->d_name));
    }
    closedir(dp);
    vector<string> dirs;
    for(int i=0; i<list.size(); i++){
        string d (dir+list[i]+"/");
        if(hasDirectory(d)){
            if(list[i].compare(".")!=0 && list[i].compare("..")!=0){
                dirs.push_back(list[i]+"/");
            }
        }
    }
    return dirs;
}

vector<string> getfiles (string dir){
    vector<string> list;
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error opening directory" << dir << endl;
        exit(1);
    }

    while ((dirp = readdir(dp)) != NULL) {
        list.push_back(string(dirp->d_name));
    }
    closedir(dp);
    vector<string> files;
    for(int i=0; i<list.size(); i++){
        if(list[i].compare(".")!=0 && list[i].compare("..")!=0){
            files.push_back(list[i]);
        }
    }
    return files;
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
    //feature values not used
    const float FEATURE_THRESHOLD = 0.75;
    const float CORNER_THRESHOLD = 0.8;

    //Highest level directories for data input and output
    string DATA_DIR ("/home/adamjri/aligned_completion_results/");
    string MESH_DIR ("/home/adamjri/completed_meshes/");
    //************************************************************

    /* Prompt user for shape directory */
    string SHAPE_NAME="";
    cout<<"Here are the shape directories:"<<endl;
    vector<string> shape_dirs = getdirs(DATA_DIR);
    for(int i=0; i<shape_dirs.size(); i++){
        cout<<shape_dirs[i]<<endl;
    }
    while(SHAPE_NAME.length()==0){
        cout<<"Please enter the name of the shape (with '/' at end): "<<endl;
        getline(cin, SHAPE_NAME);
        if(SHAPE_NAME.length()==0) cout<<"Shape not entered"<<endl;
    }
    if(!hasDirectory(DATA_DIR+SHAPE_NAME)){
        cout<<"Could not find directory: "<<SHAPE_NAME<<endl;
        exit(1);
    }

    /* Prompt user for capture angle */
    string SHAPE_CAPTURE_ANGLE="";
    cout<<"Here are the angle directories:"<<endl;
    vector<string> angle_dirs = getdirs(DATA_DIR+SHAPE_NAME);
    for(int j=0; j<angle_dirs.size(); j++){
        cout<<angle_dirs[j]<<endl;
    }
    while(SHAPE_CAPTURE_ANGLE.length()==0){
        cout<<"Please enter the capture angle directory (with '/' at end)."<<endl;
        getline(cin, SHAPE_CAPTURE_ANGLE);
        if(SHAPE_CAPTURE_ANGLE.length()==0) cout<<"Capture angle not entered"<<endl;
    }
    if(!hasDirectory(DATA_DIR+SHAPE_NAME+SHAPE_CAPTURE_ANGLE)){
        cout<<"Could not find directory: "<<SHAPE_CAPTURE_ANGLE<<endl;
        exit(1);
    }

    /* Prompt user for pcd file */
    string PCD_FILE="";
    cout<<"Here are the files:"<<endl;
    vector<string> files = getfiles(DATA_DIR+SHAPE_NAME+SHAPE_CAPTURE_ANGLE);
    for(int j=0; j<files.size(); j++){
        cout<<files[j]<<endl;
    }
    while(PCD_FILE.length()==0){
        cout<<"Please enter the filename (must be .pcd)."<<endl;
        getline(cin, PCD_FILE);
        if(PCD_FILE.length()==0) cout<<"Filename not entered"<<endl;
    }

    string PCD_PATH(DATA_DIR+SHAPE_NAME+SHAPE_CAPTURE_ANGLE+PCD_FILE);

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
    float resolution = getResolution(confPCL);
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

    /* Create directory for output file */
    createDirectory(MESH_DIR+SHAPE_NAME);
    string out_file (SHAPE_CAPTURE_ANGLE);
    out_file = out_file.substr(0, out_file.length()-1);
    string PCL_FILENAME = PCD_FILE.substr(2,PCD_FILE.length()-6);
    out_file+="_"+PCL_FILENAME+"_mesh.ply";
    string out_path (MESH_DIR+SHAPE_NAME+out_file);

    /* Extract mesh and write to file */
    mcubes(volume, surfaceMap, normals, 0.5, FEATURE_THRESHOLD, CORNER_THRESHOLD, out_path.c_str(), false);
}
