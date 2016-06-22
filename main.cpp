#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/surface/mls.h>
#include <boost/thread/thread.hpp>
#include <pcl/common/common_headers.h>
#include <pcl/features/normal_3d.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include "voxelize.h"

#include "binvoxToPcl.h"
#include "assign_confidence.h"

#include "ConstConf.h"

#include "getsqdist.h"
#include "narrowBand.h"

#include "optimize/quadprog.h"

#include <iostream>

using namespace std;
using namespace lemp;
//using namespace conf;

typedef unsigned char byte;

int main(int argc, char **argv){
    //create sphere volume
    gridPtr volume (new grid);
    volume->dims[0]=64;volume->dims[1]=64;volume->dims[2]=64;
    volume->voxels = allocGrid(volume->dims);
    float p = (float)floor(64.0*0.33);
    p=p*p;
    for(int i=0; i<64; i++){
        for(int j=0; j<64; j++){
            for(int k=0; k<64; k++){
                float dist = ((float)i-31.5)*((float)i-31.5);
                dist+= ((float)j-31.5)*((float)j-31.5);
                dist+= ((float)k-31.5)*((float)k-31.5);
                if(dist<p){
                    volume->voxels[i][j][k]=1.0;
                }
                else{
                    volume->voxels[i][j][k]=0.0;
                }
            }
        }
    }


    //get imbedding function
    gridPtr F = optimize(volume);
    //gridPtr F = volume;

    //write to file
    ofstream myfile;
    myfile.open("embed_func_sphere.txt");
    for(int i=0; i<F->dims[0]; i++){
        for(int j=0; j<F->dims[1]; j++){
            for(int k=0; k<F->dims[2]; k++){
                myfile<<i<<","<<j<<","<<k<<","<<F->voxels[i][j][k]<<endl;
            }
        }
    }
    myfile.close();


    //visualize
    visualizeGrid(F);

    return 1;
}
