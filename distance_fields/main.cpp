#include "dfields.h"

int main(int argc, char **argv){
    //testing normal finding
    Eigen::Vector3i dims;
    Eigen::Vector3i t_;
    int n=5;
    dims[0]=n;dims[1]=n;dims[2]=n;
    t_[0]=0;t_[1]=0;t_[2]=0;
    gridPtr volume (new grid(dims, t_));
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            for(int k=0; k<n; k++){
                if(i>0 && i<n-1 && j>0 && j<n-1 && k>0 && k<n-1){
                    (*volume)[i][j][k]=1.0;
                }
                else (*volume)[i][j][k]=0.0;
            }
        }
    }
    //volume->visualize();
    visualizeNormals(volume);
    vector<int> surface = getSurface(volume);
    gridPtr surfaceMap = getIndexMap(volume, surface);
    vector<Eigen::Vector3f> normals = getSurfaceNormals(volume, surface);

    gridPtr featureMap = getFeatureMap(volume, surfaceMap, normals, 0.85);

    featureMap->visualize();

    return 1;
}
