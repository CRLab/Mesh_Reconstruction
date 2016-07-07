
#include "dfields.h"

using namespace std;



//extract perimeter of binary volume
gridPtr fastPerim(gridPtr volume_grid){
    gridPtr g(new grid(volume_grid->dims, volume_grid->t_));
    for(int i=0; i<g->dims[0]; i++){
        for(int j=0; j<g->dims[1]; j++){
            for(int k=0; k<g->dims[2]; k++){
                (*g)[i][j][k]=1.0;
                int neg_i=1;int neg_j=1;int neg_k=1;
                int pos_i=1; int pos_j=1; int pos_k=1;
                if(i==0)neg_i=0;
                else if(i==g->dims[0]-1)pos_i=0;
                if(j==0)neg_j=0;
                else if(j==g->dims[1]-1)pos_j=0;
                if(k==0)neg_k=0;
                else if(k==g->dims[2]-1)pos_k=0;
                //search neighboring voxels for different value
                for(int i_=i-neg_i; i_<=i+pos_i; i_++){
                    for(int j_=j-neg_j; j_<=j+pos_j; j_++){
                        for(int k_=k-neg_k; k_<=k+pos_k; k_++){
                            if((*volume_grid)[i][j][k]!=(*volume_grid)[i_][j_][k_]){
                                (*g)[i][j][k]=0.0;
                            }
                        }
                    }
                }
                //end of neighbor search
            }
        }
    }

    return g;
}

//implementation of Saito and Toriwaki distance field
//squared euclidean distance to closest 0-value voxel
gridPtr getsqdist(gridPtr volume_grid){
    //make copies
    gridPtr dist_grid = copyGrid(volume_grid);
    gridPtr dist_grid_copy = copyGrid(dist_grid);
    //first transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                //get squared dist to closest 0 in row (i's)
                float min = 100000;
                for(int index=0; index<dist_grid->dims[0]; index++){
                    float dist = 100000;
                    if((*dist_grid)[index][j][k]==0.0){
                        dist = (index-i)*(index-i);
                    }
                    if(dist<min){
                        min=dist;
                    }
                }
                //set value of voxel to min q dist
                (*dist_grid_copy)[i][j][k]=min;
            }
        }
    }//end of first transformation


    dist_grid=copyGrid(dist_grid_copy);
    //second transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                int min = 100000;
                for(int index=0; index<dist_grid->dims[1]; index++){
                    int dist = (*dist_grid)[i][index][k] + (index-j)*(index-j);
                    if(dist<min) min=dist;
                }
                //set value of voxel to min q dist
                (*dist_grid_copy)[i][j][k]=(float)min;
            }
        }
    }//end of second transformation

    dist_grid=copyGrid(dist_grid_copy);
    //third transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                int min = 100000;
                for(int index=0; index<dist_grid->dims[2]; index++){
                    int dist = (*dist_grid)[i][j][index] + (index-k)*(index-k);
                    if(dist<min) min=dist;
                }
                //set value of voxel to min q dist
                (*dist_grid_copy)[i][j][k]=(float)min;
            }
        }
    }//end of third transformation

    return dist_grid_copy;
}

//returns linear indexes of closest 0-value voxel
gridPtr getsqdist_index(gridPtr volume_grid){
    //make copies
    gridPtr index_grid;
    gridPtr index_grid_copy = copyGrid(volume_grid);
    gridPtr dist_grid = copyGrid(volume_grid);
    gridPtr dist_grid_copy = copyGrid(dist_grid);
    //first transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                //get squared dist to closest 0 in row (i's)
                float min = 100000;
                float x_ind = -1.0;
                for(int index=0; index<dist_grid->dims[0]; index++){
                    float dist = 100000;
                    if((*dist_grid)[index][j][k]==0.0){
                        dist = (index-i)*(index-i);
                    }
                    if(dist<min){
                        min=dist;
                        x_ind = index;
                    }
                }
                //set value of voxel to min q dist
                (*dist_grid_copy)[i][j][k]=min;
                //set index value to x index
                (*index_grid_copy)[i][j][k]=x_ind;
            }
        }
    }//end of first transformation

    index_grid=copyGrid(index_grid_copy);
    dist_grid=copyGrid(dist_grid_copy);
    //second transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                int min = 100000;
                float xy_ind = -1.0;
                for(int index=0; index<dist_grid->dims[1]; index++){
                    int dist = (*dist_grid)[i][index][k] + (index-j)*(index-j);
                    if(dist<min){
                        min=dist;
                        xy_ind = (float)(index*index_grid->dims[0])+(*index_grid)[i][index][k];
                    }
                }
                //set value of voxel to min q dist
                (*dist_grid_copy)[i][j][k]=(float)min;
                (*index_grid_copy)[i][j][k] = xy_ind;
            }
        }
    }//end of second transformation

    index_grid=copyGrid(index_grid_copy);
    dist_grid=copyGrid(dist_grid_copy);
    //third transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                int min = 100000;
                float xyz_ind=-1.0;
                for(int index=0; index<dist_grid->dims[2]; index++){
                    int dist = (*dist_grid)[i][j][index] + (index-k)*(index-k);
                    if(dist<min){
                        min=dist;
                        xyz_ind = (float)(index*index_grid->dims[0]*index_grid->dims[1])+(*index_grid)[i][j][index];
                    }
                }
                //set value of voxel to min q dist
                (*dist_grid_copy)[i][j][k]=(float)min;
                (*index_grid_copy)[i][j][k] = xyz_ind;
            }
        }
    }//end of third transformation

    return index_grid_copy;

}

//get square root of grid
gridPtr getsqrt(gridPtr g){
    gridPtr s(new grid(g->dims, g->t_));
    for(int i=0; i<s->dims[0]; i++){
        for(int j=0; j<s->dims[1]; j++){
            for(int k=0; k<s->dims[2]; k++){
                (*s)[i][j][k]= (float) sqrt((double)(*g)[i][j][k]);
            }
        }
    }
    return s;
}

//get linear indices of neighboring voxels
vector<int> getNeighbors(gridPtr g, const Eigen::Vector3i &pnt){
    vector<int> indexes;
    for(int i=pnt[0]-1; i<=pnt[0]+1; i++){
        for(int j=pnt[1]-1; j<=pnt[1]+1; j++){
            for(int k=pnt[2]-1; k<=pnt[2]+1; k++){
                if((*g)[i][j][k]>0){
                    Eigen::Vector3i p;
                    p[0]=i; p[1]=j; p[2]=k;
                    indexes.push_back(g->sub2ind(p));
                }
            }
        }
    }
    return indexes;
}

//check if voxel is on surface
bool isSurface(gridPtr g, const Eigen::Vector3i &pnt){
    //search neighboring voxels for different value
    int i=pnt[0];
    int j=pnt[1];
    int k=pnt[2];
    if((*g)[i][j][k]!=1.0) return false;
    for(int i_=i-1; i_<=i+1; i_++){
        for(int j_=j-1; j_<=j+1; j_++){
            for(int k_=k-1; k_<=k+1; k_++){
                if((*g)[i_][j_][k_]==0.0){
                    return true;
                }
            }
        }
    }
    return false;
}

//get linear indices of neighboring voxels on the surface
vector<int> getSurfaceNeighbors(gridPtr g, const Eigen::Vector3i &pnt){
    vector<int> indexes;
    for(int i=pnt[0]-1; i<=pnt[0]+1; i++){
        for(int j=pnt[1]-1; j<=pnt[1]+1; j++){
            for(int k=pnt[2]-1; k<=pnt[2]+1; k++){
                if((*g)[i][j][k]>0){
                    Eigen::Vector3i p;
                    p[0]=i; p[1]=j; p[2]=k;
                    if(isSurface(g,p)){
                        indexes.push_back(g->sub2ind(p));
                    }
                }
            }
        }
    }
    return indexes;
}

//get centroid of a list of points
Eigen::Vector3f getCentroid(gridPtr g, const vector<int> &pnts){
    Eigen::Vector3f centroid;
    centroid[0]=0.0; centroid[1]=0.0; centroid[2]=0.0;
    for(int i=0; i<pnts.size(); i++){
        Eigen::Vector3i p = g->ind2sub(pnts[i]);
        centroid[0]+=(float)p[0]; centroid[1]+=(float)p[1]; centroid[2]+=(float)p[2];
    }
    centroid[0]=centroid[0]/(float)pnts.size();
    centroid[1]=centroid[1]/(float)pnts.size();
    centroid[2]=centroid[2]/(float)pnts.size();
}

//get covariance matrix of surface point
Eigen::Matrix3f getCovariance(gridPtr g, const Eigen::Vector3i &pnt){
    vector<int> indexes = getSurfaceNeighbors(g, pnt);
    Eigen::Vector3f centroid = getCentroid(g, indexes);
    Eigen::Matrix3f covariance;
    covariance<<0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    for(int i=0; i<indexes.size(); i++){
        Eigen::Vector3f vec;
        vec[0]=(float)(g->ind2sub(indexes[i]))[0] - centroid[0];
        vec[1]=(float)(g->ind2sub(indexes[i]))[1] - centroid[1];
        vec[2]=(float)(g->ind2sub(indexes[i]))[2] - centroid[2];
        for(int row=0; row<3; row++){
            for(int col=0; col<3; col++){
                covariance(row,col) = covariance(row,col)+vec[row]*vec[col];
            }
        }
    }
    for(int row=0; row<3; row++){
        for(int col=0; col<3; col++){
            covariance(row,col) = covariance(row,col)/(float)indexes.size();
        }
    }
}

//get normal vector from covariance matrix
Eigen::Vector3f getNormalVector(const Eigen::Matrix3f &covariance){
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> es;
    es.compute(covariance);
    float min=10000;
    int index=0;
    for(int i=0; i<3; i++){
        int eigenvalue=es.eigenvalues()[i];
        if(eigenvalue<min){
            min=eigenvalue;
            index=i;
        }
    }
    Eigen::Vector3f vec = es.eigenvectors().col(index);
    float len = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
    vec[0]=vec[0]/len; vec[1]=vec[1]/len; vec[2]=vec[2]/len;
    return vec;
}

//orient normal vector with respect to local centroid
void orientNormal(gridPtr g, Eigen::Vector3f &normal, const Eigen::Vector3i &pnt){
    Eigen::Vector3f centroid = getCentroid(g, getNeighbors(g, pnt));
    centroid[0]-=(float)pnt[0]; centroid[1]-=(float)pnt[1]; centroid[2]-=(float)pnt[2];
    float prod = normal[0]*centroid[0]+normal[1]*centroid[1]+normal[2]*centroid[2];
    if(prod<0.0){
        normal[0]=-normal[0]; normal[1]=-normal[1]; normal[2]=-normal[2];
    }
}

//get indexes of surface points of volume grid
vector<int> getSurface(gridPtr volume){
    vector<int> indexes;
    gridPtr surface(new grid(volume->dims, volume->t_));
    for(int i=0; i<surface->dims[0]; i++){
        for(int j=0; j<surface->dims[1]; j++){
            for(int k=0; k<surface->dims[2]; k++){
                Eigen::Vector3i pnt;
                pnt[0]=i; pnt[1]=j; pnt[2]=k;
                if(isSurface(volume, pnt)){
                    indexes.push_back(volume->sub2ind(pnt));
                }
            }
        }
    }
    return indexes;
}

//get surface normals
vector<Eigen::Vector3f> getSurfaceNormals(gridPtr volume, const vector<int> &surface){
    vector<Eigen::Vector3f> normals;
    for(int i=0; i<surface.size(); i++){
        normals.push_back(getNormalVector(getCovariance(volume, volume->ind2sub(surface[i]))));
        orientNormal(volume, normals[i], volume->ind2sub(surface[i]));
    }
}

