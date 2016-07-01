
#include "getsqdist.h"

using namespace std;

gridPtr createGrid(pcl::PointCloud<pcl::InterestPoint>::Ptr grid_cloud, VoxelGridPtr vox, int res_factor){
    Eigen::Vector3i t_ = vox->getMinBoxCoordinates();
    Eigen::Vector3i dims = vox->getNrDivisions();
    //pad the grid with empty voxels 5 on each side
    int pad= 6;
    dims[0]+=pad*2; dims[1]+=pad*2; dims[2]+=pad*2;
    t_[0]-=pad; t_[1]-=pad; t_[2]-=pad;
    gridPtr gp(new grid(dims, t_));
    for(int i=0; i<gp->dims[0]; i++){
        for(int j=0; j<gp->dims[1]; j++){
            for(int k=0; k<gp->dims[2]; k++){
                Eigen::Vector3i pnt;
                pnt[0]=i+gp->t_[0]; pnt[1]=j+gp->t_[1]; pnt[2]=k+gp->t_[2];
                int index = vox->getCentroidIndexAt(pnt);
                if(index!=-1) (*gp)[i][j][k] = grid_cloud->points[index].strength;
                else (*gp)[i][j][k]=-1.0;

                if(i<pad||i>=gp->dims[0]-pad||j<pad||j>=gp->dims[1]-pad||k<pad||k>=gp->dims[2]-pad){
                    (*gp)[i][j][k]=-1.0;
                }
            }
        }
    }

    //fill in any holes between observed and completed clouds
    gp->fillGrid(res_factor);
    return gp;
}

//create binary volume grid from confidence grid
gridPtr getBinaryVolume(gridPtr grid_cloud){
    gridPtr g(new grid(grid_cloud->dims, grid_cloud->t_));
    for(int i=0; i<g->dims[0]; i++){
        for(int j=0; j<g->dims[1]; j++){
            for(int k=0; k<g->dims[2]; k++){
                if((*grid_cloud)[i][j][k]>=0) (*g)[i][j][k]=1.0;
                else (*g)[i][j][k]=0.0;
            }
        }
    }
    return g;
}

//copy grid
gridPtr copyGrid(gridPtr in){
    gridPtr g(new grid(*in));
    return g;
}

//add 2 grids voxel by voxel
gridPtr addGrids(gridPtr in1, gridPtr in2){
    gridPtr g(new grid(in1->dims, in1->t_));
    *g = ((*in1)+(*in2));
    return g;
}

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

