#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/filters/voxel_grid.h>

#include <iostream>

#include "getsqdist.h"

using namespace std;

float*** allocGrid(Eigen::Vector3i dims){
    float*** arr = new float**[dims[0]];
    for(int i=0; i<dims[0]; i++){
        arr[i] = new float*[dims[1]];
        for(int j=0; j<dims[1]; j++){
            arr[i][j] = new float[dims[2]];
        }
    }
    return arr;
}

//create grid with voxel values=confidence values of point cloud
gridPtr createGrid(pcl::PointCloud<pcl::InterestPoint>::Ptr grid_cloud, VoxelGridPtr vox){
    gridPtr g = gridPtr(new grid());
    g->t_ = vox->getMinBoxCoordinates();
    g->dims = vox->getNrDivisions();
    //pad the grid with empty voxels 5 on each side
    int pad= 5;
    g->dims[0]+=pad*2; g->dims[1]+=pad*2; g->dims[2]+=pad*2;
    g->t_[0]-=pad; g->t_[1]-=pad; g->t_[2]-=pad;
    g->voxels = allocGrid(g->dims);
    for(int i=0; i<g->dims[0]; i++){
        for(int j=0; j<g->dims[1]; j++){
            for(int k=0; k<g->dims[2]; k++){
                Eigen::Vector3i pnt;
                pnt[0]=i+g->t_[0]; pnt[1]=j+g->t_[1]; pnt[2]=k+g->t_[2];
                int index = vox->getCentroidIndexAt(pnt);
                if(index!=-1) g->voxels[i][j][k] = grid_cloud->points[index].strength;
                else g->voxels[i][j][k]=-1.0;

                if(i<pad||i>=g->dims[0]-pad||j<pad||j>=g->dims[1]-pad||k<pad||k>=g->dims[2]-pad){
                    g->voxels[i][j][k]=-1.0;
                }
            }
        }
    }
    return g;
}

//create binary volume grid from confidence grid
gridPtr getBinaryVolume(gridPtr grid_cloud){
    gridPtr g = gridPtr(new grid());
    g->t_ = grid_cloud->t_;
    g->dims = grid_cloud->dims;
    g->voxels = allocGrid(g->dims);
    for(int i=0; i<g->dims[0]; i++){
        for(int j=0; j<g->dims[1]; j++){
            for(int k=0; k<g->dims[2]; k++){
                if(grid_cloud->voxels[i][j][k]>=0) g->voxels[i][j][k]=1.0;
                else g->voxels[i][j][k]=0.0;
            }
        }
    }
    return g;
}

//copy grid
gridPtr copyGrid(gridPtr in){
    gridPtr g = gridPtr(new grid());
    g->t_=in->t_;
    g->dims = in->dims;
    g->voxels = allocGrid(g->dims);
    for(int i=0; i<g->dims[0]; i++){
        for(int j=0; j<g->dims[1]; j++){
            for(int k=0; k<g->dims[2]; k++){
                g->voxels[i][j][k]=in->voxels[i][j][k];
            }
        }
    }
    return g;
}
//get negative of grid
gridPtr getNegGrid(gridPtr in){
    gridPtr g = gridPtr(new grid());
    g->t_ = in->t_;
    g->dims = in->dims;
    g->voxels = allocGrid(g->dims);
    for(int i=0; i<g->dims[0]; i++){
        for(int j=0; j<g->dims[1]; j++){
            for(int k=0; k<g->dims[2]; k++){
                g->voxels[i][j][k]=-in->voxels[i][j][k];
            }
        }
    }
    return g;
}

//add 2 grids voxel by voxel
gridPtr addGrids(gridPtr in1, gridPtr in2){
    gridPtr g = gridPtr(new grid());
    g->t_ = in1->t_;
    g->dims = in1->dims;
    g->voxels = allocGrid(g->dims);
    for(int i=0; i<g->dims[0]; i++){
        for(int j=0; j<g->dims[1]; j++){
            for(int k=0; k<g->dims[2]; k++){
                g->voxels[i][j][k]=in1->voxels[i][j][k]+in2->voxels[i][j][k];
            }
        }
    }
    return g;
}

//extract perimeter of binary volume
gridPtr fastPerim(gridPtr volume_grid){
    gridPtr g = gridPtr(new grid());
    g->t_ = volume_grid->t_;
    g->dims = volume_grid->dims;
    g->voxels = allocGrid(g->dims);
    for(int i=0; i<g->dims[0]; i++){
        for(int j=0; j<g->dims[1]; j++){
            for(int k=0; k<g->dims[2]; k++){
                g->voxels[i][j][k]=1.0;
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
                                if(volume_grid->voxels[i][j][k]!=volume_grid->voxels[i_][j_][k_]){
                                    g->voxels[i][j][k]=0.0;
                                }
                            }
                        }
                    }
                    //end of neighbor search
                //}
            }
        }
    }

    return g;
}

//implementation of Saito and Toriwaki distance field
gridPtr getsqdist(gridPtr volume_grid){
    //extract perimeter
    gridPtr dist_grid = fastPerim(volume_grid);
    //make a copy
    gridPtr dist_grid_copy = copyGrid(dist_grid);
    //first transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                //get squared dist to closest 0 in row (i's)
                float min = 100000;
                for(int index=0; index<dist_grid->dims[0]; index++){
                    float dist = 100000;
                    if(dist_grid->voxels[index][j][k]==0.0){
                        dist = (index-i)*(index-i);
                    }
                    if(dist<min){
                        min=dist;
                    }
                }
                //set value of voxel to min q dist
                dist_grid_copy->voxels[i][j][k]=min;
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
                    int dist = dist_grid->voxels[i][index][k] + (index-j)*(index-j);
                    if(dist<min) min=dist;
                }
                //set value of voxel to min q dist
                dist_grid_copy->voxels[i][j][k]=(float)min;
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
                    int dist = dist_grid->voxels[i][j][index] + (index-k)*(index-k);
                    if(dist<min) min=dist;
                }
                //set value of voxel to min q dist
                dist_grid_copy->voxels[i][j][k]=(float)min;
            }
        }
    }//end of third transformation

    return dist_grid_copy;

}

//get square root of grid
gridPtr getsqrt(gridPtr g){
    gridPtr s = gridPtr(new grid());
    s->t_ = g->t_;
    s->dims = g->dims;
    s->voxels = allocGrid(s->dims);
    for(int i=0; i<s->dims[0]; i++){
        for(int j=0; j<s->dims[1]; j++){
            for(int k=0; k<s->dims[2]; k++){
                s->voxels[i][j][k]= (float) sqrt((double)g->voxels[i][j][k]);
            }
        }
    }
    return s;
}

