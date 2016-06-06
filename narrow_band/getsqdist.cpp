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
    g->dims = vox->getNrDivisions();
    g->voxels = allocGrid(g->dims);
    for(int i=0; i<g->dims[0]; i++){
        for(int j=0; j<g->dims[1]; j++){
            for(int k=0; k<g->dims[2]; k++){
                Eigen::Vector3i pnt;
                pnt[0]=i; pnt[1]=j; pnt[2]=k;
                int index = vox->getCentroidIndexAt(pnt);
                cout<<index<<endl<<endl;
                if(index!=-1) g->voxels[i][j][k] = grid_cloud->points[index].strength;
                else g->voxels[i][j][k]=0.0;
            }
        }
    }
    return g;
}

//create binary volume grid from confidence grid
gridPtr getBinaryVolume(gridPtr grid_cloud){
    gridPtr g = gridPtr(new grid());
    g->dims = grid_cloud->dims;
    g->voxels = allocGrid(g->dims);
    for(int i=0; i<g->dims[0]; i++){
        for(int j=0; j<g->dims[1]; j++){
            for(int k=0; k<g->dims[2]; k++){
                if(grid_cloud->voxels[i][j][k]>0) g->voxels[i][j][k]=1.0;
                else g->voxels[i][j][k]=0.0;
            }
        }
    }
    return g;
}

//get negative of grid
gridPtr getNegGrid(gridPtr in){
    gridPtr g = gridPtr(new grid());
    g->dims = in->dims;
    g->voxels = allocGrid(g->dims);
    for(int i=0; i<g->dims[0]; i++){
        for(int j=0; j<g->dims[1]; j++){
            for(int k=0; k<g->dims[2]; k++){
                g->voxels[i][j][k]=0.0-in->voxels[i][j][k];
            }
        }
    }
    return g;
}

//add 2 grids voxel by voxel
gridPtr addGrids(gridPtr in1, gridPtr in2){
    gridPtr g = gridPtr(new grid());
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
    g->dims = volume_grid->dims;
    g->voxels = allocGrid(g->dims);
    for(int i=1; i<g->dims[0]-1; i++){
        for(int j=1; j<g->dims[1]-1; j++){
            for(int k=1; k<g->dims[2]-1; k++){
                //search neighboring voxels for different value
                for(int i_=i-1; i_<=i+1; i_++){
                    for(int j_=j-1; j_<=j+1; j_++){
                        for(int k_=k-1; k_<=k+1; k_++){
                            g->voxels[i][j][k]=1.0;
                            if(volume_grid->voxels[i][j][k]!=volume_grid->voxels[i_][j_][k_]){
                                g->voxels[i][j][k]=0.0;
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
gridPtr getsqdist(gridPtr volume_grid){
    //extract perimeter
    gridPtr dist_grid = fastPerim(volume_grid);
    //first transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                //get squared dist to closest 0 in row (i's)
                int min = 10000000;
                for(int index=0; index<dist_grid->dims[0]; index++){
                    int dist = 10000000;
                    if(volume_grid->voxels[index][j][k]==0.0){
                        dist = (index-i)*(index-i);
                    }
                    if(dist<min) min=dist;
                }
                //set value of voxel to min q dist
                dist_grid->voxels[i][j][k]=(float)min;
            }
        }
    }//end of first transformation
    //second transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                int min = 10000000;
                for(int index=0; index<dist_grid->dims[1]; index++){
                    int dist = dist_grid->voxels[i][index][k] + (index-j)*(index-j);
                    if(dist<min) min=dist;
                }
                //set value of voxel to min q dist
                dist_grid->voxels[i][j][k]=(float)min;
            }
        }
    }//end of second transformation
    //third transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                int min = 10000000;
                for(int index=0; index<dist_grid->dims[2]; index++){
                    int dist = dist_grid->voxels[i][j][index] + (index-k)*(index-k);
                    if(dist<min) min=dist;
                }
                //set value of voxel to min q dist
                dist_grid->voxels[i][j][k]=(float)min;
            }
        }
    }//end of third transformation

    return dist_grid;

}

//get square root of grid
gridPtr getsqrt(gridPtr g){
    gridPtr s = gridPtr(new grid());
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

