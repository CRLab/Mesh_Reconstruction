#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/filters/voxel_grid.h>

#include <iostream>

#include "getsqdist.h"
#include "narrowBand.h"

using namespace std;

typedef boost::shared_ptr< pcl::VoxelGrid<pcl::InterestPoint> > VoxelGridPtr;
typedef boost::shared_ptr<bands> bandsPtr;

//morphological erosion with mask: [[000;010;000],[010,111,010],[000;010;000]]
gridPtr erode_grid(gridPtr gr){
    gridPtr eroded = gridPtr(new grid());
    eroded->dims = gr->dims;
    eroded->voxels = allocGrid(eroded->dims);
    for(int i=1; i<eroded->dims[0]-1; i++){
        for(int j=1; j<eroded->dims[1]-1; j++){
            for(int k=1; k<eroded->dims[2]-1; k++){
                if(gr->voxels[i-1][j][k]==0.0){
                    eroded->voxels[i][j][k]=0.0;
                }
                else if(gr->voxels[i][j-1][k]==0.0){
                    eroded->voxels[i][j][k]=0.0;
                }
                else if(gr->voxels[i][j][k-1]==0.0){
                    eroded->voxels[i][j][k]=0.0;
                }
                else if(gr->voxels[i+1][j][k]==0.0){
                    eroded->voxels[i][j][k]=0.0;
                }
                else if(gr->voxels[i][j+1][k]==0.0){
                    eroded->voxels[i][j][k]=0.0;
                }
                else if(gr->voxels[i][j][k+1]==0.0){
                    eroded->voxels[i][j][k]=0.0;
                }
                else{
                    eroded->voxels[i][j][k]=gr->voxels[i][j][k];
                }
            }
        }
    }
    return eroded;
}

//morphological dilation with mask: [[000;010;000],[010,111,010],[000;010;000]]
gridPtr dilate_grid(gridPtr gr){
    gridPtr dilated = gridPtr(new grid());
    dilated->dims = gr->dims;
    dilated->voxels = allocGrid(dilated->dims);
    for(int i=1; i<dilated->dims[0]-1; i++){
        for(int j=1; j<dilated->dims[1]-1; j++){
            for(int k=1; k<dilated->dims[2]-1; k++){
                if(gr->voxels[i-1][j][k]==1.0){
                    dilated->voxels[i][j][k]=1.0;
                }
                else if(gr->voxels[i][j-1][k]==1.0){
                    dilated->voxels[i][j][k]=1.0;
                }
                else if(gr->voxels[i][j][k-1]==1.0){
                    dilated->voxels[i][j][k]=1.0;
                }
                else if(gr->voxels[i+1][j][k]==1.0){
                    dilated->voxels[i][j][k]=1.0;
                }
                else if(gr->voxels[i][j+1][k]==1.0){
                    dilated->voxels[i][j][k]=1.0;
                }
                else if(gr->voxels[i][j][k+1]==1.0){
                    dilated->voxels[i][j][k]=1.0;
                }
                else{
                    dilated->voxels[i][j][k]=gr->voxels[i][j][k];
                }
            }
        }
    }
    return dilated;
}

//generate band and tight band (eroded band) using dist field "margin" and band_size
bandsPtr createBands(gridPtr margin, float band_size){
    bandsPtr bnds = bandsPtr(new bands());
    bnds->band = gridPtr(new grid());
    bnds->band->dims=margin->dims;
    bnds->band->voxels=allocGrid(bnds->band->dims);
    //create band
    for(int i=0; i<bnds->band->dims[0]; i++){
        for(int j=0; j<bnds->band->dims[1]; j++){
            for(int k=0; k<bnds->band->dims[2]; k++){
                if(margin->voxels[i][j][k]<=band_size){
                    bnds->band->voxels[i][j][k]=1.0;
                }
                else{
                    bnds->band->voxels[i][j][k]=0.0;
                }
            }
        }
    }

    bnds->tight_band = erode_grid(bnds->band);
    //clear edges
    for(int i=0; i<bnds->tight_band->dims[0]; i++){
        for(int j=0; j<bnds->tight_band->dims[1]; j++){
            bnds->band->voxels[i][j][0]=0.0;
            bnds->band->voxels[i][j][bnds->tight_band->dims[2]-1]=0.0;
        }
    }
    for(int i=0; i<bnds->tight_band->dims[0]; i++){
        for(int k=0; k<bnds->tight_band->dims[2]; k++){
            bnds->band->voxels[i][0][k]=0.0;
            bnds->band->voxels[i][bnds->tight_band->dims[1]-1][k]=0.0;
        }
    }
    for(int j=0; j<bnds->tight_band->dims[1]; j++){
        for(int k=0; k<bnds->tight_band->dims[2]; k++){
            bnds->band->voxels[0][j][k]=0.0;
            bnds->band->voxels[bnds->tight_band->dims[0]-1][j][k]=0.0;
        }
    }

    bnds->band=dilate_grid(bnds->tight_band);

    return bnds;
}


