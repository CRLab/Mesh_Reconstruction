
#include "narrowBand.h"

using namespace std;

typedef boost::shared_ptr< pcl::VoxelGrid<pcl::InterestPoint> > VoxelGridPtr;
typedef boost::shared_ptr<bands> bandsPtr;

//morphological erosion with mask: [[000;010;000],[010,111,010],[000;010;000]]
gridPtr erode_grid(gridPtr gr){
    gridPtr eroded = copyGrid(gr);
    for(int i=1; i<eroded->dims[0]-1; i++){
        for(int j=1; j<eroded->dims[1]-1; j++){
            for(int k=1; k<eroded->dims[2]-1; k++){
                if((*gr)[i][j][k]==0.0){
                    (*eroded)[i][j][k]=0.0;
                }
                else if((*gr)[i-1][j][k]==0.0){
                    (*eroded)[i][j][k]=0.0;
                }
                else if((*gr)[i][j-1][k]==0.0){
                    (*eroded)[i][j][k]=0.0;
                }
                else if((*gr)[i][j][k-1]==0.0){
                    (*eroded)[i][j][k]=0.0;
                }
                else if((*gr)[i+1][j][k]==0.0){
                    (*eroded)[i][j][k]=0.0;
                }
                else if((*gr)[i][j+1][k]==0.0){
                    (*eroded)[i][j][k]=0.0;
                }
                else if((*gr)[i][j][k+1]==0.0){
                    (*eroded)[i][j][k]=0.0;
                }
                else{
                    (*eroded)[i][j][k]=1.0;
                }
            }
        }
    }
    return eroded;
}

//morphological dilation with mask: [[000;010;000],[010,111,010],[000;010;000]]
gridPtr dilate_grid(gridPtr gr){
    gridPtr dilated = copyGrid(gr);
    for(int i=1; i<dilated->dims[0]-1; i++){
        for(int j=1; j<dilated->dims[1]-1; j++){
            for(int k=1; k<dilated->dims[2]-1; k++){
                if((*gr)[i][j][k]==1.0){
                    (*dilated)[i][j][k]=1.0;
                }
                else if((*gr)[i-1][j][k]==1.0){
                    (*dilated)[i][j][k]=1.0;
                }
                else if((*gr)[i][j-1][k]==1.0){
                    (*dilated)[i][j][k]=1.0;
                }
                else if((*gr)[i][j][k-1]==1.0){
                    (*dilated)[i][j][k]=1.0;
                }
                else if((*gr)[i+1][j][k]==1.0){
                    (*dilated)[i][j][k]=1.0;
                }
                else if((*gr)[i][j+1][k]==1.0){
                    (*dilated)[i][j][k]=1.0;
                }
                else if((*gr)[i][j][k+1]==1.0){
                    (*dilated)[i][j][k]=1.0;
                }
                else{
                    (*dilated)[i][j][k]=0.0;
                }
            }
        }
    }
    return dilated;
}

//generate band and tight band (eroded band) using dist field "margin" and band_size
bandsPtr createBands(gridPtr margin, float band_size){
    bandsPtr bnds = bandsPtr(new bands());
    bnds->band = gridPtr(new grid(margin->dims, margin->t_));
    //create band
    for(int i=0; i<bnds->band->dims[0]; i++){
        for(int j=0; j<bnds->band->dims[1]; j++){
            for(int k=0; k<bnds->band->dims[2]; k++){
                if((*margin)[i][j][k]<=band_size){
                    (*(bnds->band))[i][j][k]=1.0;
                }
                else{
                    (*(bnds->band))[i][j][k]=0.0;
                }
            }
        }
    }
    bnds->tight_band = erode_grid(bnds->band);
    bnds->band=dilate_grid(bnds->tight_band);

    return bnds;
}


