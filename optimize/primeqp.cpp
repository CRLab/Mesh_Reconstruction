#include <Eigen/Sparse>
#include <iostream>

#include "narrowBand.h"
#include "primeqp.h"

namespace lemp
{

using namespace std;

//linear indexing is used throughout
//linear index = z*num_x*num_y + y*num_x + x;

//convert linear index to vector subscript
Eigen::Vector3i ind2sub(int linear_index, Eigen::Vector3i dims){
    Eigen::Vector3i subs;
    subs[2] = linear_index/(dims[0]*dims[1]);
    linear_index = linear_index%(dims[0]*dims[1]);
    subs[1] = linear_index/dims[0];
    subs[0] = linear_index%dims[0];
    return subs;
}

//convert vector subscript to linear index
int sub2ind(Eigen::Vector3i subs, Eigen::Vector3i dims){
    return subs[2]*dims[0]*dims[1]+subs[1]*dims[0]+subs[0];
}


//get linear indices of all non-zero voxels in grid
vector<int> findIndexes(gridPtr band){
    vector<int> indexes_copy (band->dims[0]*band->dims[1]*band->dims[2], 0);
    int num_ind = 0;
    for(int i=0; i<band->dims[0]; i++){
        for(int j=0; j<band->dims[1]; j++){
            for(int k=0; k<band->dims[2]; k++){
                if(band->voxels[i][j][k]!=0.0){
                    Eigen::Vector3i pnt;
                    pnt[0]=i; pnt[1]=j; pnt[2]=k;
                    indexes_copy[num_ind] = sub2ind(pnt, band->dims);
                    num_ind++;
                }
            }
        }
    }
    //resize indexes
    vector<int> indexes(num_ind, 0);
    for(int i=0; i<num_ind; i++){
        indexes[i] = indexes_copy[i];
    }
    return indexes;
}

//create index map
gridPtr getIndexMap(gridPtr band, vector<int>& indexes){
    gridPtr map_ (new grid());
    map_->dims = band->dims;
    map_->t_ = band->t_;
    map_->voxels = allocGrid(map_->dims);
    //set all values to -1
    for(int i=0; i<map_->dims[0]; i++){
        for(int j=0; j<map_->dims[1]; j++){
            for(int k=0; k<map_->dims[2]; k++){
                map_->voxels[i][j][k]=-1.0;
            }
        }
    }
    //give i.d.'s to band location
    for(int i=0; i<indexes.size(); i++){
        Eigen::Vector3i pnt = ind2sub(indexes[i], band->dims);
        map_->voxels[pnt[0]][pnt[1]][pnt[2]]=(float)i;
    }
    return map_;
}


//make H matrix
SparseMatrixPtr getHMat(gridPtr tightBand, gridPtr indexMap){
    vector<int> indexes = findIndexes(tightBand);
    //set ntight
    int ntight = indexes.size();
    //get subscript vector for indexes
    Eigen::Vector3i subs[ntight];
    for(int i=0; i<ntight; i++){
        subs[i] = ind2sub(indexes[i], tightBand->dims);
    }

    //create Hi
    vector<int> Hi (ntight*9, 0);
    int index=0;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<ntight; k++){
                Hi[index] = k+ntight*i;
                index++;
            }
        }
    }

    //create Hj
    vector<int> Hj (ntight*9, 0);
    index=0;
    //add mid
    for(int i=0; i<ntight; i++){
        Hj[index] = indexMap->voxels[subs[i][0]][subs[i][1]][subs[i][2]];
        index++;
    }
    //add left
    for(int i=0; i<ntight; i++){
        Hj[index] = indexMap->voxels[subs[i][0]-1][subs[i][1]][subs[i][2]];
        index++;
    }
    //add right
    for(int i=0; i<ntight; i++){
        Hj[index] = indexMap->voxels[subs[i][0]+1][subs[i][1]][subs[i][2]];
        index++;
    }
    //add mid
    for(int i=0; i<ntight; i++){
        Hj[index] = indexMap->voxels[subs[i][0]][subs[i][1]][subs[i][2]];
        index++;
    }
    //add top
    for(int i=0; i<ntight; i++){
        Hj[index] = indexMap->voxels[subs[i][0]][subs[i][1]-1][subs[i][2]];
        index++;
    }
    //add bottom
    for(int i=0; i<ntight; i++){
        Hj[index] = indexMap->voxels[subs[i][0]][subs[i][1]+1][subs[i][2]];
        index++;
    }
    //add mid
    for(int i=0; i<ntight; i++){
        Hj[index] = indexMap->voxels[subs[i][0]][subs[i][1]][subs[i][2]];
        index++;
    }
    //add front
    for(int i=0; i<ntight; i++){
        Hj[index] = indexMap->voxels[subs[i][0]][subs[i][1]][subs[i][2]-1];
        index++;
    }
    //add back
    for(int i=0; i<ntight; i++){
        Hj[index] = indexMap->voxels[subs[i][0]][subs[i][1]][subs[i][2]+1];
        index++;
    }

    //create Hs
    vector<int> Hs (ntight*9, 1);
    index=0;
    for(int i=0; i<3; i++){
        for(int j=0; j<ntight; j++){
            Hs[index]=-2;
            index++;
        }
        for(int j=0; j<ntight*2; j++){
            Hs[index]=1;
            index++;
        }
    }


    //create initial H
    //get max value of Hj
    int maxj = 0;
    for(int i=0; i<Hj.size(); i++){
        if (Hj[i]>maxj) maxj=Hj[i];
    }
    //create triplets from Hi,Hj,Hs
    vector<Eigen::Triplet<float> > tripletList;
    tripletList.reserve(Hi.size());
    for(int i=0; i<Hi.size(); i++){
        tripletList.push_back(Eigen::Triplet<float>(Hi[i],Hj[i],(float)Hs[i]));
    }

    SparseMatrixPtr H (new Eigen::SparseMatrix<float>(3*ntight, maxj+1));
    H->setFromTriplets(tripletList.begin(), tripletList.end());

    *H = Eigen::SparseMatrix<float>(H->transpose())*(*H);
    return H;
}

//get lower bound vector
vector<float> getlb(gridPtr margin, gridPtr volume, vector<int>& indexes){
    //make a copy of margin
    //set values outside of volume to -1000
    gridPtr lbnd (new grid());
    lbnd->dims=margin->dims;
    lbnd->voxels=allocGrid(lbnd->dims);
    for(int i=0; i<lbnd->dims[0]; i++){
        for(int j=0; j<lbnd->dims[1]; j++){
            for(int k=0; k<lbnd->dims[2]; k++){
                if(volume->voxels[i][j][k]==0.0){
                    lbnd->voxels[i][j][k]=-1000.0;
                }
                else{
                    lbnd->voxels[i][j][k]=margin->voxels[i][j][k];
                }
            }
        }
    }

    //create lb vector using linear indexes
    vector<float> lb (indexes.size(), 0);
    for(int i=0; i<lb.size(); i++){
        Eigen::Vector3i pnt = ind2sub(indexes[i], lbnd->dims);
        lb[i] = lbnd->voxels[pnt[0]][pnt[1]][pnt[2]];
    }

    return lb;
}
//get upper bound vector
vector<float> getub(gridPtr margin, gridPtr volume, vector<int>& indexes){
    //make a copy of negative margin
    //set values inside of volume to 1000
    gridPtr ubnd (new grid());
    ubnd->dims=margin->dims;
    ubnd->voxels=allocGrid(ubnd->dims);
    for(int i=0; i<ubnd->dims[0]; i++){
        for(int j=0; j<ubnd->dims[1]; j++){
            for(int k=0; k<ubnd->dims[2]; k++){
                if(volume->voxels[i][j][k]==1.0){
                    ubnd->voxels[i][j][k]=1000.0;
                }
                else{
                    ubnd->voxels[i][j][k]=-margin->voxels[i][j][k];
                }
            }
        }
    }

    //create ub vector using linear indexes
    vector<float> ub (indexes.size(), 0);
    for(int i=0; i<ub.size(); i++){
        Eigen::Vector3i pnt = ind2sub(indexes[i], ubnd->dims);
        ub[i] = ubnd->voxels[pnt[0]][pnt[1]][pnt[2]];
    }

    return ub;
}

//prepare quadratic program arguments
qp_argsPtr primeQP(gridPtr volume, gridPtr margin, bandsPtr bnds){
    //create indexes of band
    vector<int> indexes = findIndexes(bnds->band);
    //create index map
    gridPtr indexMap = getIndexMap(bnds->band, indexes);
    //create H matrix
    SparseMatrixPtr H = getHMat(bnds->tight_band, indexMap);
    //H matrix has size nband x nband
    //create upper and lower bounds
    vector<float> lb_ = getlb(margin, volume, indexes);
    vector<float> ub_ = getub(margin, volume, indexes);
    //lb_ and ub_ have length nband

    //create x vector
    vector<float> x_ (H->rows(),0);
    for(int i=0; i<x_.size(); i++){
        //set values within upper and lower bounds
        if(x_[i]<lb_[i]) x_[i]=lb_[i];
        if(x_[i]>ub_[i]) x_[i]=ub_[i];
    }

    //get invdg
    Eigen::SparseVector<float> diag = (H->diagonal()).sparseView();
    vector<float> dg (diag.size(), 0);
    vector<float> invdg_ (diag.size(), 0);
    for (Eigen::SparseVector<float>::InnerIterator it(diag); it; ++it)
    {
        dg[it.index()] = (float)it.value();
        invdg_[it.index()] = 1.0/(float)it.value();
    }

    //create 0-diagonal matrix
    vector<Eigen::Triplet<float> > diag_trips;
    diag_trips.reserve(H->rows());
    for(int i=0; i<dg.size(); i++){
        if(dg[i]!=0.0){
            diag_trips.push_back(Eigen::Triplet<float>(i,i,dg[i]));
        }
    }
    SparseMatrixPtr d (new Eigen::SparseMatrix<float>(H->rows(), H->cols()));
    d->setFromTriplets(diag_trips.begin(), diag_trips.end());
    *H = ((*H)-(*d));
    H->prune(0,0);

    //create output struct
    qp_argsPtr out (new qp_args());
    out->R = H;
    out->invdg = invdg_;
    out->lb = lb_;
    out->ub = ub_;
    out->x = x_;
    out->iter = 10;

    cout<<"quadratic program ready"<<endl;
    return out;
}


void visualizeGrid(gridPtr grid){
    //convert imbedding function to point cloud for visualization
    //create point clouds from bands for visualization
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr pcl_grid (new pcl::PointCloud<pcl::PointXYZRGB>());
    for(int i=0; i<grid->dims[0]; i++){
        for(int j=0; j<grid->dims[1]; j++){
            for(int k=0; k<grid->dims[2]; k++){
                if(grid->voxels[i][j][k]>0){
                    pcl::PointXYZRGB pnt;
                    pnt.x=(float)i; pnt.y=(float)j; pnt.z=(float)k;
                    pnt.r=0;pnt.g=0;pnt.b=255;
                    pcl_grid->push_back(pnt);
                }
                else if(grid->voxels[i][j][k]<0){
                    pcl::PointXYZRGB pnt;
                    pnt.x=(float)i; pnt.y=(float)j; pnt.z=(float)k;
                    pnt.r=0;pnt.g=0;pnt.b=0;
                    pcl_grid->push_back(pnt);
                }
                else{
                    pcl::PointXYZRGB pnt;
                    pnt.x=(float)i; pnt.y=(float)j; pnt.z=(float)k;
                    pnt.r=0;pnt.g=255;pnt.b=0;
                    pcl_grid->push_back(pnt);
                }
            }
        }
    }

    //visualize
    //display in visualizor
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->setBackgroundColor (0, 0, 0);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(pcl_grid);
    viewer->addPointCloud<pcl::PointXYZRGB> (pcl_grid, rgb, "cloud");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "cloud");
    viewer->addCoordinateSystem (1.0);
    viewer->initCameraParameters ();
    while (!viewer->wasStopped ())
    {
        viewer->spinOnce (100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }

}

}
