
#include <iostream>

#include "narrowBand.h"
#include "primeqp.h"

using namespace std;

//linear indexing is used throughout
//linear index = x*num_y*num_z + y*num_z + z

struct smartVectorInt{
    int* vec;
    int dim;
    ~smartVectorInt(){
        delete [] vec;
    }
};
struct smartVectorFloat{
    float* vec;
    int dim;
    ~smartVectorFloat(){
        delete [] vec;
    }
};
typedef boost::shared_ptr<smartVectorInt> smartVectorIntPtr;
typedef boost::shared_ptr<smartVectorFloat> smartVectorFloatPtr;

struct smartArrayInt{
    int** arr;
    int dims[2];
    ~smartArrayInt(){
        for(int i=0; i<dims[0]; i++){
            delete [] arr[i];
        }
        delete [] arr;
    }
};
struct smartArrayFloat{
    float** arr;
    int dims[2];
    ~smartArrayFloat(){
        for(int i=0; i<dims[0]; i++){
            delete [] arr[i];
        }
        delete [] arr;
    }
};
typedef boost::shared_ptr<smartArrayInt> smartArrayIntPtr;
typedef boost::shared_ptr<smartArrayFloat> smartArrayFloatPtr;

struct qp_args{
    smartArrayIntPtr R;
    smartVectorFloatPtr invdg;
    smartVectorFloatPtr lb;
    smartVectorFloatPtr ub;
    smartVectorFloatPtr x;
    float omega;

};
typedef boost::shared_ptr<qp_args> qp_argsPtr;

//transpose and mutliplication of smart arrays
smartArrayIntPtr transpose(smartArrayIntPtr in){
    smartArrayIntPtr out;
    out->dims[0]=in->dims[1];
    out->dims[1]=in->dims[0];
    out->arr = new int*[out->dims[0]];
    for(int i=0; i<out->dims[0]; i++){
        out->arr[i] = new int[out->dims[1]];
    }
    for(int i=0; i<in->dims[0]; i++){
        for(int j=0; j<in->dims[1]; j++){
            out->arr[j][i] = in->arr[i][j];
        }
    }
    return out;
}
smartArrayFloatPtr transpose(smartArrayFloatPtr in){
    smartArrayFloatPtr out;
    out->dims[0]=in->dims[1];
    out->dims[1]=in->dims[0];
    out->arr = new float*[out->dims[0]];
    for(int i=0; i<out->dims[0]; i++){
        out->arr[i] = new float[out->dims[1]];
    }
    for(int i=0; i<in->dims[0]; i++){
        for(int j=0; j<in->dims[1]; j++){
            out->arr[j][i] = in->arr[i][j];
        }
    }
    return out;
}
smartArrayIntPtr multiply(smartArrayIntPtr in1, smartArrayIntPtr in2){
    smartArrayIntPtr out;
    out->dims[0]=in1->dims[0];
    out->dims[1]=in2->dims[1];
    out->arr = new int*[out->dims[0]];
    for(int i=0; i<out->dims[0]; i++){
        out->arr[i] = new int[out->dims[1]];
    }
    for(int i=0; i<out->dims[0]; i++){
        for(int j=0; j<out->dims[1]; j++){
            //loop through row of in1 and col of in2
            int sum=0;
            for(int ind=0; ind<in1->dims[1]; ind++){
                sum+=(in1->arr[i][ind]*in2->arr[ind][j]);
            }
            out->arr[i][j]=sum;
        }
    }
    return out;
}
smartArrayFloatPtr multiply(smartArrayFloatPtr in1, smartArrayFloatPtr in2){
    smartArrayFloatPtr out;
    out->dims[0]=in1->dims[0];
    out->dims[1]=in2->dims[1];
    out->arr = new float*[out->dims[0]];
    for(int i=0; i<out->dims[0]; i++){
        out->arr[i] = new float[out->dims[1]];
    }
    for(int i=0; i<out->dims[0]; i++){
        for(int j=0; j<out->dims[1]; j++){
            //loop through row of in1 and col of in2
            float sum=0.0;
            for(int ind=0; ind<in1->dims[1]; ind++){
                sum+=(in1->arr[i][ind]*in2->arr[ind][j]);
            }
            out->arr[i][j]=sum;
        }
    }
    return out;
}


//convert linear index to vector subscript
Eigen::Vector3i ind2sub(int linear_index, Eigen::Vector3i dims){
    Eigen::Vector3i subs;
    subs[0] = linear_index/(dims[1]*dims[2]);
    linear_index = linear_index%(dims[1]*dims[2]);
    subs[1] = linear_index/dims[2];
    subs[2] = linear_index%dims[2];
    return subs;
}

//convert vector subscript to linear index
int sub2ind(Eigen::Vector3i subs, Eigen::Vector3i dims){
    return subs[0]*dims[1]*dims[2]+subs[1]*dims[2]+subs[2];
}






//get linear indices of all non-zero voxels in grid
smartVectorIntPtr findIndexes(gridPtr band){
    smartVectorIntPtr indexes (new smartVectorInt());
    indexes->vec = new int[band->dims[0]*band->dims[1]*band->dims[2]];
    int num_ind = 0;
    for(int i=0; i<band->dims[0]; i++){
        for(int j=0; j<band->dims[1]; j++){
            for(int k=0; k<band->dims[2]; k++){
                if(band->voxels[i][j][k]!=0.0){
                    Eigen::Vector3i pnt;
                    pnt[0]=i; pnt[1]=j; pnt[2]=k;
                    indexes->vec[num_ind++] = sub2ind(pnt, band->dims);
                }
            }
        }
    }
    //resize indexes
    int* indexes_= new int[num_ind];
    for(int i=0; i<num_ind; i++){
        indexes_[i] = indexes->vec[i];
    }
    delete [] indexes->vec;
    indexes->vec = indexes_;

    return indexes;
}

//create index map
gridPtr getIndexMap(gridPtr band, smartVectorIntPtr indexes){
    gridPtr map_ (new grid());
    map_->dims = band->dims;
    map_->voxels = allocGrid(map_->dims);
    //set all values to 0
    for(int i=0; i<map_->dims[0]; i++){
        for(int j=0; j<map_->dims[1]; j++){
            for(int k=0; k<map_->dims[2]; k++){
                map_->voxels[i][j][k]=-1.0;
            }
        }
    }
    //give i.d.'s to band location
    for(int i=0; i<indexes->dim; i++){
        Eigen::Vector3i pnt = ind2sub(indexes->vec[i], band->dims);
        map_->voxels[pnt[0]][pnt[1]][pnt[2]]=(float)i;
    }
    return map_;
}


//make H matrix
smartArrayIntPtr getHMat(gridPtr tightBand, gridPtr indexMap){
    smartVectorIntPtr indexes = findIndexes(tightBand);
    //set ntight
    int ntight = indexes->dim;
    //get subscript vector for indexes
    Eigen::Vector3i subs[ntight];
    for(int i=0; i<ntight; i++){
        subs[i] = ind2sub(indexes->vec[i], tightBand->dims);
    }
    //create Hi
    smartVectorIntPtr Hi;
    Hi->dim = ntight*9;
    Hi->vec = new int[Hi->dim];
    {
        int index=0;
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                for(int k=0; k<ntight; k++){
                    Hi->vec[index] = k*(i+1)+ntight*i;
                    index++;
                }
            }
        }
    }
    //create Hj
    smartVectorIntPtr Hj;
    Hj->dim = ntight*9;
    Hj->vec = new int[Hj->dim];
    {
        int index=0;
        //add mid
        for(int i=0; i<ntight; i++){
            Hj->vec[index] = indexMap->voxels[subs[i][0]][subs[i][1]][subs[i][2]];
            index++;
        }
        //add left
        for(int i=0; i<ntight; i++){
            Hj->vec[index] = indexMap->voxels[subs[i][0]-1][subs[i][1]][subs[i][2]];
            index++;
        }
        //add right
        for(int i=0; i<ntight; i++){
            Hj->vec[index] = indexMap->voxels[subs[i][0]+1][subs[i][1]][subs[i][2]];
            index++;
        }
        //add mid
        for(int i=0; i<ntight; i++){
            Hj->vec[index] = indexMap->voxels[subs[i][0]][subs[i][1]][subs[i][2]];
            index++;
        }
        //add top
        for(int i=0; i<ntight; i++){
            Hj->vec[index] = indexMap->voxels[subs[i][0]][subs[i][1]-1][subs[i][2]];
            index++;
        }
        //add bottom
        for(int i=0; i<ntight; i++){
            Hj->vec[index] = indexMap->voxels[subs[i][0]][subs[i][1]+1][subs[i][2]];
            index++;
        }
        //add mid
        for(int i=0; i<ntight; i++){
            Hj->vec[index] = indexMap->voxels[subs[i][0]][subs[i][1]][subs[i][2]];
            index++;
        }
        //add front
        for(int i=0; i<ntight; i++){
            Hj->vec[index] = indexMap->voxels[subs[i][0]][subs[i][1]][subs[i][2]-1];
            index++;
        }
        //add back
        for(int i=0; i<ntight; i++){
            Hj->vec[index] = indexMap->voxels[subs[i][0]][subs[i][1]][subs[i][2]+1];
            index++;
        }
    }
    //create Hs
    smartVectorIntPtr Hs;
    Hs->dim = ntight*9;
    Hs->vec = new int[Hs->dim];
    {
        int index=0;
        for(int i=0; i<3; i++){
            for(int j=0; j<ntight; j++){
                Hs->vec[index]=-2;
                index++;
            }
            for(int j=0; j<ntight; j++){
                Hs->vec[index]=1;
                index++;
            }
            for(int j=0; j<ntight; j++){
                Hs->vec[index]=1;
                index++;
            }
        }
    }
    //create initial H
    //get max value of Hj
    int maxj = 0;
    for(int i=0; i<Hj->dim; i++){
        if (Hj->vec[i]>maxj) maxj=Hj->vec[i];
    }
    smartArrayIntPtr H;
    H->dims[0] = 3*ntight;
    H->dims[1] = maxj;
    H->arr = new int*[H->dims[0]];
    for(int i=0; i<H->dims[0]; i++){
        H->arr[i] = new int[H->dims[1]];
    }
    //set all values to 0
    for(int i=0; i<H->dims[0]; i++){
        for(int j=0; j<H->dims[1]; j++){
            H->arr[i][j]=0;
        }
    }
    //set all nonzero values
    for(int i=0; i<ntight*9; i++){
        H->arr[Hi->vec[i]][Hj->vec[i]]=Hs->vec[i];
    }
    Hi->~smartVectorInt();
    Hj->~smartVectorInt();
    Hs->~smartVectorInt();

    //make H square hermitian
    H = multiply(transpose(H), H);
    return H;
}

//get lower bound vector
smartVectorFloatPtr getlb(gridPtr margin, gridPtr volume, smartVectorIntPtr indexes){
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
    smartVectorFloatPtr lb (new smartVectorFloat());
    lb->dim=indexes->dim;
    lb->vec=new float[lb->dim];
    for(int i=0; i<lb->dim; i++){
        Eigen::Vector3i pnt = ind2sub(indexes->vec[i], lbnd->dims);
        lb->vec[i] = lbnd->voxels[pnt[0]][pnt[1]][pnt[2]];
    }

    return lb;
}
//get upper bound vector
smartVectorFloatPtr getub(gridPtr margin, gridPtr volume, smartVectorIntPtr indexes){
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
                    ubnd->voxels[i][j][k]=0.0-margin->voxels[i][j][k];
                }
            }
        }
    }

    //create lb vector using linear indexes
    smartVectorFloatPtr ub (new smartVectorFloat());
    ub->dim=indexes->dim;
    ub->vec=new float[ub->dim];
    for(int i=0; i<ub->dim; i++){
        Eigen::Vector3i pnt = ind2sub(indexes->vec[i], ubnd->dims);
        ub->vec[i] = ubnd->voxels[pnt[0]][pnt[1]][pnt[2]];
    }

    return ub;
}

//prepare quadratic program arguments
qp_argsPtr primeQP(gridPtr volume, gridPtr margin, bandsPtr bnds){
    //create indexes of band
    smartVectorIntPtr indexes = findIndexes(bnds->band);
    //create index map
    gridPtr indexMap = getIndexMap(bnds->band, indexes);
    //create H matrix
    smartArrayIntPtr H_ = getHMat(bnds->tight_band, indexMap);
    //H matrix has size nband x nband
    //create upper and lower bounds
    smartVectorFloatPtr lb_ = getlb(margin, volume, indexes);
    smartVectorFloatPtr ub_ = getub(margin, volume, indexes);
    //lb_ and ub_ have length nband
    indexes->~smartVectorInt();
    indexMap->~grid();

    //create x vector
    smartVectorFloatPtr x_;
    x_->dim = H_->dims[0];
    x_->vec = new float[x_->dim];
    for(int i=0; i<x_->dim; i++){
        x_->vec[i]=0.0;
        //set values within upper and lower bounds
        if(x_->vec[i]<lb_->vec[i]) x_->vec[i]=lb_->vec[i];
        if(x_->vec[i]>ub_->vec[i]) x_->vec[i]=ub_->vec[i];
    }

    //create invdg
    smartVectorFloatPtr invdg_;
    invdg_->dim = H_->dims[0];
    invdg_->vec = new float[invdg_->dim];
    for(int i=0; i<invdg_->dim; i++){
        invdg_->vec[i] = 1.0/(float)H_->arr[i][i];
    }
    //create R matrix: H with 0's on diagonal
    smartArrayIntPtr R_ = H_;
    for(int i=0; i<R_->dims[0]; i++){
        R_->arr[i][i]=0;
    }

    //create output struct
    qp_argsPtr out (new qp_args());
    out->R = R_;
    out->invdg = invdg_;
    out->lb = lb_;
    out->ub = ub_;
    out->x = x_;
    out->omega = 1000.0;

    return out;
}













