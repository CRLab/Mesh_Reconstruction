
#include "primeqp.h"

using namespace std;

//make H matrix
SparseMatrixPtr getHMat(gridPtr tightBand, gridPtr indexMap){
    vector<int> indexes = findIndexes(tightBand);
    //set ntight
    int ntight = indexes.size();

    //get subscript vector for indexes
    Eigen::Vector3i subs[ntight];
    for(int i=0; i<ntight; i++){
        subs[i] = tightBand->ind2sub(indexes[i]);
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
        Hj[index] = (*indexMap)[subs[i][0]][subs[i][1]][subs[i][2]];
        index++;
    }
    //add left
    for(int i=0; i<ntight; i++){
        Hj[index] = (*indexMap)[subs[i][0]-1][subs[i][1]][subs[i][2]];
        index++;
    }
    //add right
    for(int i=0; i<ntight; i++){
        Hj[index] = (*indexMap)[subs[i][0]+1][subs[i][1]][subs[i][2]];
        index++;
    }
    //add mid
    for(int i=0; i<ntight; i++){
        Hj[index] = (*indexMap)[subs[i][0]][subs[i][1]][subs[i][2]];
        index++;
    }
    //add top
    for(int i=0; i<ntight; i++){
        Hj[index] = (*indexMap)[subs[i][0]][subs[i][1]-1][subs[i][2]];
        index++;
    }
    //add bottom
    for(int i=0; i<ntight; i++){
        Hj[index] = (*indexMap)[subs[i][0]][subs[i][1]+1][subs[i][2]];
        index++;
    }
    //add mid
    for(int i=0; i<ntight; i++){
        Hj[index] = (*indexMap)[subs[i][0]][subs[i][1]][subs[i][2]];
        index++;
    }
    //add front
    for(int i=0; i<ntight; i++){
        Hj[index] = (*indexMap)[subs[i][0]][subs[i][1]][subs[i][2]-1];
        index++;
    }
    //add back
    for(int i=0; i<ntight; i++){
        Hj[index] = (*indexMap)[subs[i][0]][subs[i][1]][subs[i][2]+1];
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
vector<float> getlb(gridPtr margin, gridPtr volume, const vector<int> &indexes){
    //make a copy of margin
    //set values outside of volume to -1000
    gridPtr lbnd (new grid(margin->dims, margin->scale, margin->shift, margin->pad));
    for(int i=0; i<lbnd->dims[0]; i++){
        for(int j=0; j<lbnd->dims[1]; j++){
            for(int k=0; k<lbnd->dims[2]; k++){
                if((*volume)[i][j][k]==0.0){
                    (*lbnd)[i][j][k]=-1000.0;
                }
                else{
                    (*lbnd)[i][j][k]=(*margin)[i][j][k];
                }
            }
        }
    }

    //create lb vector using linear indexes
    vector<float> lb (indexes.size(), 0);
    for(int i=0; i<lb.size(); i++){
        lb[i] = (*lbnd)(indexes[i]);
    }

    return lb;
}
//get upper bound vector
vector<float> getub(gridPtr margin, gridPtr volume, const vector<int> &indexes){
    //make a copy of negative margin
    //set values inside of volume to 1000
    gridPtr ubnd (new grid(margin->dims, margin->scale, margin->shift, margin->pad));
    for(int i=0; i<ubnd->dims[0]; i++){
        for(int j=0; j<ubnd->dims[1]; j++){
            for(int k=0; k<ubnd->dims[2]; k++){
                if((*volume)[i][j][k]==1.0){
                    (*ubnd)[i][j][k]=1000.0;
                }
                else{
                    (*ubnd)[i][j][k]=-(*margin)[i][j][k];
                }
            }
        }
    }

    //create ub vector using linear indexes
    vector<float> ub (indexes.size(), 0);
    for(int i=0; i<ub.size(); i++){
        ub[i] = (*ubnd)(indexes[i]);
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
    out->iter = 500;

    cout<<"quadratic program ready"<<endl;
    return out;
}

//************************************************************************************
//get feature index vector from feature map and index map
//sorted low to high
bool lowtohigh(int i, int j){
    return (i<j);
}
vector<int> getFeatureIndexes(gridPtr featureMap, gridPtr indexMap){
    vector<int> out;
    for(int i=0; i<featureMap->dims[0]; i++){
        for(int j=0; j<featureMap->dims[1]; j++){
            for(int k=0; k<featureMap->dims[2]; k++){
                if((*featureMap)[i][j][k]==1.0){
                    out.push_back((int)(*indexMap)[i][j][k]);
                }
            }
        }
    }
    sort(out.begin(), out.end(), lowtohigh);
    return out;
}
