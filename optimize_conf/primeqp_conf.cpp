
#include "primeqp_conf.h"

namespace conf
{

using namespace std;

//linear indexing is used throughout
//linear index = z*num_x*num_y + y*num_x + x;

//apply confidences to all voxels using gaussian distribution on distance field
gridPtr applyConfidence(gridPtr confGrid){
    gridPtr margin = copyGrid(confGrid);
    gridPtr margin_ind = copyGrid(confGrid);
    //set all volume points to 0 in margin
    for(int i=0; i<margin->dims[0]; i++){
        for(int j=0; j<margin->dims[1]; j++){
            for(int k=0; k<margin->dims[2]; k++){
                if((*margin)[i][j][k]>-1.0){
                    (*margin)[i][j][k]=0.0;
                }
                else{
                    (*margin)[i][j][k]=1.0;
                }
            }
        }
    }
    //get distance field
    margin = getsqrt(getsqdist(margin));
    //get distance field indexes
    margin_ind = getsqdist_index(margin);


    //apply confidences based on gaussian
    //mean=0, variance=2
    for(int i=0; i<margin->dims[0]; i++){
        for(int j=0; j<margin->dims[1]; j++){
            for(int k=0; k<margin->dims[2]; k++){
                //get confidence value of closest point
                float conf = (*confGrid)((int)(*margin_ind)[i][j][k]);
                //get new confidence value using gaussian
                float dist = (*margin)[i][j][k];
                float var = 1.0;
                float exp = -dist*dist/(2.0*var);
                float gauss = pow(2.71828f,exp);
                (*margin)[i][j][k] = conf*gauss;
                if((*margin)[i][j][k]<0.001) (*margin)[i][j][k]=0.0;
            }
        }
    }

    return margin;
}

//solve for inverse of sparse matrix
SparseMatrixPtr getInverse(SparseMatrixPtr in){
    Eigen::SparseMatrix<float> I(in->cols(),in->cols());
    I.setIdentity();

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float> > solver;
    solver.compute(*in);
    SparseMatrixPtr out (new Eigen::SparseMatrix<float> (solver.solve(I)));

    return out;
}

//create diagonal confidence matrix
SparseMatrixPtr getCMat(gridPtr confGrid, const vector<int> &indexes){
    //create triplet list
    vector<Eigen::Triplet<float> > tripletList;
    tripletList.reserve(indexes.size());
    for(int i=0; i<indexes.size(); i++){
        float conf = (*confGrid)(indexes[i]);
        tripletList.push_back(Eigen::Triplet<float>(i,i,conf));
    }
    //create sparse matrix
    SparseMatrixPtr C (new Eigen::SparseMatrix<float>(indexes.size(),indexes.size()));
    C->setFromTriplets(tripletList.begin(), tripletList.end());
    return C;
}

//get triplet list from sparse matrix
vector<Eigen::Triplet<float> > getTriplets(SparseMatrixPtr in){
    vector<Eigen::Triplet<float> > out;
    out.reserve(in->nonZeros());
    for(int k=0; k<in->outerSize(); k++){
        for(Eigen::SparseMatrix<float>::InnerIterator it(*in,k); it; ++it){
            out.push_back(Eigen::Triplet<float>(it.row(),it.col(),it.value()));
        }
    }
    return out;
}
//function for comparing triplets by row major order
bool compRow(const Eigen::Triplet<float> &in1, const Eigen::Triplet<float> &in2){
    if(in1.row()<in2.row()) return true;
    else if(in1.row()==in2.row()){
        return (in1.col()<in2.col());
    }
    else return false;
}
//reorder triplets by row major order
void sortByRow(vector<Eigen::Triplet<float> > &in){
    bool (*comp)(const Eigen::Triplet<float> &in1, const Eigen::Triplet<float> &in2);
    comp = &compRow;
    sort(in.begin(), in.end(), comp);
}
//function for multiplying sparse matrix by vector
vector<float> multiplyMatVec(SparseMatrixPtr mat, vector<float> &vec){
    if(mat->cols()!=vec.size()){
        cerr<<"Matrix and vector sizes must match"<<endl;
        return vec;
    }
    vector<float> out (vec.size(), 0.0f);
    vector<Eigen::Triplet<float> > trips = getTriplets(mat);
    sortByRow(trips);
    int row=0;
    float val=0.0f;
    for(int i=0; i<trips.size(); i++){
        if(trips[i].row()>row){
            out[row] = val;
            val=0.0f;
            row=trips[i].row();
        }
        val+=(trips[i].value()*vec[trips[i].col()]);
    }
    out[row] = val;
    return out;
}
//function for adding two vectors
vector<float> addVec(vector<float> &in1, const vector<float> &in2){
    if(in1.size()!=in2.size()){
        cerr<<"Vector sizes must match"<<endl;
        return in1;
    }
    vector<float> out (in1.size(), 0.0);
    for(int i=0; i<in1.size(); i++){
        out[i]=in1[i]+in2[i];
    }
    return out;
}
//function for subtracting two vectors
vector<float> subtractVec(vector<float> &in1, const vector<float> &in2){
    if(in1.size()!=in2.size()){
        return in1;
    }
    vector<float> out (in1.size(), 0.0);
    for(int i=0; i<in1.size(); i++){
        out[i]=in1[i]-in2[i];
    }
    return out;
}

//function for getting z vector
vector<float> getZVec(SparseMatrixPtr R, SparseMatrixPtr C, vector<float> &x_0){
    SparseMatrixPtr M (new Eigen::SparseMatrix<float>((*R)+(*C)));
    cout<<"computing inverse"<<endl;
    *M = (*getInverse(M))*(*C);
    cout<<"inverse computed"<<endl;
    M->prune(0,0);
    return multiplyMatVec(M, x_0);
}



//**********************************************************************

//create index map
gridPtr getIndexMap(gridPtr band, const vector<int> &indexes){
    gridPtr map_(new grid(band->dims, band->t_));
    //set all values to -1
    for(int i=0; i<map_->dims[0]; i++){
        for(int j=0; j<map_->dims[1]; j++){
            for(int k=0; k<map_->dims[2]; k++){
                (*map_)[i][j][k]=-1.0;
            }
        }
    }
    //give i.d.'s to band location
    for(int i=0; i<indexes.size(); i++){
        (*map_)(indexes[i])=(float)i;
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
    gridPtr lbnd (new grid(margin->dims, margin->t_));
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
    gridPtr ubnd (new grid(margin->dims, margin->t_));
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
qp_argsPtr primeQP(gridPtr confGrid, gridPtr volume, gridPtr margin, bandsPtr bnds){
    //create indexes of band
    vector<int> indexes = findIndexes(bnds->band);
    //create index map
    gridPtr indexMap = getIndexMap(bnds->band, indexes);
    //create H matrix
    SparseMatrixPtr H = getHMat(bnds->tight_band, indexMap);
    //H matrix has size nband x nband

    //create C matrix
    SparseMatrixPtr C = getCMat(applyConfidence(confGrid), indexes);

    //create upper and lower bounds
    vector<float> lb_ = getlb(margin, volume, indexes);
    vector<float> ub_ = getub(margin, volume, indexes);
    //lb_ and ub_ have length nband

    //get invdg
    Eigen::SparseVector<float> diag = (H->diagonal()).sparseView();
    vector<Eigen::Triplet<float> > diag_trips;
    diag_trips.reserve(H->rows());
    vector<Eigen::Triplet<float> > invdiag_trips;
    invdiag_trips.reserve(H->rows());
    for (Eigen::SparseVector<float>::InnerIterator it(diag); it; ++it)
    {
        diag_trips.push_back(Eigen::Triplet<float>(it.index(),it.index(),(float)it.value()));
        invdiag_trips.push_back(Eigen::Triplet<float>(it.index(),it.index(),1.0/(float)it.value()));
    }
    SparseMatrixPtr dg (new Eigen::SparseMatrix<float>(H->rows(), H->cols()));
    dg->setFromTriplets(diag_trips.begin(), diag_trips.end());
    SparseMatrixPtr invdg (new Eigen::SparseMatrix<float>(H->rows(), H->cols()));
    invdg->setFromTriplets(invdiag_trips.begin(), invdiag_trips.end());
    *H = ((*H)-(*dg));
    H->prune(0,0);

    //normalize H
    *H = ((*H)*(*invdg));
    H->prune(0,0);

    //create x vector
    vector<float> x_ (H->rows(),0);
    for(int i=0; i<x_.size(); i++){
        //set values within upper and lower bounds
        if(x_[i]<lb_[i]) x_[i]=lb_[i];
        if(x_[i]>ub_[i]) x_[i]=ub_[i];
    }

    //create z vector
    vector<float> z_ = getZVec(H, C, x_);

    vector<float> x_in = subtractVec(x_, z_);

    //create M matrix
    *H = ((*H)+(*C));
    H->prune(0,0);

    //create output struct
    qp_argsPtr out (new qp_args());
    out->M = H;
    out->lb = lb_;
    out->ub = ub_;
    out->x = x_in;
    out->z = z_;
    //out->z = x_;
    out->iter = 500;

    cout<<"quadratic program ready"<<endl;
    return out;
}

}
