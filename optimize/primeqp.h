#ifndef PRIMEQP_H
#define PRIMEQP_H

#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include <Eigen/Eigen>
#include <Eigen/Sparse>

#include "narrowBand.h"

using namespace std;

//linear indexing is used throughout
//linear index = x*num_y*num_z + y*num_z + z

typedef boost::shared_ptr<Eigen::SparseMatrix<float> > SparseMatrixPtr;

struct qp_args{
    SparseMatrixPtr R;
    vector<float> invdg;
    vector<float> lb;
    vector<float> ub;
    vector<float> x;
    int iter;
};
typedef boost::shared_ptr<qp_args> qp_argsPtr;

//make H matrix
SparseMatrixPtr getHMat(gridPtr tightBand, gridPtr indexMap);

//get lower bound vector
vector<float> getlb(gridPtr margin, gridPtr volume, const vector<int>& indexes);
//get upper bound vector
vector<float> getub(gridPtr margin, gridPtr volume, const vector<int>& indexes);

//prepare quadratic program arguments
qp_argsPtr primeQP(gridPtr volume, gridPtr margin, bandsPtr bnds);

//************************************************************************************
//get feature index vector from feature map and index map
bool lowtohigh(int i, int j);
vector<int> getFeatureIndexes(gridPtr featureMap, gridPtr indexMap);

#endif
