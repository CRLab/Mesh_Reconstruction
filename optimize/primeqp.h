#ifndef PRIMEQP_H
#define PRIMEQP_H

#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <iostream>

#include "narrowBand.h"

namespace lemp
{

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

//convert linear index to vector subscript
Eigen::Vector3i ind2sub(int linear_index, Eigen::Vector3i dims);
//convert vector subscript to linear index
int sub2ind(Eigen::Vector3i subs, Eigen::Vector3i dims);

//get linear indices of all non-zero voxels in grid
vector<int> findIndexes(gridPtr band);

//create index map
gridPtr getIndexMap(gridPtr band, vector<int>& indexes);

//make H matrix
SparseMatrixPtr getHMat(gridPtr tightBand, gridPtr indexMap);

//get lower bound vector
vector<float> getlb(gridPtr margin, gridPtr volume, vector<int>& indexes);
//get upper bound vector
vector<float> getub(gridPtr margin, gridPtr volume, vector<int>& indexes);

//prepare quadratic program arguments
qp_argsPtr primeQP(gridPtr volume, gridPtr margin, bandsPtr bnds);

//visualize grid data
void visualizeGrid(gridPtr grid);

}


#endif
