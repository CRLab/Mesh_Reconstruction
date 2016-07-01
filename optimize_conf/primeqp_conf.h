#ifndef PRIMEQP_CONF_H
#define PRIMEQP_CONF_H

#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include "narrowBand.h"

namespace conf
{

using namespace std;

//linear indexing is used throughout
//linear index = x*num_y*num_z + y*num_z + z

typedef boost::shared_ptr<Eigen::SparseMatrix<float> > SparseMatrixPtr;

//apply confidences to all voxels using gaussian distribution on distance field
gridPtr applyConfidence(gridPtr confGrid);

//returns linear indexes of closes 0-value voxel
gridPtr getsqdist_index(gridPtr volume_grid);


//solve for inverse of sparse matrix
SparseMatrixPtr getInverse(SparseMatrixPtr in);

//create diagonal confidence matrix
SparseMatrixPtr getCMat(gridPtr confGrid, vector<int>& indexes);

//get triplet list from sparse matrix
vector<Eigen::Triplet<float> > getTriplets(SparseMatrixPtr in);
//function for comparing triplets by row major order
bool compRow(Eigen::Triplet<float> in1, Eigen::Triplet<float> in2);
//reorder triplets by row major order
void sortByRow(vector<Eigen::Triplet<float> >& in);
//function for multiplying sparse matrix by vector
vector<float> multiplyMatVec(SparseMatrixPtr mat, vector<float>& vec);
//function for adding two vectors
vector<float> addVec(vector<float>& in1, vector<float>& in2);
//function for subtracting two vectors
vector<float> subtractVec(vector<float>& in1, vector<float>& in2);

//function for getting z vector
vector<float> getZVec(SparseMatrixPtr R, SparseMatrixPtr C, vector<float> x_0);

struct qp_args{
    SparseMatrixPtr M;
    vector<float> lb;
    vector<float> ub;
    vector<float> x;
    vector<float> z;
    int iter;
};
typedef boost::shared_ptr<qp_args> qp_argsPtr;

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
qp_argsPtr primeQP(gridPtr confGrid, gridPtr volume, gridPtr margin, bandsPtr bnds);

//visualize grid data
void visualizeGrid(gridPtr grid);

}


#endif
