#ifndef QP_CONF_H
#define QP_CONF_H

#include "primeqp_conf.h"

namespace conf
{

using namespace std;
using namespace conf;

//more formatting for rapid algorithm
//takes advantage of sparseness
vector<int> getIr(SparseMatrixPtr R);
vector<int> getJc(SparseMatrixPtr R);
vector<float> getPr(SparseMatrixPtr R);

//row increment
void doStep(int row, const vector<float> &in, vector<float> &out,
    const vector<int> &ir, const vector<int> &jc, const vector<float> &pr,
    const vector<float> &lb, const vector<float> &ub);

//quadratic programming optimization algorithm
vector<float> runQP(qp_argsPtr args);


//Function for computing weighted voxel grid for marching cubes
//takes as input a binary volume
gridPtr optimize(gridPtr confGrid, gridPtr volume);

}

#endif
