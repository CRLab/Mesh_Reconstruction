#ifndef QUADPROG_H
#define QUADPROG_H

#include <math.h>

#include "primeqp.h"

using namespace std;

//more formatting for rapid algorithm
//takes advantage of sparseness
vector<int> getIr(SparseMatrixPtr R);
vector<int> getJc(SparseMatrixPtr R);
vector<float> getPr(SparseMatrixPtr R);

//row increment
void doStep(int row, const vector<float>& in, vector<float>& out,
    const vector<int>& ir, const vector<int>& jc, const vector<float>& pr, const vector<float>& invdg,
    const vector<float>& lb, const vector<float>& ub);

//quadratic programming optimization algorithm
vector<float> runQP(qp_argsPtr args);


//Function for computing weighted voxel grid for marching cubes
//takes as input a binary volume
gridPtr optimize(gridPtr volume);

#endif
