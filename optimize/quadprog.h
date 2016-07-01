#ifndef QUADPROG_H
#define QUADPROG_H

#include <math.h>

#include "primeqp.h"

namespace lemp
{

using namespace std;
using namespace lemp;

//more formatting for rapid algorithm
//takes advantage of sparseness
vector<int> getIr(SparseMatrixPtr R);
vector<int> getJc(SparseMatrixPtr R);
vector<float> getPr(SparseMatrixPtr R);

//row increment
void doStep(int row, vector<float>& in, vector<float>& out,
    vector<int>& ir, vector<int>& jc, vector<float>& pr, vector<float>& invdg,
    vector<float>& lb, vector<float>& ub);

//quadratic programming optimization algorithm
vector<float> runQP(qp_argsPtr args);


//Function for computing weighted voxel grid for marching cubes
//takes as input a binary volume
gridPtr optimize(gridPtr volume);

}

#endif
