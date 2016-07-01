
#include "quadprog.h"

namespace lemp
{

using namespace std;
using namespace lemp;

//more formatting for rapid algorithm
//takes advantage of sparseness
vector<int> getIr(SparseMatrixPtr R){
    vector<int> ir (R->nonZeros());
    int count=0;
    for (int k=0; k<R->outerSize(); ++k){
        for (Eigen::SparseMatrix<float>::InnerIterator it(*R,k); it; ++it){
            ir[count] = it.row();
            count++;
        }
    }
    return ir;

}
vector<int> getJc(SparseMatrixPtr R){
    vector<int> jc ((R->cols())+1);
    jc[0]=0;
    int sum=0;
    for (int k=0; k<R->outerSize(); ++k){
        for (Eigen::SparseMatrix<float>::InnerIterator it(*R,k); it; ++it){
            sum++;
        }
        jc[k+1]=sum;
    }
    return jc;
}
vector<float> getPr(SparseMatrixPtr R){
    vector<float> pr (R->nonZeros(), 0);
    int count=0;
    for (int k=0; k<R->outerSize(); ++k){
        for (Eigen::SparseMatrix<float>::InnerIterator it(*R,k); it; ++it){
            pr[count]=it.value();
            count++;
        }
    }
    return pr;
}

//row increment
void doStep(int row, vector<float>& in, vector<float>& out,
    vector<int>& ir, vector<int>& jc, vector<float>& pr, vector<float>& invdg,
    vector<float>& lb, vector<float>& ub)
{
    float res = 0;
    int start = jc[row];
    int end   = jc[row+1];
    int i;

    for(i = start; i < end; i++)
        res +=  pr[i]*in[ir[i]];

    res = (in[row]-res*invdg[row])*0.5;
    if(res < lb[row]) res = lb[row];
    if(res > ub[row]) res = ub[row];
    out[row] = res;
}

//quadratic programming optimization algorithm
vector<float> runQP(qp_argsPtr args){
    //get ir and jc
    vector<int> ir = getIr(args->R);
    vector<int> jc = getJc(args->R);
    vector<float> pr = getPr(args->R);

    int size = args->R->rows();

    vector<float>* buf1 = new vector<float>(args->x);
    vector<float>* buf2 = new vector<float>(size, 0);
    vector<float>* in;
    vector<float>* out;

    //perform algorithm
    for(int i = 0; i < args->iter; i++){
        if(i % 2){
            in = buf2;
            out = buf1;
        }
        else{
            in = buf1;
            out = buf2;
        }

        for(int r = 0; r < size; r++){
            doStep(r, *in, *out, ir, jc, pr, args->invdg, args->lb, args->ub);
        }
    }

    cout<<"quadratic program finished"<<endl;
    return *out;
}



//Function for computing weighted voxel grid for marching cubes
//takes as input a binary volume
gridPtr optimize(gridPtr volume){
    int BAND_SIZE=4.0;
    //prime quadratic programming arguments
    //prepare margin
    gridPtr margin = getsqrt(getsqdist(fastPerim(volume)));
    cout<<"margin calculated"<<endl;
    //prepare bands
    bandsPtr bnds = createBands(margin, BAND_SIZE);
    cout<<"bands created"<<endl;
    //get band point indexes
    vector<int> indexes = findIndexes(bnds->band);
    cout<<"indexes stored"<<endl;

    //prepare qp_args
    qp_argsPtr args = primeQP(volume, margin, bnds);

    //run quadratic programming
    vector<float> x = runQP(args);


    //prepare new voxel grid with imbedding function
    gridPtr F = copyGrid(volume);
    for(int i=0; i<F->dims[0]; i++){
        for(int j=0; j<F->dims[1]; j++){
            for(int k=0; k<F->dims[2]; k++){
                (*F)[i][j][k]=(2.0*(*F)[i][j][k]-1.0)*(BAND_SIZE+1.0);
            }
        }
    }
    for(int i=0; i<indexes.size(); i++){
        (*F)(indexes[i]) = x[i];
    }

    return F;
}

}


