
#include "quadprog.h"
#include <cuda.h>

using namespace std;

#define GPU_CHECKERROR( err ) (gpuCheckError( err, __FILE__, __LINE__ ))
static void gpuCheckError( cudaError_t err,
                           const char *file,
                           int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}

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

__global__ void doStepDevice(
    float *in,
    float *out,
    int* ir,
    int* jc,
    float* pr,
    float* invdg,
    float* lb,
    float* ub,
    int nRows
)
{
    int row = blockIdx.x*blockDim.x+threadIdx.x;
    if(row < nRows) {

        float res = 0;
        int start = jc[row];
        int end = jc[row + 1];
        int i;

        for (i = start; i < end; i++)
            res += pr[i] * in[ir[i]];

        res = (in[row] - res * invdg[row]) * 0.5;
        if (res < lb[row]) res = lb[row];
        if (res > ub[row]) res = ub[row];
        out[row] = res;
    }
}

//row increment
void doStep(int row, const vector<float>& in, vector<float>& out,
    const vector<int>& ir, const vector<int>& jc, const vector<float>& pr, const vector<float>& invdg,
    const vector<float>& lb, const vector<float>& ub)
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

vector<float> runQPGPU(qp_argsPtr args,
                       const vector<int> &featureIndexes,
                       const bool USING_FEATURES)
{
    //get ir, jc, and pr
    vector<int> ir = getIr(args->R);
    vector<int> jc = getJc(args->R);
    vector<float> pr = getPr(args->R);

    int size = args->R->rows();

    vector<float>* out = new vector<float>(size, 0);

    int *irCu, *jcCu;
    float *prCu, *invdgCu, *lbCu, *ubCu, *inCu, *outCu;

    cout << "Creating memory spaces on GPU" << endl;
    GPU_CHECKERROR(cudaMalloc((void**)&inCu,    size * sizeof(float)));
    GPU_CHECKERROR(cudaMalloc((void**)&outCu,   size * sizeof(float)));
    GPU_CHECKERROR(cudaMalloc((void**)&irCu,    ir.size() * sizeof(int)));
    GPU_CHECKERROR(cudaMalloc((void**)&jcCu,    jc.size() * sizeof(int)));
    GPU_CHECKERROR(cudaMalloc((void**)&prCu,    pr.size() * sizeof(float)));
    GPU_CHECKERROR(cudaMalloc((void**)&invdgCu, args->invdg.size() * sizeof(float)));
    GPU_CHECKERROR(cudaMalloc((void**)&lbCu,    args->lb.size() * sizeof(float)));
    GPU_CHECKERROR(cudaMalloc((void**)&ubCu,    args->ub.size() * sizeof(float)));

    cout << "Copying memory to GPU" << endl;
    GPU_CHECKERROR(cudaMemcpy(inCu, &(args->x)[0],        size * sizeof(float), cudaMemcpyHostToDevice));
    GPU_CHECKERROR(cudaMemcpy(irCu, &ir[0],               ir.size() * sizeof(int), cudaMemcpyHostToDevice));
    GPU_CHECKERROR(cudaMemcpy(jcCu, &jc[0],               jc.size() * sizeof(int), cudaMemcpyHostToDevice));
    GPU_CHECKERROR(cudaMemcpy(prCu, &pr[0],               pr.size() * sizeof(float), cudaMemcpyHostToDevice));
    GPU_CHECKERROR(cudaMemcpy(invdgCu, &(args->invdg)[0], args->invdg.size() * sizeof(float), cudaMemcpyHostToDevice));
    GPU_CHECKERROR(cudaMemcpy(lbCu, &(args->lb)[0],       args->lb.size() * sizeof(float), cudaMemcpyHostToDevice));
    GPU_CHECKERROR(cudaMemcpy(ubCu, &(args->ub)[0],       args->ub.size() * sizeof(float), cudaMemcpyHostToDevice));

    int BLOCK_SIZE = 256;
    int NBLOCKS = (int)ceil(double(size)/BLOCK_SIZE);

    //perform algorithm
    cout << "Perform GPU algorithm" << endl;
    for(int i = 0; i < args->iter; i++){
        if(i > 0) {
            float* swap = inCu;
            inCu = outCu;
            outCu = swap;
        }
        doStepDevice<<<NBLOCKS,BLOCK_SIZE>>>(inCu, outCu, irCu, jcCu, prCu, invdgCu, lbCu, ubCu, size);
//        cudaMemcpy(inCu, outCu, size * sizeof(float), cudaMemcpyDeviceToDevice);
    }

    cudaThreadSynchronize();

    cudaMemcpy(&((*out)[0]), outCu, size * sizeof(float), cudaMemcpyDeviceToHost);

    cudaThreadSynchronize();

    cout << "Free all GPU Memory" << endl;
    cudaFree(irCu);
    cudaFree(jcCu);
    cudaFree(prCu);
    cudaFree(invdgCu);
    cudaFree(lbCu);
    cudaFree(ubCu);
    cudaFree(inCu);
    cudaFree(outCu);

    if(USING_FEATURES){
        //reset values at feature points
        int count=0;
        for(int i=0; i<args->x.size(); i++){
            if(i==featureIndexes[count]){
                (*out)[i]=0.2; //<--------------------------------------could be related to band size
                count++;
            }
        }
    }

    cout<<"quadratic program finished"<<endl;
    return *out;
}

//quadratic programming optimization algorithm
vector<float> runQP(qp_argsPtr args, const vector<int> &featureIndexes, const bool USING_FEATURES){
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
    if(USING_FEATURES){
        //reset values at feature points
        int count=0;
        for(int i=0; i<args->x.size(); i++){
            if(i==featureIndexes[count]){
                (*out)[i]=0.2; //<--------------------------------------could be related to band size
                count++;
            }
        }
    }

    cout<<"quadratic program finished"<<endl;
    return *out;
}



//Function for computing weighted voxel grid for marching cubes
//takes as input a binary volume
gridPtr optimize(gridPtr volume, gridPtr featureMap, const bool USING_FEATURES, const bool USING_CUDA){
    int BAND_SIZE=4.0;
    //prime quadratic programming arguments
    //prepare margin
    gridPtr margin;
    if(USING_CUDA)
        margin = dfield_gpu(volume);
    else
        margin = getsqrt(getsqdist(fastPerim(volume)));
    cout<<"margin calculated"<<endl;
    //prepare bands
    bandsPtr bnds = createBands(margin, BAND_SIZE);
    cout<<"bands created"<<endl;
    //get band point indexes
    vector<int> indexes = findIndexes(bnds->band);
    //create index map
    gridPtr indexMap = getIndexMap(bnds->band, indexes);
    vector<int> featureIndexes = getFeatureIndexes(featureMap, indexMap);
    cout<<"indexes stored"<<endl;

    //prepare qp_args
    qp_argsPtr args = primeQP(volume, margin, bnds);

    //run quadratic programming
    vector<float> x;
    if(USING_CUDA)
        x = runQPGPU(args, featureIndexes, USING_FEATURES);
    else
        x = runQP(args, featureIndexes, USING_FEATURES);

    //prepare new voxel grid with embedding function
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
