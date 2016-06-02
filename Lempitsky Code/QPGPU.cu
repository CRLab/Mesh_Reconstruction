#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include "mex.h"
                     

texture<float> inTex;

__global__ void ProcessRow(float *out,
    int *ir, float *pr, float *invdg,
    float *lb, float *ub, int nRows, int nCols)
{
    int row = blockIdx.x*blockDim.x+threadIdx.x;

    if(row < nRows) {
        float res = 0;
        int i;

        for(i = 0; i < nCols; i++) 
        {
            int a = nRows*i+row;
            float val = pr[a];
            float x = tex1Dfetch(inTex, ir[a]);
            if(val)
                res +=  val*x;
        }

        res = (tex1Dfetch(inTex, row)-res*invdg[row])*0.5f;
        if(res < lb[row]) res = lb[row];
        if(res > ub[row]) res = ub[row];
        out[row] = res;
    }
}


//nvmex -f nvmexopts.bat QPGPU.cu -Ic:\cuda\include -Lc:\cuda\lib -lcudart

void mexFunction(
		 int		nlhs,
		 mxArray	*plhs[],
		 int		nrhs,
		 const mxArray	*prhs[]
		 )
{
    if (nrhs != 7) 
      mexErrMsgIdAndTxt("MATLAB:WrongNumberOfInputs", "7 input arguments required");
    if (nlhs != 1) 
      mexErrMsgIdAndTxt("MATLAB:WrongNumberOfOutputs", "1 output argument required");
    int *ir = mxGetIr(prhs[0]);
    int *jc = mxGetJc(prhs[0]);
    double *pr = mxGetPr(prhs[0]);
    int size = mxGetM(prhs[0]);
    
    float *invdg = (float *)mxGetData(prhs[1]);
    float *lb = (float *)mxGetData(prhs[2]);
    float *ub = (float *)mxGetData(prhs[3]);
    float *x0 = (float *)mxGetData(prhs[4]);
    int iter = int(*mxGetPr(prhs[5]));

    float *prSingle;
    int *ir_;
    int ncols;

    ncols = (int)*mxGetPr(prhs[6]);
    prSingle = new float[ncols*size];
    ir_ = new int[ncols*size];
    memset(prSingle, 0, sizeof(float)*ncols*size);
    memset(ir_, 0, sizeof(int)*ncols*size);
    for(int i = 0; i < size; i++)
        for(int j = jc[i], k = 0; j < jc[i+1]; j++, k++)
        {
            prSingle[k*size+i] = (float)pr[j]; 
            ir_[k*size+i] = ir[j];
        }
    

    int *irCu;
    float *prCu, *invdgCu, *lbCu, *ubCu, *inCu, *outCu;

    int sizeMat = size*ncols;

    cudaMalloc((void**)&irCu, sizeMat * sizeof(int));
    cudaMalloc((void**)&prCu, sizeMat * sizeof(float));  
    cudaMalloc((void**)&invdgCu, size * sizeof(float));  
    cudaMalloc((void**)&lbCu, size * sizeof(float));  
    cudaMalloc((void**)&ubCu, size * sizeof(float));  
    cudaMalloc((void**)&inCu, size * sizeof(float));  
    cudaMalloc((void**)&outCu, size * sizeof(float));  

    cudaMemcpy(irCu, ir_, sizeMat * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(prCu, prSingle, sizeMat * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(invdgCu, invdg, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(lbCu, lb, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(ubCu, ub, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(inCu, x0, size * sizeof(float), cudaMemcpyHostToDevice);

   
    int BLOCK_SIZE = 256;
    int NBLOCKS = (int)ceil(double(size)/BLOCK_SIZE);

    cudaBindTexture(0, inTex, inCu, size*sizeof(float));

    clock_t startTime = clock();

    for(int i = 0; i < iter; i++)
    {
        ProcessRow<<<NBLOCKS,BLOCK_SIZE>>>(outCu, irCu, prCu, invdgCu, lbCu, ubCu, size, ncols);                                            
        cudaMemcpy(inCu, outCu, size * sizeof(float), cudaMemcpyDeviceToDevice);
    }
    cudaThreadSynchronize();
    mexPrintf("Elapsed time for the iterations = %lf\n", (double(clock())-startTime)/CLOCKS_PER_SEC);    
        
    int dims[] = {size,1};
    plhs[0] = mxCreateNumericArray( 2, dims, mxSINGLE_CLASS, mxREAL );
    float *result = (float *)mxGetData(plhs[0]);
    cudaMemcpy(result, outCu, size * sizeof(float), cudaMemcpyDeviceToHost);
    
    delete ir_;
    delete prSingle;
    cudaFree(irCu);
    cudaFree(prCu);
    cudaFree(invdgCu);
    cudaFree(lbCu);
    cudaFree(ubCu);
    cudaFree(inCu);
    cudaFree(outCu);
}