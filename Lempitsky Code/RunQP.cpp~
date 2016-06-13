#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include "mex.h"


void doStep(int row, double *in, double *out,
    int *ir, int *jc, double *pr, double *invdg,
    double *lb, double *ub)
{
    double res = 0;
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

/*void ProcessRow(int row, double *in, double *out,
    int *ir, int *jc, double *pr, double *invdg,
    double *lb, double *ub)
{
    double res = 0;
    int start = jc[row];
    int end   = jc[row+1];
    int i;
    
    for(i = start; i < end; i++) 
        res +=  pr[i]*in[ir[i]];
        
    res = (in[row]-res*invdg[row])*0.5;
    if(res < lb[row]) res = lb[row];
    if(res > ub[row]) res = ub[row];
    out[row] = res;
}*/


void mexFunction(
		 int		nlhs,
		 mxArray	*plhs[],
		 int		nrhs,
		 const mxArray	*prhs[]
		 )
{
    float    *data, *colors, *val, *noshare, *verbose;
    const mwSize *dims;
    mwSize	   ndims;
    mwSize  xSize, ySize, zSize;
    float	*vOut, *fOut, *cOut;
    
    
    /* Check for proper number of arguments */
    
    if (nrhs != 6) 
      mexErrMsgIdAndTxt("MATLAB:WrongNumberOfInputs",
                        "6 input arguments required");
    else if (nlhs > 1) 
      mexErrMsgIdAndTxt("MATLAB:WrongNumberOfOutputs",
                        "1 output argument required");
    int *ir = mxGetIr(prhs[0]);
    int *jc = mxGetJc(prhs[0]);
    double *pr = mxGetPr(prhs[0]);
    
    double *invdg = mxGetPr(prhs[1]);
    double *lb = mxGetPr(prhs[2]);
    double *ub = mxGetPr(prhs[3]);
    double *x0 = mxGetPr(prhs[4]);
    int iter = int(*mxGetPr(prhs[5]));
    
    
    int size = mxGetM(prhs[0]);
    
    double *buf1 = new double[size];
    double *buf2 = new double[size];    
    double *in, *out;
    
    memcpy(buf1, x0, sizeof(double)*size);
    
    for(int i = 0; i < iter; i++)
    {
        if(i % 2)        {
            in = buf2;
            out = buf1;
        }
        else        {
            in = buf1;
            out = buf2;
        }   
        
        for(int r = 0; r < size; r++)
            doStep(r, in, out, ir, jc, pr, invdg, lb, ub);                       
    }
    
    plhs[0] = mxCreateDoubleMatrix( size, 1, mxREAL );
    double *result = mxGetPr(plhs[0]);
    memcpy(result, out, sizeof(double)*size);
    
    delete buf1;
    delete buf2;
}