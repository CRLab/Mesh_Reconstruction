#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{ 
   if(nrhs != 1)
        mexErrMsgTxt("1 input parameters required!");
   if(nlhs != 1)
        mexErrMsgTxt("1 output parameters required!");
   
   int nx = (mxGetDimensions(prhs[0]))[0];
   int ny = (mxGetDimensions(prhs[0]))[1];   
   int nz;
   if(mxGetNumberOfDimensions(prhs[0]) == 3)
	   nz = (mxGetDimensions(prhs[0]))[2];      
   else
       nz = 1;
   
   int nIter;
   unsigned char *im = (unsigned char *)mxGetData(prhs[0]); 
   
   
   int dims[] = {nx, ny, nz};
   plhs[0] =  mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
   unsigned char *out = (unsigned char *)mxGetData(plhs[0]);    
   
   int x,y,z,x_,y_,z_, addr,addr2;
   for(z = 1; z < nz-1; z++)
       for(y = 1; y < ny-1; y++)
           for(x = 1; x < nx-1; x++)
               for(z_=z-1; z_ <= z+1; z_++)
                   for(y_=y-1; y_ <= y+1; y_++)
                       for(x_=x-1; x_ <= x+1; x_++)
                       {
                           addr = z*nx*ny+y*nx+x;
                           addr2 = z_*nx*ny+y_*nx+x_;
                           if(im[addr] != im[addr2])
                                out[addr] = 1;
                       }
}       