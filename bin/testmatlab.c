#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSize mrows,mcols;
    char *name;
    double *x,*y;
    mxArray *rhs[1], *lhs[1];
    
    /* check for proper number of arguments */
    if(nrhs!=2) 
      mexErrMsgTxt("Two input required.");
    else if(nlhs > 1) 
      mexErrMsgTxt("One output arguments.");
    
    /* check to make sure the first input argument is a scalar */
    mrows = mxGetM(prhs[0]);
    mcols = mxGetN(prhs[0]);
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
          mrows*mcols!=1 ) {
        mexErrMsgTxt("Input x must be a scalar.");
    }
  
    /*  get the scalar input x */
    //x = mxGetScalar(prhs[0]);

    /* input must be a string */
    if ( mxIsChar(prhs[1]) != 1)
      mexErrMsgTxt("Input must be a string.");

    /* input must be a row vector */
    if (mxGetM(prhs[1])!=1)
      mexErrMsgTxt("Input must be a row vector.");
    
    /* copy the string data from prhs[0] into a C string input_ buf.    */
    name= mxArrayToString(prhs[1]);
    
    if(name == NULL) 
      mexErrMsgTxt("Could not convert input to string.");
    
    mexCallMATLAB(nlhs, plhs, 1, prhs, name);

    //plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);

    mxFree(name);
    return;
}

