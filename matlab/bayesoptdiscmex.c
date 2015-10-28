/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU Affero General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Affero General Public License for more details.

   You should have received a copy of the GNU Affero General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/

#include "bayesoptextras.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  double *xptr;
  mxArray *xopt;
  const mxArray *func_name, *params;
  user_function_data udata;
  size_t nDim,nPoints;
  double* xset;
  bopt_params parameters;
  double fmin = 0.0;
  int error_code;
    
  /* Check correct number of parameters */
  CHECK0(nlhs != 2 || nrhs != 3, "wrong number of arguments");
    
  /* TODO: Change This */
  udata.neval = 0;
  udata.verbose = 0;
   
  /* First term is the function handle or name */
  func_name = prhs[0];

  if (mxIsChar(func_name))
    {
      CHECK0(mxGetString(func_name, udata.f, FLEN) == 0,
	     "error reading function name string (too long?)");
      udata.nrhs = 1;
      udata.xrhs = 0;
    }
#ifndef HAVE_OCTAVE
  else if (mxIsFunctionHandle(func_name))
    {
      udata.prhs[0] = (mxArray *)func_name;
      strcpy(udata.f, "feval");
      udata.nrhs = 2;
      udata.xrhs = 1;
      }
#endif
  else
    {
      mexErrMsgTxt("First term should be a function name "
		   "(Matlab/Octave) or function handle (Matlab)");
    }

  /* Second parameter. Set of values. */
  CHECK0(mxIsDouble(prhs[1]) && !mxIsComplex(prhs[1]) &&
	 mxGetNumberOfDimensions(prhs[1]) == 2,
	 "The set of values must be a 2D real matrix.");

  nDim = mxGetM(prhs[1]);
  nPoints = mxGetN(prhs[1]);
  xset = mxGetPr(prhs[1]);
  mexPrintf("Loading set of values. nDims=%i, nPoints=%i\n",nDim,nPoints);

  udata.prhs[udata.xrhs] = mxCreateDoubleMatrix(1, nDim, mxREAL);

  xopt = mxCreateDoubleMatrix(1, nDim, mxREAL);
  xptr = mxGetPr(xopt);
    

  /* Third term. Parameters  */
  if (nrhs != 2)
    {
      CHECK0(mxIsStruct(prhs[2]), "3rd element must be a struct");
      params = prhs[2];
    }
  else
    {
      params = mxCreateStructMatrix(1,1,0,NULL);
    }

  parameters = load_parameters(params);
  
  error_code = bayes_optimization_disc(nDim,user_function,&udata,xset,nPoints,
				       xptr,&fmin,parameters);

  mxDestroyArray(udata.prhs[udata.xrhs]);
  plhs[0] = xopt;
  if (nlhs > 1) 
    {
      plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
      *(mxGetPr(plhs[1])) = fmin;
    }
    
  if (nlhs > 2)
    {
      plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
      *(mxGetPr(plhs[2])) = (double)(error_code);
    }
    

}
