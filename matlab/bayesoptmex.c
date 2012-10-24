/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2012 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
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
  size_t nDim;
  bopt_params parameters = initialize_parameters_to_default();
     
     
  /* Check correct number of parameters */
  CHECK0(nlhs != 2 || nrhs != 3 || nrhs != 5, 
	 "wrong number of arguments");
    
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
      mexErrMsgTxt("First term should be a function name or function handle");
    }

  /* Second parameter. nDim */
  CHECK0(mxIsNumeric(prhs[1]) && !mxIsComplex(prhs[1]) 
	 && mxGetM(prhs[1]) * mxGetN(prhs[1]) == 1,
	 "nDim must be a scalar");
  nDim = (unsigned int) mxGetScalar(prhs[1]);

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

  struct_size(params,"iterations", &parameters.n_iterations);
  struct_size(params,"inner_iterations", &parameters.n_inner_iterations);
  struct_size(params, "init_iterations", &parameters.n_init_samples);
  struct_size(params, "verbose_level", &parameters.verbose_level);

  struct_value(params, "alpha", &parameters.alpha);
  struct_value(params, "beta",  &parameters.beta);
  struct_value(params, "delta", &parameters.delta);
  struct_value(params, "noise", &parameters.noise);

  struct_array(params, "theta", &parameters.n_theta, 
	       &parameters.theta[0]);

  struct_array(params, "mu", &parameters.n_mu, 
	       &parameters.mu[0]);

  /* Extra configuration
  /  See parameters.h for the available options */
  char c_str[100], s_str[100], k_str[100], m_str[100];
  
  strcpy( c_str, crit2str(parameters.c_name));
  strcpy( s_str, surrogate2str(parameters.s_name));
  strcpy( k_str, kernel2str(parameters.k_name));
  strcpy( m_str, mean2str(parameters.m_name));

  struct_string(params, "c_name", c_str);
  parameters.c_name = str2crit(c_str);

  struct_string(params, "s_name", s_str);
  parameters.s_name = str2surrogate(s_str);
  
  struct_string(params, "k_name", k_str);
  parameters.k_name = str2kernel(k_str);

  struct_string(params, "m_name", m_str);
  parameters.m_name = str2mean(m_str);


  double *ub, *lb;    /* Upper and lower bound */
  double fmin;

  if(nrhs == 5)
    {
      /* Load bounds */
      CHECK0(mxIsDouble(prhs[3]) && !mxIsComplex(prhs[3])
	     && (mxGetM(prhs[3]) == 1    || mxGetN(prhs[3]) == 1)
	     && (mxGetM(prhs[3]) == nDim || mxGetN(prhs[3]) == nDim),
	     "lowerBound must be real row or column vector");

      lb = mxGetPr(prhs[3]);
      mexPrintf("Loading bounds \n");

      CHECK0(mxIsDouble(prhs[4]) && !mxIsComplex(prhs[4])
	     && (mxGetM(prhs[4]) == 1    || mxGetN(prhs[4]) == 1)
	     && (mxGetM(prhs[4]) == nDim || mxGetN(prhs[4]) == nDim),
	     "upperBound must be real row or column vector");

      ub = mxGetPr(prhs[4]);
    }
  else
    {
      lb = mxCalloc(nDim,sizeof(double));
      ub = mxCalloc(nDim,sizeof(double));
	 
      unsigned int ii;
      
      for (ii = 0; ii < nDim; ++ii) 
	{
	  lb[ii] = 0.;    
	  ub[ii] = 1.;
	}
    }

  bayes_optimization(nDim,user_function,&udata,lb,ub,xptr,
		     &fmin,parameters);

  if(nrhs != 5)
    {
      mxFree(lb); 
      mxFree(ub);
    }

  mxDestroyArray(udata.prhs[udata.xrhs]);
  plhs[0] = xopt;
  if (nlhs > 1) 
    {
      plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
      *(mxGetPr(plhs[1])) = fmin;
    }
    
}
