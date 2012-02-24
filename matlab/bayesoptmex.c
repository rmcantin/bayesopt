/*
-----------------------------------------------------------------------------
   This file is part of BayesOptimization, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOptimization is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOptimization is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with BayesOptimization.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>

#include "defaults.h"
#include "bayesoptwpr.h"

#define CHECK0(cond, msg) if (!(cond)) mexErrMsgTxt(msg);

static double struct_val_default(const mxArray *s, const char *name, double dflt)
{
  mxArray *val = mxGetField(s, 0, name);
  if (val) 
    {
      if(!(mxIsNumeric(val) && !mxIsComplex(val) 
	   && mxGetM(val) * mxGetN(val) == 1))
	{
	  mexErrMsgTxt("param fields must be real scalars");
	  return dflt;
	}
      return mxGetScalar(val);
    }
  return dflt;
}

static void struct_str_default(const mxArray *s, 
			       const char *name, 
			       char* dflt,
			       char* result)
{
  mxArray *val = mxGetField(s, 0, name);
  char *str;
  result = dflt;
  if (!val)  
    return;
  else
    {
      mwSize strlen = mxGetM(val) * mxGetN(val) + 1;
      if( !mxIsChar(val) )
	{
	  mexErrMsgTxt("Method name must be a string");
	  return;
	}
      if (mxGetString(val, result, strlen))
	{
	  mexErrMsgTxt("Error reading string");
	  return;
	}
      result = str;
    }
}


#define FLEN 128 /* max length of user function name */
#define MAXRHS 2 /* max nrhs for user function */
typedef struct {
     char f[FLEN];
     mxArray *plhs[2];
     mxArray *prhs[MAXRHS];
     int xrhs, nrhs;
     int verbose, neval;
} user_function_data;



static double user_function(unsigned n, double *x,
			    double *gradient, /* NULL if not needed */
			    void *d_)
{
  user_function_data *d = (user_function_data *) d_;
  double f;

  
  d->plhs[0] = d->plhs[1] = NULL;

  memcpy(mxGetPr(d->prhs[d->xrhs]), x, n * sizeof(double));

  CHECK0(0 == mexCallMATLAB(gradient ? 2 : 1, d->plhs, 
			   d->nrhs, d->prhs, d->f),
	"error calling user function");

  CHECK0(mxIsNumeric(d->plhs[0]) && !mxIsComplex(d->plhs[0]) 
	&& mxGetM(d->plhs[0]) * mxGetN(d->plhs[0]) == 1,
	"user function must return real scalar");
  f = mxGetScalar(d->plhs[0]);
  mxDestroyArray(d->plhs[0]);

  if (gradient) {
     CHECK0(mxIsDouble(d->plhs[1]) && !mxIsComplex(d->plhs[1])
	   && (mxGetM(d->plhs[1]) == 1 || mxGetN(d->plhs[1]) == 1)
	   && mxGetM(d->plhs[1]) * mxGetN(d->plhs[1]) == n,
	   "gradient vector from user function is the wrong size");
     memcpy(gradient, mxGetPr(d->plhs[1]), n * sizeof(double));
     mxDestroyArray(d->plhs[1]);
  }
  d->neval++;

  if (d->verbose) mexPrintf("Optimize eval #%d: %g\n", d->neval, f);
  return f;
}



void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  double *xptr;
  mxArray *xopt;
  const mxArray *func_name, *params;
  user_function_data udata;
  unsigned int nDim, nIterations;    
  gp_params par;
     
     
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
  else if (mxIsFunctionHandle(func_name))
    {
      udata.prhs[0] = func_name;
      strcpy(udata.f, "feval");
      udata.nrhs = 2;
      udata.xrhs = 1;
    }
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

  par.theta = struct_val_default(params, "theta", KERNEL_THETA);
  par.alpha = struct_val_default(params, "alpha", PRIOR_ALPHA);
  par.beta = struct_val_default(params, "beta", PRIOR_BETA);
  par.delta = struct_val_default(params, "delta", PRIOR_DELTA_SQ);
  par.noise = struct_val_default(params, "noise", DEF_REGULARIZER);
  nIterations = (unsigned int) struct_val_default(params, "iterations", 300);

  /* Extra configuration
  /  See ctypes.h for the available options */
  criterium_name c_name;
  surrogate_name s_name;
  kernel_name k_name;
  char *c_str, *s_str, *k_str;

  struct_str_default(params, "criteria", "ei", c_str);
  c_name = str2crit(c_str);

  struct_str_default(params, "surrogate", "gp", s_str);
  s_name = str2surrogate(s_str);
  
  struct_str_default(params, "kernel", "materniso", k_str);
  k_name = str2kernel(s_str);
  

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
		     &fmin,nIterations,par,c_name,s_name,k_name);

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
