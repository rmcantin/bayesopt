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

#include "parameters.h"
#include "bayesoptwpr.h"

#define CHECK0(cond, msg) if (!(cond)) mexErrMsgTxt(msg);

static void struct_value(const mxArray *s, const char *name, double *result)
{
  mxArray *val = mxGetField(s, 0, name);
  if (val) 
    {
      if(!(mxIsNumeric(val) && !mxIsComplex(val) 
	   && mxGetM(val) * mxGetN(val) == 1))
	{
	  mexErrMsgTxt("param fields must be real scalars");
	}
      else
	{
	  *result = mxGetScalar(val);
	}
    }
  return;
}

static void struct_array(const mxArray *s, const char *name, size_t *n, double *result)
{
  mxArray *val = mxGetField(s, 0, name);
  if (val) 
    {
      if(!(mxIsNumeric(val) && !mxIsComplex(val)))
	{
	  mexErrMsgTxt("Param fields must be vector");
	}
      else
	{	   
	  *n = mxGetM(val) * mxGetN(val);
	  memcpy(result, mxGetPr(val), *n * sizeof(double));
	}
    }
  return;
}


static void struct_size(const mxArray *s, const char *name, size_t *result)
{
  mxArray *val = mxGetField(s, 0, name);
  if (val) 
    {
      if(!(mxIsNumeric(val) && !mxIsComplex(val) 
	   && mxGetM(val) * mxGetN(val) == 1))
	{
	  mexErrMsgTxt("param fields must be real scalars");
	}
      else
	{
	  *result = (size_t) mxGetScalar(val);
	}
    }
  return;
}


static void struct_string(const mxArray *s, const char *name, char* result)
{
  mxArray *val = mxGetField(s, 0, name);
  char *str;
  if (!val)  
    {
      mexPrintf("Field not found. Default not modified.");
    }
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
    }
  return;
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



static double user_function(unsigned n, const double *x,
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
	 && (mxGetM(d->plhs[0]) * mxGetN(d->plhs[0]) == 1),
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
  /  See ctypes.h for the available options */
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
