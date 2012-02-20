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
     if (val) {
	  CHECK0(mxIsNumeric(val) && !mxIsComplex(val) 
		&& mxGetM(val) * mxGetN(val) == 1,
		"param fields must be real scalars");
	  return mxGetScalar(val);
     }
     return dflt;
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
     unsigned n;
     double *x, *x0;
     mxArray *x_mx, *func_name;
     user_function_data d;

     /* TODO: Check correct number of parameters */

     CHECK0(nlhs < 2 && (nrhs == 3 || nrhs == 5), 
	    "wrong number of arguments")

     /* TODO: Change This */
     d.neval = 0;
     d.verbose = 0;
   
     /* First term is the function handle or name */
     func_name = prhs[0];

     if (mxIsChar(func_name))
       {
	 CHECK0(mxGetString(func_name, d.f, FLEN) == 0,
		"error reading function name string (too long?)");
	 d.nrhs = 1;
	 d.xrhs = 0;
       }
     else if (mxIsFunctionHandle(func_name))
       {
	 d.prhs[0] = func_name;
	 strcpy(d.f, "feval");
	 d.nrhs = 2;
	 d.xrhs = 1;
       }
     else
       {
	 mexErrMsgTxt("First term should be a function name or function handle");
       }

     CHECK0(mxIsNumeric(prhs[1]) && !mxIsComplex(prhs[1]) 
	    && mxGetM(prhs[1]) * mxGetN(prhs[1]) == 1,
	    "nDim must be a real scalars");

     n = (int) mxGetScalar(prhs[1]);

     /*     n = mxGetM(prhs[1]) * mxGetN(prhs[1]);
	    x0 = mxGetPr(prhs[1]);*/

     d.prhs[d.xrhs] = mxCreateDoubleMatrix(1, n, mxREAL);


     x_mx = mxCreateDoubleMatrix(1, n, mxREAL);
     x = mxGetPr(x_mx);
     /*     memcpy(x, x0, sizeof(double) * n);*/
     
     /* Configure C interface */
     int nIterations;       /* Number of iterations */
     gp_params par;
     
     CHECK0(mxIsStruct(prhs[2]), "3rd element must be a struct");

     par.theta = struct_val_default(prhs[2], "theta", KERNEL_THETA);
     par.alpha = struct_val_default(prhs[2], "alpha", PRIOR_ALPHA);
     par.beta = struct_val_default(prhs[2], "beta", PRIOR_BETA);
     par.delta = struct_val_default(prhs[2], "delta", PRIOR_DELTA_SQ);
     par.noise = struct_val_default(prhs[2], "noise", DEF_REGULARIZER);
     nIterations = (int) struct_val_default(prhs[2], "iterations", 300);

     /* Common configuration
     /  See ctypes.h for the available options */
     criterium_name c_name = c_ei;
     surrogate_name s_name = s_gaussianProcess;

     double *u, *l;
     double fmin;

     if(nrhs == 5)
       {
	 /* Load limits */
	 CHECK0(mxIsDouble(prhs[3]) && !mxIsComplex(prhs[3])
		&& (mxGetM(prhs[3]) == 1 || mxGetN(prhs[3]) == 1)
		&& (mxGetM(prhs[3]) == n || mxGetN(prhs[3]) == n),
		"lowerBound must be real row or column vector");

	 l = mxGetPr(prhs[3]);

	 CHECK0(mxIsDouble(prhs[4]) && !mxIsComplex(prhs[4])
		&& (mxGetM(prhs[4]) == 1 || mxGetN(prhs[4]) == 1)
		&& (mxGetM(prhs[4]) == n || mxGetN(prhs[4]) == n),
		"upperBound must be real row or column vector");

	 u = mxGetPr(prhs[4]);
       }
     else
       {
	 l = mxCalloc(n,sizeof(double));
	 u = mxCalloc(n,sizeof(double));
	 int i;

	 for (i = 0; i < n; ++i) 
	   {
	     l[i] = 0.;    
	     u[i] = 1.;
	   }
       }

     bayes_optimization(n,user_function,&d,l,u,x,&fmin,
			nIterations,par,c_name,s_name);

     if(nrhs == 5)
       {
	 mxFree(l); mxFree(u);
       }

     mxDestroyArray(d.prhs[d.xrhs]);
     plhs[0] = x_mx;
     if (nlhs > 1) {
	  plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	  *(mxGetPr(plhs[1])) = fmin;
     }
    
}
