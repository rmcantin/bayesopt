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

#include "bayesoptwpr.h"

#define CHECK0(cond, msg) if (!(cond)) mexErrMsgTxt(msg);

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
     double *x, *x0, opt_f;
     nlopt_result ret;
     mxArray *x_mx, *func_name;
     user_function_data d, *dfc = NULL, *dh = NULL;

     // First term is the function handle
     func_name = mxGetPr(prhs[0]);

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

     //TODO: Change This
     d.neval = 0;
     d.verbose = 0;
     CHECK(mxIsDouble(prhs[1]) && !mxIsComplex(prhs[1])
	   && (mxGetM(prhs[1]) == 1 || mxGetN(prhs[1]) == 1),
	   "x must be real row or column vector");
     n = mxGetM(prhs[1]) * mxGetN(prhs[1]),
     x0 = mxGetPr(prhs[1]);

     x_mx = mxCreateDoubleMatrix(mxGetM(prhs[1]), mxGetN(prhs[1]), mxREAL);
     x = mxGetPr(x_mx);
     memcpy(x, x0, sizeof(double) * n);
     
     //TODO:
     ret = nlopt_optimize(opt, x, &opt_f);
}
