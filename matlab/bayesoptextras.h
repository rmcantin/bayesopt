/**  \file bayesoptextras.h 
              \brief Helper functions to Matlab/Octave wrappers. */
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

#ifndef __BAYESOPTEXTRAS_H__
#define __BAYESOPTEXTRAS_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>

#include "bayesopt/bayesopt.h"

/* Declaration */

static void struct_value(const mxArray *s, const char *name, double *result);
static void struct_array(const mxArray *s, const char *name, size_t *n, double *result);
static void struct_size(const mxArray *s, const char *name, size_t *result);
static void struct_int(const mxArray *s, const char *name, int *result);
static void struct_string(const mxArray *s, const char *name, char* result);

static double user_function(unsigned n, const double *x,
		     double *gradient, /* NULL if not needed */
		     void *d_);

/* See parameters.h for the available options */
static bopt_params load_parameters(const mxArray* params);

#define FLEN 128 /* max length of user function name */
#define MAXRHS 2 /* max nrhs for user function */
typedef struct {
     char f[FLEN];
     mxArray *plhs[2];
     mxArray *prhs[MAXRHS];
     int xrhs, nrhs;
     int verbose, neval;
} user_function_data;

/* Implementation */

#define CHECK0(cond, msg) if (!(cond)) mexErrMsgTxt(msg);


void struct_value(const mxArray *s, const char *name, double *result)
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
  else
    {
      mexPrintf("Field %s not found. Default not modified.\n", name);
    }
  return;
}

void struct_array(const mxArray *s, const char *name, size_t *n, double *result)
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
  else
    {
      mexPrintf("Field %s not found. Default not modified.\n", name);
    }
  return;
}


void struct_size(const mxArray *s, const char *name, size_t *result)
{
  mxArray *val = mxGetField(s, 0, name);
  if (val) 
    {
      if(!(mxIsNumeric(val) && !mxIsComplex(val) 
	   && mxGetM(val) * mxGetN(val) == 1))
	{
	  mexErrMsgTxt("param fields must be scalar");
	}
      else
	{
	  *result = (size_t)(mxGetScalar(val));
	}
    }
  else
    {
      mexPrintf("Field %s not found. Default not modified.\n", name);
    }
  return;
}

void struct_int(const mxArray *s, const char *name, int *result)
{
  mxArray *val = mxGetField(s, 0, name);
  if (val) 
    {
      if(!(mxIsNumeric(val) && !mxIsComplex(val) 
	   && mxGetM(val) * mxGetN(val) == 1))
	{
	  mexErrMsgTxt("param fields must be scalar");
	}
      else
	{
	  *result = (int)(mxGetScalar(val));
	}
    }
  else
    {
      mexPrintf("Field %s not found. Default not modified.\n", name);
    }
  return;
}



void struct_string(const mxArray *s, const char *name, char* result)
{
  mxArray *val = mxGetField(s, 0, name);

  if (val) 
    {
      if( mxIsChar(val) ) 
	{
	  /* Using MSVC 2010, data is lost unless we copy it. It seems
	     that mxGetField does not get the pointer and create a new
	     array or mxArrayToString fails to copy the string.   */
	  strcpy(result,mxArrayToString(val));
	} 
      else 
	{
	  mexErrMsgTxt("Method name must be a string");
	}
    } 
  else 
    {
      mexPrintf("Field %s not found. Default not modified.\n", name);
    }
  return;
}



double user_function(unsigned n, const double *x,
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

bopt_params load_parameters(const mxArray* params)
{  
  
  /* See parameters.h for the available options */
  
  char l_str[100], sc_str[100], name[100];
  size_t n_hp_test, n_coef_test;

  bopt_params parameters = initialize_parameters_to_default();

  n_hp_test = parameters.kernel.n_hp;
  n_coef_test = parameters.mean.n_coef;

  struct_size(params,"n_iterations", &parameters.n_iterations);
  struct_size(params,"n_inner_iterations", &parameters.n_inner_iterations);
  struct_size(params, "n_init_samples", &parameters.n_init_samples);
  struct_size(params, "n_iter_relearn", &parameters.n_iter_relearn);

  struct_size(params, "init_method", &parameters.init_method);
  struct_int(params, "random_seed", &parameters.random_seed);
  
  struct_int(params, "verbose_level", &parameters.verbose_level);
  struct_string(params, "log_filename", parameters.log_filename);

  struct_size(params, "load_save_flag", &parameters.load_save_flag);
  struct_string(params, "load_filename", parameters.load_filename);
  struct_string(params, "save_filename", parameters.save_filename);

  
  struct_string(params, "surr_name", parameters.surr_name);

  struct_value(params, "sigma_s", &parameters.sigma_s);
  struct_value(params, "noise", &parameters.noise);
  struct_value(params, "alpha", &parameters.alpha);
  struct_value(params, "beta",  &parameters.beta);
  

  strcpy( l_str, learn2str(parameters.l_type));
  struct_string(params, "l_type", l_str);
  parameters.l_type = str2learn(l_str);

  strcpy( sc_str, score2str(parameters.sc_type));
  struct_string(params, "sc_type", sc_str);
  parameters.sc_type = str2score(sc_str);


  struct_value(params, "epsilon",  &parameters.epsilon);
  struct_size(params, "force_jump",  &parameters.force_jump);

  struct_string(params, "crit_name", parameters.crit_name);
  struct_array(params, "crit_params", &parameters.n_crit_params, 
	       &parameters.crit_params[0]);

  /* Kernel parameters */
  struct_string(params, "kernel_name", parameters.kernel.name);
  struct_array(params, "kernel_hp_mean", &parameters.kernel.n_hp, 
	       &parameters.kernel.hp_mean[0]);
  struct_array(params, "kernel_hp_std", &n_hp_test, 
	       &parameters.kernel.hp_std[0]);

  CHECK0(parameters.kernel.n_hp == n_hp_test, 
	 "Error processing kernel parameters");

  /* Mean function parameters */
  struct_string(params, "mean_name", parameters.mean.name);
  struct_array(params, "mean_coef_mean", &parameters.mean.n_coef, 
	       &parameters.mean.coef_mean[0]);

  struct_array(params, "mean_coef_std", &n_coef_test, 
	       &parameters.mean.coef_std[0]);

  CHECK0(parameters.mean.n_coef == n_coef_test, 
	 "Error processing mean parameters");


  return parameters;
}


#endif
