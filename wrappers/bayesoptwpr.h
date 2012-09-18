/** -*- c -*- \file bayesoptwpr.hpp \brief C-wrapper for Bayesian optimization */
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

#ifndef _BAYESOPTWPR_H_
#define _BAYESOPTWPR_H_

#include "ctypes.h"

/** \addtogroup BayesOptimization */
/*@{*/

#ifdef __cplusplus
extern "C" {
#endif 

  /* TODO: make it const double *x */
  typedef double (*eval_func)(unsigned int n, const double *x,
			      double *gradient, /* NULL if not needed */
			      void *func_data);


/** 
 * @brief C functional wrapper for the Bayesian optimization algorithm. 
 * This is an efficient, C/C++ implementation of the Bayesian optimization.
 * Basically, it uses the active learning strategy to optimize an "arbitrary" 
 * funtion using few iterations.
 * 
 */
  int bayes_optimization(int nDim, eval_func f, void* f_data,
			 const double *lb, const double *ub, /* bounds */
			 double *x, /* out: minimizer */
			 double *minf, /* out: minimum */
			 bopt_params parameters);
  
#ifdef __cplusplus
}
#endif 

/**@}*/

#endif
