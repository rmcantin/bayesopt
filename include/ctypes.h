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

#ifndef __C_TYPES_H__
#define __C_TYPES_H__

#include <string.h>
#include "defaults.h"

#ifdef __cplusplus
extern "C" {
#endif 

  typedef enum {
    k_materniso1,
    k_materniso3,
    k_materniso5,
    k_seiso,
    k_seard,
    k_error
  } kernel_name;
  
  typedef enum {  
    c_ei,
    c_lcb,
    c_poi,
    c_gp_hedge,
    c_greedyAOptimality,
    c_expectedReturn,
    c_optimisticSampling,
    c_error
  } criterium_name;

  const unsigned int nAlgorithmsInGPHedge = 5;
  const criterium_name algorithmsInGPHedge[] = {c_ei, c_lcb, c_poi, 
					      c_expectedReturn,
					      c_optimisticSampling };

  typedef enum {  
    s_gaussianProcess,
    s_gaussianProcessHyperPriors,
    s_studentTProcess,
    s_error
  } surrogate_name;

  typedef enum {
    m_zero,
    m_one,
    m_linear,
    m_error
  } mean_name;

  /** Parameters */
  typedef struct {
    /** Maximum SKO evaluations (budget) */
    unsigned int n_iterations;
    /** Number of samples before optimization */
    unsigned int n_init_samples;
    /** Verbose level */
    unsigned int verbose_level;
    /** Surrogate function parameters */
    double theta;  
    double alpha, beta, delta;
    /** Observation noise */
    double noise;
    /** Name of the surrogate function */
    surrogate_name s_name;
    /** Name of the kernel */
    kernel_name k_name;
    /** Name of the criteria */
    criterium_name c_name;
  } sko_params;

  const sko_params DEFAULT_PARAMS = {
    300, 30, 1,
    KERNEL_THETA, 
    PRIOR_ALPHA, PRIOR_BETA, PRIOR_DELTA_SQ,
    DEF_REGULARIZER,
    s_gaussianProcess,
    k_materniso3,
    c_ei
  };

  /* These functions are added to simplify wrapping code */
  kernel_name str2kernel(const char* name);
  criterium_name str2crit(const char* name);
  surrogate_name str2surrogate(const char* name);
  mean_name str2mean(const char* name);
  sko_params initialize_parameters_to_default(void);

#ifdef __cplusplus
}
#endif 


#endif
