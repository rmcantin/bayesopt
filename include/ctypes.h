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
//#include "defaults.h"

#ifdef __cplusplus
extern "C" {
#endif 

  /*** Type definitions                                       **/
  /*************************************************************/
  
  typedef enum {
    K_MATERN_ISO1,
    K_MATERN_ISO3,
    K_MATERN_ISO5,
    K_SE_ISO,
    K_SE_ARD,
    K_ERROR = -1
  } kernel_name;
  
  typedef enum {
    M_ZERO,
    M_ONE,
    M_LINEAR,
    M_ERROR = -1
  } mean_name;

  typedef enum {  
    C_EI,
    C_LCB,
    C_POI,
    C_GP_HEDGE,
    C_GREEDY_A_OPTIMALITY,
    C_EXPECTED_RETURN,
    C_OPTIMISTIC_SAMPLING,
    C_ERROR = -1
  } criterium_name;

  typedef enum {  
    S_GAUSSIAN_PROCESS,
    S_GAUSSIAN_PROCESS_INV_GAMMA_NORMAL,
    S_STUDENT_T_PROCESS_JEFFREYS,
    S_ERROR = -1
  } surrogate_name;

  /** SKO Parameters */
  typedef struct {
    unsigned int n_iterations;   /**< Maximum SKO evaluations (budget) */
    unsigned int n_init_samples; /**< Number of samples before optimization */
    unsigned int verbose_level;  /**< Verbose level */
    double theta;                /**< Kernel hyperparameters */
    double alpha, beta, delta;   /**< Inv-Gamma-Normal hyperparameters */
    double noise;                /**< Observation noise */
    surrogate_name s_name;       /**< Name of the surrogate function */
    kernel_name k_name;          /**< Name of the kernel function */
    criterium_name c_name;       /**< Name of the criterion */
    mean_name m_name;            /**< Name of the mean function */
  } sko_params;


  /* Default values                                            */
  /*************************************************************/

  /* Nonparametric process "parameters" */
  const double KERNEL_THETA    = 0.06;
  const double PRIOR_ALPHA     = 1.0;
  const double PRIOR_BETA      = 1.0;
  const double PRIOR_DELTA_SQ  = 1000.0;
  const double DEFAULT_NOISE   = 1e-4;

  /* Algorithm parameters */
  const size_t DEFAULT_ITERATIONS  = 300;
  const size_t DEFAULT_SAMPLES     = 30;
  const size_t DEFAULT_VERBOSE     = 1;

  /* Algorithm limits */
  const size_t MAX_ITERATIONS  = 1000;       /* Not used */
  const size_t MAX_DIM         = 40;         /* Not used */

  /* INNER Optimizer default values */
  const size_t MAX_INNER_EVALUATIONS = 3000;
  const size_t MAX_INNER_ITERATIONS  = 3000; /* Not used */

  /* Latin Hypercube Sampling (LHS) default values */
  const size_t N_LHS_EVALS_PER_DIM = 30;     /* Not used */
  const size_t MAX_LHS_EVALUATIONS = 100;    /* Not used */

  const size_t nAlgorithmsInGPHedge = 5;
  const criterium_name algorithmsInGPHedge[] = { C_EI, C_LCB, C_POI,
						 C_EXPECTED_RETURN,
						 C_OPTIMISTIC_SAMPLING };


  const sko_params DEFAULT_PARAMS = {
    DEFAULT_ITERATIONS, DEFAULT_SAMPLES, DEFAULT_VERBOSE,
    KERNEL_THETA, 
    PRIOR_ALPHA, PRIOR_BETA, PRIOR_DELTA_SQ,
    DEFAULT_NOISE,
    S_GAUSSIAN_PROCESS,
    K_MATERN_ISO3,
    C_EI
  };

  /* These functions are added to simplify wrapping code */
  kernel_name    str2kernel    (const char* name);
  criterium_name str2crit      (const char* name);
  surrogate_name str2surrogate (const char* name);
  mean_name      str2mean      (const char* name);

  sko_params initialize_parameters_to_default(void);

#ifdef __cplusplus
}
#endif 


#endif
