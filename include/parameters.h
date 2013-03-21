/**  \file parameters.h \brief Parameter definitions. */
/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2013 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
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


#ifndef __BOPT_PARAMETERS_H__
#define __BOPT_PARAMETERS_H__

#include <string.h>
#include "dll_stuff.h"

#ifdef __cplusplus
extern "C" {
#endif 

  /*************************************************************/
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
    M_CONSTANT,
    M_LINEAR,
    M_LINEAR_CONSTANT,
    M_ERROR = -1
  } mean_name;

  typedef enum {  
    C_EI,
    C_EI_A,
    C_LCB,
    C_LCB_A,
    C_POI,
    C_GREEDY_A_OPTIMALITY,
    C_EXPECTED_RETURN,
    C_OPTIMISTIC_SAMPLING,
    C_THOMPSON_SAMPLING,
    C_GP_HEDGE,
    C_GP_HEDGE_RANDOM,
    C_ERROR = -1
  } criterium_name;

  typedef enum {  
    S_GAUSSIAN_PROCESS,
    S_GAUSSIAN_PROCESS_INV_GAMMA_NORMAL,
    S_STUDENT_T_PROCESS_JEFFREYS,
    S_ERROR = -1
  } surrogate_name;

  /** \brief Configuration parameters */
  typedef struct {

    size_t n_iterations;         /**< Maximum BayesOpt evaluations (budget) */
    size_t n_inner_iterations;   /**< Maximum inner optimizer evaluations */
    size_t n_init_samples;       /**< Number of samples before optimization */

    size_t verbose_level;        /**< 1-Error,2-Warning,3-Info. 4-6 log file*/
    char* log_filename;          /**< Log file path (if applicable) */

    double theta[128];           /**< Kernel hyperparameters prior (mean) */
    double s_theta[128];         /**< Kernel hyperparameters prior (std) */
    size_t n_theta;              /**< Number of kernel hyperparameters */

    double mu[128];              /**< Mean function hyperparameters */
    size_t n_mu;                 /**< Number of mean funct. hyperparameters */

    double alpha;                /**< Inverse Gamma prior for signal var */
    double beta;                 /**< Inverse Gamma prior for signal var*/
    double delta;                /**< Normal prior for mean hyperparameters */
    double noise;                /**< Observation noise (and nugget) */

    surrogate_name s_name;       /**< Name of the surrogate function */
    kernel_name k_name;          /**< Name of the kernel funct. -DEPRECATED-*/
    char* k_s_name;              /**< Name of the kernel function */
    criterium_name c_name;       /**< Name of the criterion */
    mean_name m_name;            /**< Name of the mean funct. -DEPRECATED-*/
    char* m_s_name;              /**< Name of the mean function */

  } bopt_params;


  /*************************************************************/
  /*** Default values                                        ***/
  /*************************************************************/

  /* Nonparametric process "parameters" */
  const double KERNEL_THETA    = 1.0;
  const double KERNEL_SIGMA    = 10.0;
  const double MEAN_MU         = 1.0;
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
  const size_t MAX_INNER_EVALUATIONS = 500;   /*Used per dimmension */
  const size_t MAX_INNER_ITERATIONS  = 3000; /* Not used */

  /* Latin Hypercube Sampling (LHS) default values */
  const size_t N_LHS_EVALS_PER_DIM = 30;     /* Not used */
  const size_t MAX_LHS_EVALUATIONS = 100;    /* Not used */

  const size_t N_ALGORITHMS_IN_GP_HEDGE = 5;
  const criterium_name ALGORITHMS_IN_GP_HEDGE[] = { C_EI, C_LCB, C_POI,
						    C_EXPECTED_RETURN,
						    C_OPTIMISTIC_SAMPLING };

  /*************************************************************/
  /* These functions are added to simplify wrapping code       */
  /*************************************************************/
  kernel_name    str2kernel    (const char* name);
  criterium_name str2crit      (const char* name);
  surrogate_name str2surrogate (const char* name);
  mean_name      str2mean      (const char* name);

  BAYESOPT_API const char* kernel2str(kernel_name name);
  BAYESOPT_API const char* crit2str(criterium_name name);
  BAYESOPT_API const char* surrogate2str(surrogate_name name);
  BAYESOPT_API const char* mean2str(mean_name name);

  BAYESOPT_API bopt_params initialize_parameters_to_default(void);

#ifdef __cplusplus
}
#endif 


#endif
