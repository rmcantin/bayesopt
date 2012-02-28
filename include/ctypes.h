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

#ifdef __cplusplus
extern "C" {
#endif 

  typedef struct {
    double theta;
    double alpha;
    double beta;
    double delta;
    double noise;
  } gp_params;


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


  typedef enum {  
    s_gaussianProcess,
    s_gaussianProcessHyperPriors,
    s_studentTProcess,
    s_error
  } surrogate_name;

  kernel_name str2kernel(const char* name);
  criterium_name str2crit(const char* name);
  surrogate_name str2surrogate(const char* name);

#ifdef __cplusplus
}
#endif 


#endif
