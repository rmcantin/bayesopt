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

#ifndef __CTYPES_H__
#define __CTYPES_H__

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
    k_matterniso,
    k_seiso,
    k_seard,
  } kernel_name;


  typedef enum {  
    c_ei,
    c_lcb,
    c_poi,
    c_gp_hedge,
    c_greedyAOptimality,
    c_expectedReturn,
    c_optimisticSampling
  } criterium_name;

  typedef enum {  
    s_gaussianProcess,
    s_gaussianProcessHyperPriors,
    s_studentTProcess
  } surrogate_name;

#ifdef __cplusplus
}
#endif 


#endif
