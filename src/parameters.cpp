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
#include <iostream>
#include "parameters.h"



/*surrogate_name str2surrogate(const char* name)
{
  if      (!strcmp(name,  "GAUSSIAN_PROCESS"))
    return S_GAUSSIAN_PROCESS;
  if      (!strcmp(name,  "GAUSSIAN_PROCESS_ML"))
    return S_GAUSSIAN_PROCESS_ML;
  else if (!strcmp(name,  "STUDENT_T_PROCESS_NORMAL_INV_GAMMA"))
    return S_STUDENT_T_PROCESS_NORMAL_INV_GAMMA;
  else if (!strcmp(name,  "STUDENT_T_PROCESS_JEFFREYS"))
    return S_STUDENT_T_PROCESS_JEFFREYS;
  else return S_ERROR;
}

const char* surrogate2str(surrogate_name name)
{
  switch(name)
    {
    case S_GAUSSIAN_PROCESS: return "GAUSSIAN_PROCESS"; 
    case S_GAUSSIAN_PROCESS_ML: return "GAUSSIAN_PROCESS_ML"; 
    case S_STUDENT_T_PROCESS_NORMAL_INV_GAMMA: return "STUDENT_T_PROCESS_NORMAL_INV_GAMMA";
    case S_STUDENT_T_PROCESS_JEFFREYS: return "STUDENT_T_PROCESS_JEFFREYS"; 
    case S_ERROR:
    default: return "ERROR!";
    }
    }*/


learning_type str2learn(const char* name)
{
  if      (!strcmp(name,  "L_ML"))
    return L_ML;
  else if (!strcmp(name,  "L_MAP"))
    return L_MAP;
  else if (!strcmp(name,  "L_LOO"))
    return L_LOO;
  else return L_ERROR;
}




const char* learn2str(learning_type name)
{
  switch(name)
    {
    case L_ML: return "L_ML"; 
    case L_MAP: return "L_MAP"; 
    case L_LOO: return "L_LOO"; 
    case L_ERROR:
    default: return "ERROR!";
    }
}


char DEF_LOG_FILE[128] = "bayesopt.log";
char DEF_SUR_NAME[128] = "sGaussianProcess";
char DEF_KERNEL_NAME[128] = "kMaternISO3";
char DEF_MEAN_NAME[128] = "mOne";
char DEF_CRITERIA_NAME[128] = "cEI";

static const kernel_parameters DEFAULT_KERNEL = {
  DEF_KERNEL_NAME, {KERNEL_THETA}, {KERNEL_SIGMA}, 1 
};

static const mean_parameters DEFAULT_MEAN = {
  DEF_MEAN_NAME, {MEAN_MU}, {MEAN_SIGMA}, 1
};


static const bopt_params DEFAULT_PARAMS = {
  DEFAULT_ITERATIONS, MAX_INNER_EVALUATIONS, 
  DEFAULT_SAMPLES, 0, 
  DEFAULT_VERBOSE, DEF_LOG_FILE,
  DEF_SUR_NAME,
  DEFAULT_SIGMA, DEFAULT_NOISE,
  PRIOR_ALPHA, PRIOR_BETA, 
  L_MAP, 0.0,
  DEFAULT_KERNEL, DEFAULT_MEAN,
  DEF_CRITERIA_NAME, {1.0}, 1
};


bopt_params initialize_parameters_to_default(void)
{
  return DEFAULT_PARAMS;
}
