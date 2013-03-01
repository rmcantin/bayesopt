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

static const char* const BOPT_KERNEL_STRINGS[] = {
  "MATERN_ISO1", "MATERN_ISO3", "MATERN_ISO5", "SE_ISO", "SE_ARD"
};

// static const char* const BOPT_CRITERIA_STRINGS[] = {
// };


kernel_name str2kernel(const char* name)
{
  size_t nk = 5;
  for(size_t i= 0;i<nk;++i)
    if (strcmp(name,BOPT_KERNEL_STRINGS[i]) == 0) 
      return static_cast<kernel_name>(i);

  return K_ERROR;
}

const char* kernel2str(kernel_name name)
{
  if(name == K_ERROR)
    return "ERROR!";
  else
    return BOPT_KERNEL_STRINGS[name];
}

criterium_name str2crit(const char* name)
{
  if      (!strcmp(name,  "EI"))
    return C_EI;
  else if (!strcmp(name,  "LCB"))
    return C_LCB;
  else if (!strcmp(name,  "EI_A"))
    return C_EI_A;
  else if (!strcmp(name,  "LCB_A"))
    return C_LCB_A;
  else if (!strcmp(name,  "POI"))
    return C_POI;
  else if (!strcmp(name,  "GREEDY_A_OPTIMALITY"))
    return C_GREEDY_A_OPTIMALITY;
  else if (!strcmp(name,  "EXPECTED_RETURN"))
    return C_EXPECTED_RETURN;
  else if (!strcmp(name,  "THOMPSON_SAMPLING"))
    return C_THOMPSON_SAMPLING;
  else if (!strcmp(name,  "OPTIMISTIC_SAMPLING"))
    return C_OPTIMISTIC_SAMPLING;
  else if (!strcmp(name,  "GP_HEDGE"))
    return C_GP_HEDGE;
  else if (!strcmp(name,  "GP_HEDGE_RANDOM"))
    return C_GP_HEDGE_RANDOM;
  else return C_ERROR;
}


const char* crit2str(criterium_name name)
{
  switch(name)
    {
    case C_EI: return "EI"; 
    case C_LCB: return "LBC"; 
    case C_EI_A: return "EI_A"; 
    case C_LCB_A: return "LBC_A"; 
    case C_POI: return "POI"; 
    case C_GREEDY_A_OPTIMALITY: return "GREEDY_A_OPTIMALITY";
    case C_EXPECTED_RETURN: return "EXPECTED_RETURN";
    case C_THOMPSON_SAMPLING: return "THOMPSON_SAMPLING";
    case C_OPTIMISTIC_SAMPLING: return "OPTIMISTIC_SAMPLING";
    case C_GP_HEDGE: return "GP_HEDGE"; 
    case C_GP_HEDGE_RANDOM: return "GP_HEDGE_RANDOM"; 
    case C_ERROR:
    default: return "ERROR!";
    }
}


surrogate_name str2surrogate(const char* name)
{
  if      (!strcmp(name,  "GAUSSIAN_PROCESS"))
    return S_GAUSSIAN_PROCESS;
  else if (!strcmp(name,  "GAUSSIAN_PROCESS_INV_GAMMA_NORMAL"))
    return S_GAUSSIAN_PROCESS_INV_GAMMA_NORMAL;
  else if (!strcmp(name,  "STUDENT_T_PROCESS_JEFFREYS"))
    return S_STUDENT_T_PROCESS_JEFFREYS;
  else return S_ERROR;
}


const char* surrogate2str(surrogate_name name)
{
  switch(name)
    {
    case S_GAUSSIAN_PROCESS: return "GAUSSIAN_PROCESS"; 
    case S_GAUSSIAN_PROCESS_INV_GAMMA_NORMAL: return "GAUSSIAN_PROCESS_INV_GAMMA_NORMAL"; 
    case S_STUDENT_T_PROCESS_JEFFREYS: return "S_STUDENT_T_PROCESS_JEFFREYS"; 
    case S_ERROR:
    default: return "ERROR!";
    }
}



mean_name str2mean(const char* name)
{
  if      (strcmp(name, "ZERO")   == 0)
    return M_ZERO;
  else if (strcmp(name, "ONE")    == 0)
    return M_CONSTANT;
  else if (strcmp(name, "CONSTANT")    == 0)
    return M_CONSTANT;
  else if (strcmp(name, "LINEAR") == 0)
    return M_LINEAR;
  else if (strcmp(name, "LINEAR_CONSTANT")    == 0)
    return M_CONSTANT;
  else return M_ERROR;
}

const char* mean2str(mean_name name)
{
  switch(name)
    {
    case M_ZERO: return "ZERO"; 
    case M_ONE: return "ONE"; 
    case M_CONSTANT: return "CONSTANT"; 
    case M_LINEAR: return "LINEAR"; 
    case M_LINEAR_CONSTANT: return "LINEAR_CONSTANT"; 
    case M_ERROR:
    default: return "ERROR!";
    }
}

char DEF_LOG_FILE[] = "bayesopt.log";
char DEF_KERNEL[] = "kMaternISO3";

static const bopt_params DEFAULT_PARAMS = {
  DEFAULT_ITERATIONS, MAX_INNER_EVALUATIONS, DEFAULT_SAMPLES, 
  DEFAULT_VERBOSE, DEF_LOG_FILE,
  {KERNEL_THETA}, {KERNEL_SIGMA}, 1, 
  {MEAN_MU}, 1, 
  PRIOR_ALPHA, PRIOR_BETA, PRIOR_DELTA_SQ,
  DEFAULT_NOISE,
  S_GAUSSIAN_PROCESS,
  K_MATERN_ISO3,
  DEF_KERNEL,
  C_EI, M_ONE
};


bopt_params initialize_parameters_to_default(void)
{
  return DEFAULT_PARAMS;
}
