#include <iostream>
#include "ctypes.h"


// TODO: All !strcmp should be negated!!!
kernel_name str2kernel(const char* name)
{
  if      (!strcmp(name,  "MATERN_ISO1"))
      return K_MATERN_ISO1; 
  else if (!strcmp(name,  "MATERN_ISO3"))
      return K_MATERN_ISO3; 
  else if (!strcmp(name,  "MATERN_ISO5"))
      return K_MATERN_ISO5; 
  else if (!strcmp(name,  "SE_ISO"))
      return K_SE_ISO;
  else if (!strcmp(name,  "SE_ARD"))
    return K_SE_ARD;
  else
    return K_ERROR;
}

criterium_name str2crit(const char* name)
{
  if      (!strcmp(name,  "EI"))
    return C_EI;
  else if (!strcmp(name,  "LCB"))
    return C_LCB;
  else if (!strcmp(name,  "POI"))
    return C_POI;
  else if (!strcmp(name,  "GP_HEDGE"))
    return C_GP_HEDGE;
  else if (!strcmp(name,  "GREEDY_A_OPTIMALITY"))
    return C_GREEDY_A_OPTIMALITY;
  else if (!strcmp(name,  "EXPECTED_RETURN"))
    return C_EXPECTED_RETURN;
  else if (!strcmp(name,  "OPTIMISTIC_SAMPLING"))
    return C_OPTIMISTIC_SAMPLING;
  else return C_ERROR;
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

mean_name str2mean(const char* name)
{
  if      (strcmp(name, "ZERO")   == 0)
    return M_ZERO;
  else if (strcmp(name, "ONE")    == 0)
    return M_ONE;
  else if (strcmp(name, "LINEAR") == 0)
    return M_LINEAR;
  else return M_ERROR;
}

sko_params initialize_parameters_to_default(void)
{
  return DEFAULT_PARAMS;
}
