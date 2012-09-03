#include <iostream>
#include "ctypes.h"


// TODO: All !strcmp should be negated!!!
kernel_name str2kernel(const char* name)
{
  if (!strcmp(name,  "materniso1"))
      return k_materniso1; 
  else if (!strcmp(name,  "materniso3"))
      return k_materniso3; 
  else if (!strcmp(name,  "materniso5"))
      return k_materniso5; 
  else if (!strcmp(name,  "seiso"))
      return k_seiso;
  else if (!strcmp(name,  "seard"))
    return k_seard;
  else
    return k_error;
}

criterium_name str2crit(const char* name)
{
  if (!strcmp(name, "ei"))
    return c_ei;
  else if (!strcmp(name,  "lcb"))
    return c_lcb;
  else if (!strcmp(name,  "poi"))
    return c_poi;
  else if (!strcmp(name,  "gphedge"))
    return c_gp_hedge;
  else if (!strcmp(name,  "aopt"))
    return c_greedyAOptimality;
  else if (!strcmp(name,  "expreturn"))
    return c_expectedReturn;
  else if (!strcmp(name,  "optsampling"))
    return c_optimisticSampling;
  else return c_error;
}

surrogate_name str2surrogate(const char* name)
{
  if (!strcmp(name,  "gp"))
    return s_gaussianProcess;
  else if (!strcmp(name,  "gp_ign"))
    return s_gaussianProcessHyperPriors;
  else if (!strcmp(name,  "stp_jef"))
    return s_studentTProcess;
  else return s_error;
}

mean_name str2mean(const char* name)
{
  if (strcmp(name, "one") == 0)
    return m_one;
  else if (strcmp(name, "zero") == 0)
    return m_zero;
  else if (strcmp(name, "linear") == 0)
    return m_linear;
  else return m_error;
}

sko_params initialize_parameters_to_default(void){
  return DEFAULT_PARAMS;
}
