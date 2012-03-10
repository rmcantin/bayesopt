#
#   Pyrex wrapper for the bayesian optimization API
#

import numpy as np
cimport numpy as np
#from python_ref cimport Py_INCREF, Py_DECREF
from cpython cimport Py_INCREF, Py_DECREF

###########################################################################
cdef extern from "ctypes.h":
    ctypedef struct gp_params:
        double theta
        double alpha
        double beta
        double delta
        double noise

    ctypedef enum kernel_name:
        k_materniso1,
        k_materniso3,
        k_materniso5,
        k_seiso,
        k_seard,
        k_error

    ctypedef enum criterium_name:
        c_ei,
        c_lcb,
        c_poi,
        c_gp_hedge,
        c_greedyAOptimality,
        c_expectedReturn,
        c_error

    ctypedef enum surrogate_name:
        s_gaussianProcess,
        s_gaussianProcessHyperPriors,
        s_studentTProcess,
        s_error

    kernel_name str2kernel(char* name)
    criterium_name str2crit(char* name)
    surrogate_name str2surrogate(char* name)


###########################################################################
cdef extern from "bayesoptwpr.h":
    ctypedef double (*eval_func)(unsigned int n, double *x,
                                 double *gradient, void *func_data)

    int bayes_optimization(int nDim, eval_func f, void* f_data,
                           double *lb, double *ub, double *x,
                           double *minf, int maxiniteval, int maxeval,
                           gp_params params,
                           criterium_name c_name,
                           surrogate_name gp_name,
                           kernel_name k_name)


###########################################################################
cdef dict2structparams(dict dparams, gp_params *params):
    params.theta = dparams['theta']
    params.alpha = dparams['alpha']
    params.beta = dparams['beta']
    params.delta = dparams['delta']
    params.noise = dparams['noise']




cdef double callback(unsigned int n, double *x,
                     double *gradient, void *func_data):
    x_np = np.zeros(n)
    for i in range(0,n):
        x_np[i] = x[i]
    result = (<object>func_data)(x_np)
    return result



def optimize(f, int nDim, np.ndarray[np.double_t] np_lb,
             np.ndarray[np.double_t] np_ub,
             int maxiniteval, 
             int maxeval, dict dparams,
             bytes criteria, bytes surrogate, bytes kernel):

    cdef gp_params params
    dict2structparams(dparams,&params)

    cdef criterium_name crit
    crit = str2crit(criteria)

    cdef surrogate_name surr
    surr = str2surrogate(surrogate)

    cdef kernel_name ker
    cdef char *k_name = kernel
    ker = str2kernel(k_name)
    
    cdef double minf[1000]
    cdef np.ndarray np_x = np.zeros([nDim, 1], dtype=np.double)
    
    Py_INCREF(f)
    error_code = bayes_optimization(nDim, callback, <void *> f,
                                    <double *>np_lb.data, <double *>np_ub.data,
                                    <double *>np_x.data, minf, maxiniteval,
                                    maxeval, params, crit, surr, ker)
    Py_DECREF(f)
    min_value = minf[0]
    return min_value,np_x,error_code
