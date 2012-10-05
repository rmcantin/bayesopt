
#
#   Pyrex wrapper for the bayesian optimization API
#

import numpy as np
cimport numpy as np
#from python_ref cimport Py_INCREF, Py_DECREF
from cpython cimport Py_INCREF, Py_DECREF

cdef extern from *:
    ctypedef double* const_double_ptr "const double*"

###########################################################################
cdef extern from "parameters.h":

    ctypedef enum kernel_name:
        K_MATERN_ISO1,
        K_MATERN_ISO3,
        K_MATERN_ISO5,
        K_SE_ISO,
        K_SE_ARD,
        K_ERROR = -1


    ctypedef enum mean_name:
        M_ZERO,
        M_CONSTANT,
        M_LINEAR,
        M_ERROR = -1
        
        
    ctypedef enum criterium_name:
        C_EI,
        C_LCB,
        C_POI,
        C_GP_HEDGE,
        C_GREEDY_A_OPTIMALITY,
        C_EXPECTED_RETURN,
        C_OPTIMISTIC_SAMPLING,
        C_ERROR = -1
    

    ctypedef enum surrogate_name:
        S_GAUSSIAN_PROCESS,
        S_GAUSSIAN_PROCESS_INV_GAMMA_NORMAL,
        S_STUDENT_T_PROCESS_JEFFREYS,
        S_ERROR = -1


 
    ctypedef struct bopt_params:
        unsigned int n_iterations, n_init_samples, verbose_level
        double theta
        double alpha, beta, delta
        double noise
        surrogate_name s_name
        kernel_name k_name
        criterium_name c_name
        mean_name m_name

    kernel_name str2kernel(char* name)
    criterium_name str2crit(char* name)
    surrogate_name str2surrogate(char* name)
    mean_name str2mean(char* name)

    char* kernel2str(kernel_name name)
    char* crit2str(criterium_name name)
    char* surrogate2str(surrogate_name name)
    char* mean2str(mean_name name)
    
    bopt_params initialize_parameters_to_default()

###########################################################################
cdef extern from "bayesoptwpr.h":
    ctypedef double (*eval_func)(unsigned int n, const_double_ptr x,
                                 double *gradient, void *func_data)

    int bayes_optimization(int nDim, eval_func f, void* f_data,
                           double *lb, double *ub, double *x,
                           double *minf,
                           bopt_params params)

    
###########################################################################
cdef bopt_params dict2structparams(dict dparams):

    params = initialize_parameters_to_default()

    params.n_iterations = dparams.get('n_iterations',params.n_iterations)
    params.n_init_samples = dparams.get('n_init_samples',params.n_init_samples)
    params.verbose_level = dparams.get('verbose_level',params.verbose_level)

    params.theta = dparams.get('theta',params.theta)
    params.alpha = dparams.get('alpha',params.alpha)
    params.beta = dparams.get('beta',params.beta)
    params.delta = dparams.get('delta',params.delta)
    params.noise = dparams.get('noise',params.noise)
    

    criteria = dparams.get('c_name',None)
    if criteria is not None:
        params.c_name = str2crit(criteria)
        
    surrogate = dparams.get('s_name', None)
    if criteria is not None:
        params.s_name = str2surrogate(surrogate)

    kernel = dparams.get('k_name',None)
    if kernel is not None:
        params.k_name = str2kernel(kernel)

    mean = dparams.get('m_name',None)
    if mean is not None:
        params.m_name = str2mean(mean)
    
    return params

cdef double callback(unsigned int n, const_double_ptr x,
                     double *gradient, void *func_data):
    x_np = np.zeros(n)
    for i in range(0,n):
        x_np[i] = <double>x[i]
        result = (<object>func_data)(x_np)
    return result


def initialize_params():
    params = {
        "theta"  : 1.0,
        "alpha"  : 1.0,
        "beta"   : 1.0,
        "delta"  : 1000.0,
        "noise"  : 0.001,
        "c_name" : "EI",        
        "s_name" : "GAUSSIAN_PROCESS" ,
        "k_name" : "MATERN_ISO3",
        "m_name" : "ZERO",
        "n_iterations"   : 300,
        "n_init_samples" : 30,
        "verbose_level"  : 1,
        }
    return params

def optimize(f, int nDim, np.ndarray[np.double_t] np_lb,
             np.ndarray[np.double_t] np_ub, dict dparams):

    cdef bopt_params params = dict2structparams(dparams)    
    cdef double minf[1000]
    cdef np.ndarray np_x = np.zeros([nDim, 1], dtype=np.double)
    
    Py_INCREF(f)
    
    error_code = bayes_optimization(nDim, callback, <void *> f,
                                    <double *>np_lb.data, <double *>np_ub.data,
                                    <double *>np_x.data, minf, params)
    Py_DECREF(f)
    min_value = minf[0]
    return min_value,np_x,error_code
