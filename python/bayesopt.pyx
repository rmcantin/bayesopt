#
#   Pyrex wrapper for the bayesian optimization API
#

import numpy as np
cimport numpy as np

cdef extern from "ctypes.h":
    ctypedef struct gp_params:
        double theta
        double alpha
        double beta
        double delta
        double noise

    ctypedef enum criterium_name:
        c_ei,
        c_lcb,
        c_poi,
        c_gp_hedge,
        c_greedyAOptimality,
        c_expectedReturn

    ctypedef enum surrogate_name:
        s_gaussianProcess,
        s_gaussianProcessHyperPriors,
        s_studentTProcess



cdef extern from "krigwpr.h":
    ctypedef double (*eval_func)(unsigned int n, double *x, double *gradient, void *func_data)
    int krigging_optimization(int nDim, eval_func f, void* f_data, double *lb, double *ub, double *x, double *minf, int maxeval, gp_params params, criterium_name c_name, surrogate_name gp_name)


cdef dict2structparams(dict dparams, gp_params *params):
    params.theta = dparams['theta']
    params.alpha = dparams['alpha']
    params.beta = dparams['beta']
    params.delta = dparams['delta']
    params.noise = dparams['noise']

cdef criterium_name find_criterium(str criteria):
    if criteria == 'ei':
        return c_ei


cdef surrogate_name find_surrogate(str surrogate):
    if surrogate == 'gp':
        return s_gaussianProcess



cdef ndarray2pointer(varray, double* p_array):
    # We assume is a vector. We just copy the first row
    n = varray.shape[0];
    for i in range(0,n):
        p_array[i] = varray[i]

        
cdef pointer2ndarray(int n, double* p_array):
    # We assume is a vector. We just copy the first row
    varray = np.zeros(n)
    for i in range(0,n):
        varray[i] = p_array[i]
    return varray
    

cdef double callback(unsigned int n, double *x, double *gradient, void *func_data):
    invector = pointer2ndarray(n,x)
    result = (<object>func_data)(invector)
    return result



def optimize(object f, int nDim, np.ndarray[np.double_t] np_lb, np.ndarray[np.double_t] np_ub, np.ndarray[np.double_t] np_x, int maxeval, dict dparams, str criteria, str surrogate):
    cdef gp_params zero_params
    zero_params.theta = 0
    zero_params.alpha = 0
    zero_params.beta = 0
    zero_params.delta = 0
    zero_params.noise = 0
    
    cdef gp_params params
    dict2structparams(dparams,&zero_params)
    params = zero_params

    cdef criterium_name crit
    crit = find_criterium(criteria)
    cdef surrogate_name surr
    surr = find_surrogate(surrogate)
    
    
    #cdef double lb[1000]
    #cdef double ub[1000]
    #cdef double c_x[1000]
    cdef double minf[1000]

    #ndarray2pointer(np_lb,lb)
    #ndarray2pointer(np_ub,ub)
    #ndarray2pointer(np_x,c_x)
    
    error_code = krigging_optimization(nDim, callback, <void *> f, <double *>np_lb.data, <double *>np_ub.data, <double *>np_x.data, minf, maxeval, params, crit, surr)
    min_value = minf[0]
    #point_min_value = pointer2ndarray(nDim,c_x)
    return min_value,np_x,error_code
