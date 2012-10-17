## \file bayesopt.pyx \brief Cython wrapper for the BayesOpt Python API

# ----------------------------------------------------------------------------
#    This file is part of BayesOptimization, an efficient C++ library for 
#    Bayesian optimization.
#
#    Copyright (C) 2011-2012 Ruben Martinez-Cantin <rmcantin@unizar.es>
#
#    BayesOptimization is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    BayesOptimization is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with BayesOptimization. If not, see <http://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------------

import numpy as np
cimport numpy as np
#from python_ref cimport Py_INCREF, Py_DECREF
from cpython cimport Py_INCREF, Py_DECREF

cdef extern from *:
    ctypedef double* const_double_ptr "const double*"

###########################################################################
cdef extern from "parameters.h":

    ctypedef enum kernel_name:
        pass

    ctypedef enum mean_name:
        pass
        
    ctypedef enum criterium_name:
        pass
    
    ctypedef enum surrogate_name:
        pass

    ctypedef struct bopt_params:
        unsigned int n_iterations, n_init_samples, verbose_level
        double* theta
        unsigned int n_theta
        double* mu
        unsigned int n_mu
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

    params.alpha = dparams.get('alpha',params.alpha)
    params.beta = dparams.get('beta',params.beta)
    params.delta = dparams.get('delta',params.delta)
    params.noise = dparams.get('noise',params.noise)
    
    
    theta = dparams.get('theta',None)
    if theta is not None:
        params.n_theta = len(theta)
        for i in range(0,params.n_theta):
            params.theta[i] = theta[i]

    mu = dparams.get('mu',None)
    if mu is not None:
        params.n_mu = len(mu)
        for i in range(0,params.n_mu):
            params.mu[i] = mu[i]

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
        "theta"  : [1.0],
        "n_theta": 1,
        "mu"     : [1.0],
        "n_mu"   : 1,
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
