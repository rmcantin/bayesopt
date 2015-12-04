## \file bayesopt.pyx \brief Cython wrapper for the BayesOpt Python API
# -------------------------------------------------------------------------
#    This file is part of BayesOpt, an efficient C++ library for
#    Bayesian optimization.
#
#    Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>
#
#    BayesOpt is free software: you can redistribute it and/or modify it
#    under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    BayesOpt is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------


import numpy as np
cimport numpy as np
#from python_ref cimport Py_INCREF, Py_DECREF
from cpython cimport Py_INCREF, Py_DECREF

# Cython > 0.20
#from libc.math cimport HUGE_VAL
# Cython <= 0.19
cdef extern from "math.h" nogil:
    double HUGE_VAL

cdef extern from *:
    ctypedef double* const_double_ptr "const double*"

###########################################################################
cdef extern from "bayesopt/parameters.h":

    ctypedef enum learning_type:
        pass

    ctypedef enum score_type:
        pass

    ctypedef struct kernel_parameters:
        char*  name
        double* hp_mean
        double* hp_std
        unsigned int n_hp

    ctypedef struct mean_parameters:
        char* name
        double* coef_mean
        double* coef_std
        unsigned int n_coef

    ctypedef struct bopt_params:
        unsigned int n_iterations
        unsigned int n_inner_iterations
        unsigned int n_init_samples
        unsigned int n_iter_relearn
        unsigned int init_method
        int random_seed
        int verbose_level
        char* log_filename
        unsigned int load_save_flag
        char* load_filename
        char* save_filename
        char* surr_name
        double sigma_s
        double noise
        double alpha, beta
        score_type sc_type
        learning_type l_type
        double epsilon
        unsigned int force_jump
        kernel_parameters kernel
        mean_parameters mean
        char* crit_name
        double* crit_params
        unsigned int n_crit_params

    learning_type str2learn(char* name)
    char* learn2str(learning_type name)

    score_type str2score(char* name)
    char* score2str(score_type name)

    void set_kernel(bopt_params* params, char* name)
    void set_mean(bopt_params* params, char* name)
    void set_criteria(bopt_params* params, char* name)
    void set_surrogate(bopt_params* params, char* name)
    void set_log_file(bopt_params* params, char* name)
    void set_load_file(bopt_params* params, char* name)
    void set_save_file(bopt_params* params, char* name)
    void set_learning(bopt_params* params, const char* name)
    void set_score(bopt_params* params, const char* name)
    
    bopt_params initialize_parameters_to_default()

###########################################################################
cdef extern from "bayesopt/bayesopt.h":
    ctypedef double (*eval_func)(unsigned int n, const_double_ptr x,
                                 double *gradient, void *func_data)

    int bayes_optimization(int nDim, eval_func f, void* f_data,
                           double *lb, double *ub, double *x,
                           double *minf,
                           bopt_params params)

    int bayes_optimization_disc(int nDim, eval_func f, void* f_data,
                                double *valid_x, size_t n_points,
                                double *x, double *minf,
                                bopt_params parameters)

    int bayes_optimization_categorical(int nDim, eval_func f, void* f_data,
                                        int *categories, double *x,
                                        double *minf, bopt_params parameters)

###########################################################################
cdef bopt_params dict2structparams(dict dparams):

    params = initialize_parameters_to_default()

    params.n_iterations = dparams.get('n_iterations',params.n_iterations)
    params.n_inner_iterations = dparams.get('n_inner_iterations',
                                            params.n_inner_iterations)
    params.n_init_samples = dparams.get('n_init_samples',params.n_init_samples)
    params.n_iter_relearn = dparams.get('n_iter_relearn',params.n_iter_relearn)

    params.init_method = dparams.get('init_method',params.init_method)
    params.random_seed = dparams.get('random_seed',params.random_seed)

    params.verbose_level = dparams.get('verbose_level',params.verbose_level)
    name = dparams.get('log_filename',params.log_filename)
    set_log_file(&params,name)

    params.load_save_flag = dparams.get('load_save_flag',params.load_save_flag)
    l_name = dparams.get('load_filename',params.load_filename)
    set_load_file(&params,l_name)
    s_name = dparams.get('save_filename',params.save_filename)
    set_save_file(&params,s_name)
        
    name = dparams.get('surr_name',params.surr_name)
    set_surrogate(&params,name)

    params.sigma_s = dparams.get('sigma_s',params.sigma_s)
    params.noise = dparams.get('noise',params.noise)
    params.alpha = dparams.get('alpha',params.alpha)
    params.beta = dparams.get('beta',params.beta)

    learning = dparams.get('l_type', None)
    if learning is not None:
        set_learning(&params,learning)

    score = dparams.get('sc_type', None)
    if score is not None:
        set_score(&params,score)
    
    params.epsilon = dparams.get('epsilon',params.epsilon)
    params.force_jump= dparams.get('force_jump',params.force_jump)

    name = dparams.get('kernel_name',params.kernel.name)
    set_kernel(&params,name)

    theta = dparams.get('kernel_hp_mean',None)
    stheta = dparams.get('kernel_hp_std',None)
    if theta is not None and stheta is not None:
        params.kernel.n_hp = len(theta)
        for i in range(0,params.kernel.n_hp):
            params.kernel.hp_mean[i] = theta[i]
            params.kernel.hp_std[i] = stheta[i]

    name = dparams.get('mean_name',params.mean.name)
    set_mean(&params,name)
    
    mu = dparams.get('mean_coef_mean',None)
    smu = dparams.get('mean_coef_std',None)
    if mu is not None and smu is not None:
        params.mean.n_coef = len(mu)
        for i in range(0,params.mean.n_coef):
            params.mean.coef_mean[i] = mu[i]
            params.mean.coef_std[i] = smu[i]

    name = dparams.get('crit_name',params.crit_name)
    set_criteria(&params,name)

    cp = dparams.get('crit_params',None)
    if cp is not None:
        params.n_crit_params = len(cp)
        for i in range(0,params.n_crit_params):
            params.crit_params[i] = cp[i]

    return params

cdef double callback(unsigned int n, const_double_ptr x,
                     double *gradient, void *func_data):
    try:
        x_np = np.zeros(n)

        for i in range(0,n):
            x_np[i] = <double>x[i]

        result = (<object>func_data)(x_np)
        return result
    except:
        return HUGE_VAL

def raise_problem(error_code):
    # This is a little bit hacky, but we lose track of the C++
    # exception since we use the C wrapper for interface:
    # C++ (excep) <-> C (error codes) <-> Python (excep)

    # From bayesoptwpr.cpp
    #static const int BAYESOPT_FAILURE = -1; 
    #static const int BAYESOPT_INVALID_ARGS = -2;
    #static const int BAYESOPT_OUT_OF_MEMORY = -3;
    #static const int BAYESOPT_RUNTIME_ERROR = -4;
    
    if error_code == -1: raise Exception('Unknown error');
    elif error_code == -2: raise ValueError('Invalid argument');
    elif error_code == -3: raise MemoryError;
    elif error_code == -4: raise RuntimeError;

def optimize(f, int nDim, np.ndarray[np.double_t] np_lb,
             np.ndarray[np.double_t] np_ub, dict dparams):

    cdef bopt_params params = dict2structparams(dparams)
    cdef double minf[1]
    cdef np.ndarray np_x = np.ones([nDim], dtype=np.double)*0.5

    cdef np.ndarray[np.double_t, ndim=1, mode="c"] lb
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] ub
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] x

    lb = np.ascontiguousarray(np_lb,dtype=np.double)
    ub = np.ascontiguousarray(np_ub,dtype=np.double)
    x  = np.ascontiguousarray(np_x,dtype=np.double)

    Py_INCREF(f)

    error_code = bayes_optimization(nDim, callback, <void *> f,
                                    &lb[0], &ub[0], &x[0], minf, params)

    
    Py_DECREF(f)

    raise_problem(error_code)
    
    min_value = minf[0]
    return min_value,np_x,error_code


def optimize_discrete(f, np.ndarray[np.double_t,ndim=2] np_valid_x,
                      dict dparams):

    cdef bopt_params params = dict2structparams(dparams)
    
    nDim = np_valid_x.shape[1]

    cdef double minf[1]
    cdef np.ndarray np_x = np.zeros([nDim], dtype=np.double)

    cdef np.ndarray[np.double_t, ndim=1, mode="c"] x
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] valid_x

    x  = np.ascontiguousarray(np_x,dtype=np.double)
    valid_x = np.ascontiguousarray(np_valid_x,dtype=np.double)

    Py_INCREF(f)

    error_code = bayes_optimization_disc(nDim, callback, <void *> f,
                                         &valid_x[0,0], np_valid_x.shape[0],
                                         &x[0], minf, params)

    Py_DECREF(f)

    raise_problem(error_code)
    
    min_value = minf[0]
    return min_value,np_x,error_code

def optimize_categorical(f, np.ndarray[np.int_t,ndim=1] np_categories,
                         dict dparams):

    cdef bopt_params params = dict2structparams(dparams)
    
    nDim = np_categories.shape[0]

    cdef double minf[1]
    cdef np.ndarray np_x = np.zeros([nDim], dtype=np.double)

    cdef np.ndarray[np.double_t, ndim=1, mode="c"] x
    cdef np.ndarray[np.int_t, ndim=1, mode="c"] categories

    x  = np.ascontiguousarray(np_x,dtype=np.double)
    categories = np.ascontiguousarray(np_categories,dtype=np.int)

    Py_INCREF(f)

    error_code = bayes_optimization_categorical(nDim, callback, <void *> f,
                                                <int *>&categories[0], &x[0],
                                                minf, params)

    Py_DECREF(f)

    raise_problem(error_code)
    
    min_value = minf[0]
    return min_value,np_x,error_code
