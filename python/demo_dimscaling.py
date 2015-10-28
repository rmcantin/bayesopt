#!/usr/bin/env python
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

# This example was provided by Janto Dreijer <jantod@gmail.com>

import math
import numpy as np
import bayesopt

def quad(x,mu):
    return ((np.asarray(x) - mu)**2).mean()

def func(x):
    #print "x", x
    #~ target = np.ones(len(x))*0.3
    target = np.arange(1,1+len(x))
    target2 = np.ones(len(x))*10
    #print "target", target
    e = quad(x,target)
    return e

# Initialize the parameters by default
params = {} #bayesopt.initialize_params()

# We decided to change some of them
params['n_init_samples'] = 30
params['n_iter_relearn'] = 20
params['noise'] = 1e-10
params['kernel_name'] = "kMaternISO5"
params['kernel_hp_mean'] = [1]
params['kernel_hp_std'] = [5]
params['surr_name'] = "sStudentTProcessNIG"
#params['crit_name'] = "cMI"

dim = 20
lb = np.ones((dim,))*0
ub = np.ones((dim,))*20

mvalue, x_out, error = bayesopt.optimize(func, dim, lb, ub, params)

print "Result", mvalue, x_out

print "Global optimal", 0, np.arange(1,1+dim)

print "Y Gap", mvalue
print "X Gap", math.sqrt(mvalue*dim)
