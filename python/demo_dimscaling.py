#!/usr/bin/env python
# -------------------------------------------------------------------------
#    This file is part of BayesOpt, an efficient C++ library for 
#    Bayesian optimization.
#
#    Copyright (C) 2011-2013 Ruben Martinez-Cantin <rmcantin@unizar.es>
# 
#    BayesOpt is free software: you can redistribute it and/or modify it 
#    under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    BayesOpt is distributed in the hope that it will be useful, but 
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------

# This example was provided by Janto Dreijer <jantod@gmail.com>

import sys
sys.path.append('/usr/local/lib')

import numpy as np
import bayesopt

def func(x):
	#print "x", x
	#~ target = np.ones(len(x))*0.3
	target = np.arange(1,1+len(x))
	#print "target", target
	e = ((np.asarray(x) - target)**2).mean()
	#print "e", e
	return e

# Initialize the parameters by default
params = bayesopt.initialize_params()

# We decided to change some of them
params['n_init_samples'] = 150
params['n_iter_relearn'] = 20
#params['noise'] = 0.01
params['kernel_name'] = "kMaternISO3"
params['kernel_hp_mean'] = [1]
params['kernel_hp_std'] = [5]
params['surr_name'] = "sStudentTProcessNIG"

dim = 20
lb = np.ones((dim,))*0
ub = np.ones((dim,))*20

mvalue, x_out, error = bayesopt.optimize(func, dim, lb, ub, params)

print mvalue, x_out, error
