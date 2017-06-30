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

import sys
import bayesopt
from bayesoptmodule import BayesOptContinuous
import numpy as np

from time import clock

# Python3 compat
if sys.version_info[0] == 3:
    raw_input = input

# Function for testing.
def testfunc(Xin):
    total = 5.0
    for value in Xin:
        total = total + (value -0.33)*(value-0.33)

    return total

# Class for OO testing.
class BayesOptTest(BayesOptContinuous):
    def evaluateSample(self,Xin):
        return testfunc(Xin)


# Let's define the parameters
# For different options: see parameters.h and cpp
# If a parameter is not define, it will be automatically set
# to a default value.
params = {}
params['n_iterations'] = 50
params['n_iter_relearn'] = 5
params['n_init_samples'] = 2

print("Callback implementation")

n = 5                     # n dimensions
lb = np.zeros((n,))
ub = np.ones((n,))

start = clock()
mvalue, x_out, error = bayesopt.optimize(testfunc, n, lb, ub, params)

print("Result", mvalue, "at", x_out)
print("Running time:", clock() - start, "seconds")
raw_input('Press INTRO to continue')

print("OO implementation")
bo_test = BayesOptTest(n)
bo_test.parameters = params
bo_test.lower_bound = lb
bo_test.upper_bound = ub

start = clock()
mvalue, x_out, error = bo_test.optimize()

print("Result", mvalue, "at", x_out)
print("Running time:", clock() - start, "seconds")
raw_input('Press INTRO to continue')

print("Callback discrete implementation")
x_set = np.random.rand(100,n)
start = clock()

mvalue, x_out, error = bayesopt.optimize_discrete(testfunc, x_set, params)

print("Result", mvalue, "at", x_out)
print("Running time:", clock() - start, "seconds")

value = np.array([testfunc(i) for i in x_set])
print("Optimum", value.min(), "at", x_set[value.argmin()])
