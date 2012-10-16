#!/usr/bin/env python
import sys
import numpy as np
import time as tm
import bayesopt as kp
import bayesoptmodule as bopt

# Function for testing.
def testfunc(Xin):
    total = 5.0
    for value in Xin:
        total = total + (value -0.33)*(value-0.33)

    return total

# Class for OO testing.
class BayesOptTest(bopt.BayesOptModule):
    def evalfunc(self,Xin):
        return testfunc(Xin)




# Let's define the parameters
# For different options: see parameters.h and cpp
# If a parameter is not define, it will be automatically set
# to a default value.
params = kp.initialize_params()
params['n_iterations'] = 100
params['s_name'] = "GAUSSIAN_PROCESS_INV_GAMMA_NORMAL"
params['c_name'] = "GP_HEDGE"

n = 5                     # n dimensions
lb = np.zeros((n,))
ub = np.ones((n,))

start = tm.clock()

mvalue, x_out, error = kp.optimize(testfunc, n, lb, ub, params)

print "Callback implementation"
print "Result", x_out
print "Seconds", tm.clock() - start


bo = BayesOptTest()
bo.params = params
bo.n = n
bo.lb = lb
bo.ub = ub

start = tm.clock()
mvalue, x_out, error = bo.optimize()

print "OO implementation"
print "Result", x_out
print "Seconds", tm.clock() - start

