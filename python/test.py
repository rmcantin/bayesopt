#!/usr/bin/env python
import sys
import numpy as np
import time as tm
import bayesopt as kp
import bayesoptmodule as bopt

# Function for testing
def testfunc(Xin):
    total = 5.0
    for value in Xin:
        total = total + (value -0.33)*(value-0.33)

    return total

class BayesOptTest(bopt.BayesOptModule):
    def evalfunc(self,Xin):
        return testfunc(Xin)

    
# Let's define the parameters
params = {"theta": 0.11, "alpha": 1.0, "beta": 0.1,
          "delta": 1000.0, "noise": 0.001}

# options: see ctypes.cpp
crit = "ei"        
surr = "gp" 
kernel = "materniso3"

n = 5                     # n dimensions
ninit = 40                # n initial iterations
niter = 100               # n iterations

lb = np.zeros((n,))
ub = np.ones((n,))

start = tm.clock()

mvalue, x_out, error = kp.optimize(testfunc, n, lb, ub, ninit,
                                   niter, params, crit, surr,
                                   kernel)

print "Result", x_out
print "Seconds", tm.clock() - start


bo = BayesOptTest()
bo.params = params
bo.crit = crit
bo.surr = surr
bo.kernel = kernel
bo.n = n
bo.niter = niter
bo.lb = lb
bo.ub = ub
bo.ninititer = ninit

start = tm.clock()
mvalue, x_out, error = bo.optimize()

print "Result", x_out
print "Seconds", tm.clock() - start

