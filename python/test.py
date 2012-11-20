#!/usr/bin/env python
import bayesopt
import bayesoptmodule
import numpy as np

from time import clock


# Function for testing.
def testfunc(Xin):
    total = 5.0
    for value in Xin:
        total = total + (value -0.33)*(value-0.33)

    return total

# Class for OO testing.
class BayesOptTest(bayesoptmodule.BayesOptContinuous):
    def evalfunc(self,Xin):
        return testfunc(Xin)




# Let's define the parameters
# For different options: see parameters.h and cpp
# If a parameter is not define, it will be automatically set
# to a default value.
params = bayesopt.initialize_params()
params['n_iterations'] = 50
params['n_init_samples'] = 20
params['s_name'] = "GAUSSIAN_PROCESS_INV_GAMMA_NORMAL"
params['c_name'] = "GP_HEDGE"

print "Callback implementation"

n = 5                     # n dimensions
lb = np.zeros((n,))
ub = np.ones((n,))

start = clock()

mvalue, x_out, error = bayesopt.optimize(testfunc, n, lb, ub, params)

print "Result", x_out
print "Seconds", clock() - start


print "OO implementation"
bo_test = BayesOptTest()
bo_test.params = params
bo_test.n = n
bo_test.lb = lb
bo_test.ub = ub

start = clock()
mvalue, x_out, error = bo_test.optimize()

print "Result", x_out
print "Seconds", clock() - start


print "Callback discrete implementation"
x_set = np.random.rand(100,n)
start = clock()

mvalue, x_out, error = bayesopt.optimize_discrete(testfunc, x_set, params)

print "Result", x_out
print "Seconds", clock() - start

value = np.array([testfunc(i) for i in x_set])
print "Optimun", x_set[value.argmin()]
