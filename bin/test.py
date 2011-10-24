import bayesopt as kp
import numpy as np
import time as tm

# Function for testing
def testfunc(Xin):
    total = 10.0
    for value in Xin:
        total = total + (value -0.53)*(value-0.53)

    return total

# Let's define the parameters
params = {"theta": 0.11, "alpha": 1.0, "beta": 0.1,
          "delta": 10.0, "noise": 0.001}
crit = "ei"     # options: ei, lcb, poi, hedge, aoptimal, expmean
surr = "gp"     # options: gp, gpwpriors, stp

n = 5                      # n dimensions
niter = 200                # n iterations

lb = np.zeros((n,))
ub = np.ones((n,))
x = np.zeros((n,))

out = testfunc(x)

start = tm.clock()

mvalue, x_out, error = kp.optimize(testfunc, n, lb, ub, x,
                                   niter, params, crit, surr)

print "Result", x_out
print "Seconds", tm.clock() - start

