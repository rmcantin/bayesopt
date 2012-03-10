#!/usr/bin/env python
import numpy as np
import bayesopt as kp

class BayesOptModule:
    def __init__(self):
        # Let's define the parameters
        self.params = {"theta": 0.11,
                       "alpha": 1.0, "beta": 0.1,
                       "delta": 10.0, "noise": 0.001}

        # options: see ctypes.cpp
        self.crit = "ei"        
        self.surr = "gp_ign" 
        self.kernel = "materniso3"

        self.n = 5                      # n dimensions
        self.ninititer = 10             # n initial iterations
        self.niter = 100                # n iterations

        self.lb = np.zeros((self.n,))
        self.ub = np.ones((self.n,))

    # Function for testing
    def evalfunc(self,Xin):
        total = 10.0
        for value in Xin:
            total = total + (value -0.53)*(value-0.53)
        return total

    def optimize(self):
        mvalue, x_out, error = kp.optimize(self.evalfunc, self.n,
                                           self.lb, self.ub, self.ninititer,
                                           self.niter, self.params,
                                           self.crit, self.surr, self.kernel)
        
        return mvalue, x_out, error



if __name__ == "__main__":
    bopt = BayesOptModule()
    mvalue, x_out, error = bopt.optimize()
    print "Result", x_out
