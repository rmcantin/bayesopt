#!/usr/bin/env python
## \file bayesoptmodule.py \brief Python BayesOpt module

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
import bayesopt as kp

## @addtogroup BayesOptPython
# Python Module for BayesOpt.
#
# Check test.py and test2.py with examples about how to run this code.
# Those files also provides examples about how to run the optimization
# using a callback pattern. (see bayesopt.pyx)

## Python Module for BayesOpt
#
# Python module to get run BayesOpt library in a OO pattern.
# The objective module should inherit this one and override evalfunc.
class BayesOptModule:
    
    ## Let's define the parameters.
    #
    # For different options: see parameters.h and parameters.cpp .
    # If a parameter is not define, it will be automatically set
    # to a default value.
    def __init__(self):
        ## Library parameters 
        self.params = kp.initialize_params()
        ## n dimensions
        self.n = 5
        ## Lower bounds
        self.lb = np.zeros((self.n,))
        ## Upper bounds
        self.ub = np.ones((self.n,))

    ## Function for testing.
    # It should be overriden.
    def evalfunc(self,Xin):
        total = 10.0
        for value in Xin:
            total = total + (value -0.53)*(value-0.53)
        return total

    ## Main function. Starts the optimization process.
    def optimize(self):
        mvalue, x_out, error = kp.optimize(self.evalfunc, self.n,
                                           self.lb, self.ub,self.params)
        
        return mvalue, x_out, error



if __name__ == "__main__":
    bopt = BayesOptModule()
    mvalue, x_out, error = bopt.optimize()
    print "Result", x_out
