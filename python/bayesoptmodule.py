#!/usr/bin/env python
## \file bayesoptmodule.py \brief BayesOpt wrapper for Python interface (OOP)

# ----------------------------------------------------------------------------
#    This file is part of BayesOptimization, an efficient C++ library for 
#    Bayesian optimization.
#
#    Copyright (C) 2011-2014 Ruben Martinez-Cantin <rmcantin@unizar.es>
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
import bayesopt as bo

## @addtogroup BayesOptPython
# Python Module for BayesOpt.
#
# Check test.py and test2.py with examples about how to run this code.
# Those files also provides examples about how to run the optimization
# using a callback pattern. (see bayesopt.pyx)

## Python Module for BayesOptContinuous
#
# Python module to get run BayesOpt library in a OO pattern.
# The objective module should inherit this one and override evalfunc.
class BayesOptContinuous:
    
    ## Let's define the parameters.
    #
    # For different options: see parameters.h and parameters.cpp .
    # If a parameter is not defined, it will be automatically set
    # to a default value.
    def __init__(self):
        ## Library parameters 
        self.params = {} #bo.initialize_params()
        ## n dimensions
        self.n_dim = 5
        ## Lower bounds
        self.lower_bound = np.zeros((self.n_dim,))
        ## Upper bounds
        self.upper_bound = np.ones((self.n_dim,))

    ## Function for testing.
    # It should be overriden.
    def evalfunc(self, x_in):
        total = 10.0
        for value in x_in:
            total = total + (value -0.53)*(value-0.53)
        return total

    ## Main function. Starts the optimization process.
    def optimize(self):
        min_val, x_out, error = bo.optimize(self.evalfunc, self.n_dim,
                                            self.lower_bound, self.upper_bound,
                                            self.params)
        
        return min_val, x_out, error


## Python Module for BayesOptDiscrete
#
# Python module to get run BayesOpt library in a OO pattern.
# The objective module should inherit this one and override evalfunc.
class BayesOptDiscrete:
    
    ## Let's define the parameters.
    #
    # For different options: see parameters.h and parameters.cpp .
    # If a parameter is not define, it will be automatically set
    # to a default value.
    def __init__(self):
        ## Library parameters 
        self.params = {} #bo.initialize_params()
        n_dim = 5    ## n dimensions
        n_sam = 100  ## n samples
        ## Set of discrete points
        self.x_set = np.random.rand(n_sam, n_dim)
        
    ## Function for testing.
    # It should be overriden.
    def evalfunc(self, x_in):
        total = 10.0
        for value in x_in:
            total = total + (value -0.53)*(value-0.53)
        return total

    ## Main function. Starts the optimization process.
    def optimize(self):
        min_val, x_out, error = bo.optimize_discrete(self.evalfunc,
                                                    self.x_set,
                                                    self.params)
        
        return min_val, x_out, error


if __name__ == "__main__":
    BO = BayesOptContinuous()
    __value__, __x__, __err__ = BO.optimize()
    print "Result", __x__
