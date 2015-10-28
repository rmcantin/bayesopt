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

import numpy as np
from bayesoptmodule import BayesOptContinuous
from multiprocessing import Process, Pipe

# Function for testing.
def testfunc(Xin):
    total = 5.0
    for value in Xin:
        total = total + (value -0.33)*(value-0.33)
    return total


def worker(pipe):
    x = None
    while True:
        x = pipe.recv()
        if x == 'STOP':
            break

        result = testfunc(x)
        pipe.send(result)


class BayesOptProcess(Process,BayesOptContinuous):

    def __init__ (self, pipe, n_dim):
        Process.__init__(self)
        BayesOptContinuous.__init__(self, n_dim)
        self.pipe = pipe

    def run(self):
        mvalue, x_out, error = self.optimize()
        self.pipe.send('STOP')
        self.mvalue = mvalue
        self.x_out = x_out
        self.error = error

        return

    def evaluateSample(self, x):
        self.pipe.send(x)
        result = self.pipe.recv()
        return result


if __name__ == '__main__':
    params = {
        'n_iterations' : 50,
        'n_init_samples' : 20,
        's_name' : "sGaussianProcessNormal",
        'c_name' : "cHedge(cEI,cLCB,cExpReturn,cOptimisticSampling)"
    } 

    pipe_par, pipe_child = Pipe()

    bo = BayesOptProcess(pipe_child,n_dim=5)
    bo.parameters = params

    p = Process(target=worker, args=(pipe_par,))

    bo.start()
    p.start()
    
    bo.join()
    p.join()
