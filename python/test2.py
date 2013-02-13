#Example by Janto Dreijer

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
params['n_init_samples'] = 300
params['noise'] = 1
params['k_name'] = "MATERN_ISO3"
params['s_name'] = "GAUSSIAN_PROCESS_INV_GAMMA_NORMAL"

dim = 20
lb = np.ones((dim,))*0
ub = np.ones((dim,))*20

mvalue, x_out, error = bayesopt.optimize(func, dim, lb, ub, params)

print mvalue, x_out, error
