#Example by Janto Dreijer

from numpy import *

import sys
sys.path.append('/usr/local/lib')
import bayesopt

def func(x):
	#print "x", x
	#~ target = ones(len(x))*0.3
	target = arange(1,1+len(x))
	#print "target", target
	e = ((asarray(x) - target)**2).mean()
	#print "e", e
	return e

# Initialize the parameters by default
params = bayesopt.initialize_params()

# We decided to change some of them
params['n_init_samples'] = 300
params['noise'] = 1
params['k_name'] = "MATERN_ISO1"

dim = 10
lb = ones((dim,))*0
ub = ones((dim,))*20

mvalue, x_out, error = bayesopt.optimize(func, dim, lb, ub, params)

print mvalue, x_out, error
