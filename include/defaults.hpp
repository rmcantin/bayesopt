#ifndef __DEFAULTS_HPP__
#define __DEFAULTS_HPP__
// Default values

// TODO: This is really a bad hack

//#define FIBONACCI_SEED  123u
//#define MT_SEED         156u
#define KERNEL_THETA    0.06
#define KERNEL_P        1.6
#define PRIOR_ALPHA     1.0
#define PRIOR_BETA      0.1
#define PRIOR_DELTA_SQ  10.0
#define DEF_REGULARIZER 1e-4
#define MAX_ITERATIONS  300
#define MAX_DIM         20

// Default values

// DIRECT default values
#define MAX_DIRECT_EVALUATIONS  1000
#define MAX_DIRECT_ITERATIONS   300

// Latin Hypercube Sampling (LHS) default values
#define N_LHS_EVALS_PER_DIM     30
#define MAX_LHS_EVALUATIONS     100
	

#endif
