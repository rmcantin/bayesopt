#ifndef __DEFAULTS_H__
#define __DEFAULTS_H__
/* Default values */

/* Nonparametric process "parameters" */
const double KERNEL_THETA    = 0.06;
const double PRIOR_ALPHA     = 1.0;
const double PRIOR_BETA      = 1.0;
const double PRIOR_DELTA_SQ  = 1000.0;
const double DEF_REGULARIZER = 1e-4;

/* Algorithm limits */
const size_t MAX_ITERATIONS  = 1000;
const size_t MAX_DIM         = 40; /* Not used */

/* INNER Optimizer default values */
const size_t MAX_INNER_EVALUATIONS = 3000;
const size_t MAX_INNER_ITERATIONS  = 3000; /* Not used */

/* Latin Hypercube Sampling (LHS) default values */
const size_t N_LHS_EVALS_PER_DIM = 30;  /* Not used */
const size_t MAX_LHS_EVALUATIONS = 100; /* Not used */
	

#endif
