#ifndef __DEFAULTS_HPP__
#define __DEFAULTS_HPP__
// Default values

// TODO: This is really a hack

const double KERNEL_THETA    = 0.06;
const double PRIOR_ALPHA     = 1.0;
const double PRIOR_BETA      = 1.0;
const double PRIOR_DELTA_SQ  = 1000.0;
const double DEF_REGULARIZER = 1e-4;

const size_t MAX_ITERATIONS  = 1000;
const size_t MAX_DIM         = 40;

// Default values

// INNER Optimizer default values
const size_t MAX_INNER_EVALUATIONS = 1000;
const size_t MAX_INNER_ITERATIONS  = 1000;

// Latin Hypercube Sampling (LHS) default values

const size_t N_LHS_EVALS_PER_DIM = 30;
const size_t MAX_LHS_EVALUATIONS = 100;
	

#endif
