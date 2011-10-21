#ifndef _KRIGWPR_H_
#define _KRIGWPR_H_

/** \addtogroup BayesOptimization */
/*@{*/

#ifdef __cplusplus
extern "C" {
#endif 

typedef struct {
  double theta;
  double alpha;
  double beta;
  double delta;
  double noise;
} gp_params;


  typedef double (*eval_func)(unsigned int n, double *x,
			      double *gradient, /* NULL if not needed */
			      void *func_data);


/** 
 * @brief C functional wrapper for the Bayesian optimization algorithm. 
 * This is an efficient, C/C++ implementation of the Bayesian optimization.
 * Basically, it uses the active learning strategy to optimize an "arbitrary" 
 * funtion using few iterations.
 * 
 */
  int krigging_optimization(int nDim, eval_func f, void* f_data,
			    const double *lb, const double *ub, /* bounds */
			    double *x, /* in: initial guess, out: minimizer */
			    double *minf, /* out: minimum */
			    int maxeval, gp_params params,
			    int useEI,
			    int use_cooling_scheme);

#ifdef __cplusplus
}
#endif 

/**@}*/

#endif
