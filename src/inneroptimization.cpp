#include "krig_config.h"

#ifdef USE_DIRECT_FORTRAN
  #include "direct.hpp"
#else
  // NLOPT
  #include <nlopt.h>
  #include "nloptwpr.hpp"
#endif

#include "inneroptimization.hpp"


int InnerOptimization::innerOptimize(vectord &Xnext)
{   
    double x[128];
    void *objPointer = dynamic_cast<void *>(this);
    int n = static_cast<int>(Xnext.size());
    int error;

    if (objPointer == 0)
      std::cout << "Error casting the current object!" << std::endl;

    error = innerOptimize(x, n, objPointer);

    // There should be a clever way to do this.
    boost::numeric::ublas::array_adaptor<double> shared(n, x);
    boost::numeric::ublas::vector<double, 
				  boost::numeric::ublas::array_adaptor<double> 
				  > Xshared(n, shared); 

    Xnext = Xshared;
    
    return error;
} // nextPoint (uBlas)

int InnerOptimization::innerOptimize(double* x, int n, void* objPointer)
{
    double u[128], l[128];
    double fmin = 1;
    int maxf = MAX_DIRECT_EVALUATIONS;    
    int ierror;

    for (int i = 0; i < n; ++i) {
	l[i] = 0.;
	u[i] = 1.;
    }
 
#ifdef USE_DIRECT_FORTRAN

    if (alg != direct)
      {
	std::cout << "Not supported. Using direct instead." << std::endl;
      }

    int (*fpointer)(int *, double *, double *, 
		    int *, int *,int *, double *,
		    int *, char *, int *, int);
    fpointer = &(DIRECT::evaluate_wrap_);

    int maxT = MAX_DIRECT_ITERATIONS;
    DIRECT::direct(fpointer, x, &n, &fmin, l, u, 
		   &ierror, &maxf, &maxT, objPointer);	


#else /* USE_DIRECT_FORTRAN */
    double (*fpointer)(unsigned int, const double *, double *, void *);
    fpointer = &(NLOPT_WPR::evaluate_nlopt);

    double coef = 0.8;  //Percentaje of resources used in local optimization

    if (alg != combined)  coef = 1.0;

    nlopt_opt opt;

    if ((alg == direct) || (alg == combined))
      {
	opt = nlopt_create(NLOPT_GN_ORIG_DIRECT_L, n); /* algorithm and dims */
      }
    else
      {
	opt = nlopt_create(NLOPT_LD_LBFGS, n); /* algorithm and dims */
      }
    nlopt_set_lower_bounds(opt, l);
    nlopt_set_upper_bounds(opt, u);
    nlopt_set_min_objective(opt, fpointer, objPointer);
    nlopt_set_maxeval(opt, round(maxf*coef) ) ;

    nlopt_result errortype = nlopt_optimize(opt, x, &fmin);

    if ((alg == combined) && (coef < 1)) 
      {
	//opt = nlopt_create(NLOPT_LN_BOBYQA, n); /* algorithm and dims */
	opt = nlopt_create(NLOPT_LN_NELDERMEAD, n); /* algorithm and dims */
	nlopt_set_lower_bounds(opt, l);
	nlopt_set_upper_bounds(opt, u);
	nlopt_set_min_objective(opt, fpointer, objPointer);
	nlopt_set_maxeval(opt, maxf-round(maxf*coef));
	
	errortype = nlopt_optimize(opt, x, &fmin);
      }
      
    if (errortype < 0)
      std::cout << "Error:" << errortype << std::endl;

    ierror = static_cast<int>(errortype);
#endif /* USE_DIRECT_FORTRAN */

    return ierror;

} // nextPoint (C array)


