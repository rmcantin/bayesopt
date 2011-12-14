/*
-----------------------------------------------------------------------------
   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/
#include "krig_config.h"

#ifdef USE_DIRECT_FORTRAN
  #include "directwpr.h"
#else  // NLOPT
  #include <nlopt.h>
  #include "nloptwpr.h"
#endif

#include "inneroptimization.hpp"


int InnerOptimization::innerOptimize(vectord &Xnext)
{   
    double x[128];
    void *objPointer = static_cast<void *>(this);
    int n = static_cast<int>(Xnext.size());
    int error;

    if (objPointer == 0)
      std::cout << "Error casting the current object!" << std::endl;

    for (int i = 0; i < n; ++i) 
	x[i] = Xnext(i);
 
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
    int maxf = MAX_INNER_EVALUATIONS;    
    int ierror;

    for (int i = 0; i < n; ++i) {
	l[i] = mDown;	u[i] = mUp;
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

    int maxT = MAX_INNER_ITERATIONS;
    DIRECT::direct(fpointer, x, &n, &fmin, l, u, 
		   &ierror, &maxf, &maxT, objPointer);	


#else /* USE_DIRECT_FORTRAN */
    double (*fpointer)(unsigned int, const double *, double *, void *);
    fpointer = &(NLOPT_WPR::evaluate_nlopt);

    double coef = 0.8;  //Percentaje of resources used in local optimization
    if (alg != combined)  coef = 1.0;

    nlopt_opt opt;

    /* algorithm and dims */
    switch(alg)
      {
      case direct:      /* same as combined */
      case combined: 	opt = nlopt_create(NLOPT_GN_ORIG_DIRECT_L, n); break;
      case bobyqa: 	opt = nlopt_create(NLOPT_LN_BOBYQA, n); break;
      case lbfgs:       opt = nlopt_create(NLOPT_LD_LBFGS, n); break;
      default: std::cout << "Algorithm not supported" << std::endl; return -1;
      }

    nlopt_set_lower_bounds(opt, l);
    nlopt_set_upper_bounds(opt, u);
    nlopt_set_min_objective(opt, fpointer, objPointer);
    nlopt_set_maxeval(opt, round(maxf*coef) ) ;


    nlopt_result errortype = nlopt_optimize(opt, x, &fmin);

    // Local refinement
    if ((alg == combined) && (coef < 1)) 
      {
	nlopt_destroy(opt);  // Destroy previous one
	opt = nlopt_create(NLOPT_LN_SBPLX, n); /* algorithm and dims */
	nlopt_set_lower_bounds(opt, l);
	nlopt_set_upper_bounds(opt, u);
	nlopt_set_min_objective(opt, fpointer, objPointer);
	nlopt_set_maxeval(opt, maxf-round(maxf*coef));
	
	errortype = nlopt_optimize(opt, x, &fmin);
      }
      
    nlopt_destroy(opt);  // Destroy opt
    if (errortype < 0)
      std::cout << "Error:" << errortype << std::endl;

    ierror = static_cast<int>(errortype);
#endif /* USE_DIRECT_FORTRAN */

    return ierror;

} // nextPoint (C array)


