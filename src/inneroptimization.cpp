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
#include <cmath>
#include <nlopt.h>
#include "nloptwpr.h"

#include "inneroptimization.hpp"

void checkNLOPTerror(nlopt_result errortype)
{
  //TODO: Raise exceptions.
  switch(errortype)
      {
      case -1: std::cout << "NLOPT: General failure" << std::endl; break;
      case -2: std::cout << "NLOPT: Invalid arguments. Check bounds." << std::endl; break;
      case -3: std::cout << "NLOPT: Out of memory" << std::endl; break;
      case -4: std::cout << "NLOPT Warning: Potential roundoff error. " 
			 << "In general, this can be ignored." 
			 << std::endl; break;
      case -5: std::cout << "NLOPT: Force stop." << std::endl; break;
      }
}


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
	// What if x is undefined?
	if (x[i] < l[i] || x[i] > u[i])
	  x[i]=(l[i]+u[i])/2.0;
    }

    double (*fpointer)(unsigned int, const double *, double *, void *);
    fpointer = &(NLOPT_WPR::evaluate_nlopt);

    double coef = 0.8;  //Percentaje of resources used in local optimization
    if (alg != combined)  coef = 1.0;

    nlopt_opt opt;

    /* algorithm and dims */
    switch(alg)
      {
      case direct:      /* same as combined */
      case combined: 	opt = nlopt_create(NLOPT_GN_DIRECT_L, n); break;
      case bobyqa: 	opt = nlopt_create(NLOPT_LN_BOBYQA, n); break;
      case lbfgs:       opt = nlopt_create(NLOPT_LD_LBFGS, n); break;
      default: std::cout << "Algorithm not supported" << std::endl; return -1;
      }

    nlopt_set_lower_bounds(opt, l);
    nlopt_set_upper_bounds(opt, u);
    nlopt_set_min_objective(opt, fpointer, objPointer);
    int nfeval = static_cast<int>(static_cast<double>(maxf)*coef);
    nlopt_set_maxeval(opt, nfeval) ;


    nlopt_result errortype = nlopt_optimize(opt, x, &fmin);
    checkNLOPTerror(errortype);

    // Local refinement
    if ((alg == combined) && (coef < 1)) 
      {
	nlopt_destroy(opt);  // Destroy previous one
	opt = nlopt_create(NLOPT_LN_SBPLX, n); /* algorithm and dims */
	nlopt_set_lower_bounds(opt, l);
	nlopt_set_upper_bounds(opt, u);
	nlopt_set_min_objective(opt, fpointer, objPointer);
	nlopt_set_maxeval(opt, maxf-nfeval);
	
	errortype = nlopt_optimize(opt, x, &fmin);
	checkNLOPTerror(errortype);
      }
      
    nlopt_destroy(opt);  // Destroy opt
    
    ierror = static_cast<int>(errortype);

    return ierror;

} // nextPoint (C array)


