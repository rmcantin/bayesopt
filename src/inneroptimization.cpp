/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2013 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/
#include <cmath>
#include <nlopt.h>
#include "parameters.h"
#include "log.hpp"
#include "inneroptimization.hpp"

namespace bayesopt
{
  void checkNLOPTerror(nlopt_result errortype)
  {
    switch(errortype)
      {
      case -1: FILE_LOG(logERROR) << "NLOPT: General failure"; break;
      case -2: FILE_LOG(logERROR) << "NLOPT: Invalid arguments. Check bounds."; break;
      case -3: FILE_LOG(logERROR) << "NLOPT: Out of memory"; break;
      case -4: FILE_LOG(logERROR) << "NLOPT Warning: Potential roundoff error. " 
				  << "In general, this can be ignored."; break;
      case -5: FILE_LOG(logERROR) << "NLOPT: Force stop."; break;
      default: ;
      }
  }

  NLOPT_Optimization::NLOPT_Optimization(RBOptimizable* rbo, size_t dim)
  { 
    rbobj = new RBOptimizableWrapper(rbo);       rgbobj = NULL;
    alg = DIRECT;                  maxEvals = MAX_INNER_EVALUATIONS;
    mDown = zvectord(dim);         mUp = svectord(dim,1.0);  
  };

  NLOPT_Optimization::NLOPT_Optimization(RGBOptimizable* rgbo, size_t dim)
  { 
    rbobj = NULL;             rgbobj = new RGBOptimizableWrapper(rgbo);
    alg = DIRECT;             maxEvals = MAX_INNER_EVALUATIONS;
    mDown = zvectord(dim);    mUp = svectord(dim,1.0);  
  };

  NLOPT_Optimization::~NLOPT_Optimization()
  {
    if(rbobj != NULL) delete rbobj;
    if(rgbobj != NULL) delete rgbobj;
  }

  int NLOPT_Optimization::run(vectord &Xnext)
  {   
    assert(mDown.size() == Xnext.size());
    assert(mUp.size() == Xnext.size());

    nlopt_opt opt;
    double (*fpointer)(unsigned int, const double *, double *, void *);
    void *objPointer;

    int n = static_cast<int>(Xnext.size());
    double fmin = 1;
    int maxf1 = maxEvals*n;
    int maxf2 = 0;    // For a second pass
    double coef_local = 0.2;
    int ierror;

    switch(alg)
      {
      case DIRECT: // Pure global. No gradient
	opt = nlopt_create(NLOPT_GN_DIRECT_L, n); 
	fpointer = &(NLOPT_Optimization::evaluate_nlopt);
	objPointer = static_cast<void *>(rbobj);
	break;
      case COMBINED: // Combined local-global (80% DIRECT -> 20% BOBYQA). No gradient
	opt = nlopt_create(NLOPT_GN_DIRECT_L, n); 
	maxf2 = static_cast<int>(static_cast<double>(maxf1)*coef_local);
	maxf1 -= maxf2;  // That way, the number of evaluations is the same in all methods.
	fpointer = &(NLOPT_Optimization::evaluate_nlopt);
	objPointer = static_cast<void *>(rbobj);
	break;
      case BOBYQA:  // Pure local. No gradient
	opt = nlopt_create(NLOPT_LN_BOBYQA, n); 
	fpointer = &(NLOPT_Optimization::evaluate_nlopt);
	objPointer = static_cast<void *>(rbobj);
	break;
      case LBFGS:  // Pure local. Gradient based
	opt = nlopt_create(NLOPT_LD_LBFGS, n); 	
	fpointer = &(NLOPT_Optimization::evaluate_nlopt_grad);
	objPointer = static_cast<void *>(rgbobj);
	break;
      default: 
	FILE_LOG(logERROR) << "Algorithm not supported"; 
	return -1;
      }

    if (objPointer == NULL)
      {
	FILE_LOG(logERROR) << "Wrong object model (gradient/no gradient)"; 
	return -1;
      }

    nlopt_set_lower_bounds(opt, &mDown(0));
    nlopt_set_upper_bounds(opt, &mUp(0));
    nlopt_set_min_objective(opt, fpointer, objPointer);
    nlopt_set_maxeval(opt, maxf1) ;

    nlopt_result errortype = nlopt_optimize(opt, &Xnext(0), &fmin);
    checkNLOPTerror(errortype);

    if (maxf2)
      {
	nlopt_destroy(opt);  // Destroy previous one
	opt = nlopt_create(NLOPT_LN_BOBYQA, n); /* algorithm and dims */
	nlopt_set_lower_bounds(opt, &mDown(0));
	nlopt_set_upper_bounds(opt, &mUp(0));
	nlopt_set_min_objective(opt, fpointer, objPointer);
	nlopt_set_maxeval(opt, maxf2) ;
	
	errortype = nlopt_optimize(opt, &Xnext(0), &fmin);
	checkNLOPTerror(errortype);
      }
      
    nlopt_destroy(opt);  // Destroy opt

    ierror = static_cast<int>(errortype);
    return ierror;
  } // innerOptimize (uBlas)

  double NLOPT_Optimization::evaluate_nlopt (unsigned int n, const double *x,
					     double *grad, void *my_func_data)

  {
    vectord vx(n);
    std::copy(x,x+n,vx.begin());

    void *objPointer = my_func_data;
    RBOptimizableWrapper* OPTIMIZER = static_cast<RBOptimizableWrapper*>(objPointer);
    
    return OPTIMIZER->evaluate(vx);
  } /* evaluate_criteria_nlopt */


  double NLOPT_Optimization::evaluate_nlopt_grad (unsigned int n, const double *x,
						  double *grad, void *my_func_data)

  {
    vectord vx(n);
    std::copy(x,x+n,vx.begin());
    
    void *objPointer = my_func_data;
    RGBOptimizableWrapper* OPTIMIZER = static_cast<RGBOptimizableWrapper*>(objPointer);
    

    vectord vgrad = zvectord(n);
    double f =  OPTIMIZER->evaluate(vx,vgrad);
    if (grad && n)  std::copy(vgrad.begin(),vgrad.end(),grad);

    return f;
  } /* evaluate_criteria_nlopt */



}// namespace bayesopt

