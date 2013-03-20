/**  \file nloptwpr.h \brief Simple wrapper to optimize object methods with NLOPT */
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
#ifndef _NLOPTWPR_HPP_
#define _NLOPTWPR_HPP_


namespace NLOPT_WPR
{
  extern "C" {

  /** 
   * Wrapper of inner optimization to be evaluated by NLOPT
   * 
   * @param n # of dimensions
   * @param x input point
   * @param grad (only for compatibily, not used or changed)
   * @param my_func_data pointer to the InnerOptimization object
   * 
   * @return function evaluation
   */  
  double evaluate_nlopt (unsigned int n, const double *x,
			 double *grad, void *my_func_data);

  /** 
   * Wrapper of inner optimization to be evaluated by NLOPT
   * 
   * @param n # of dimensions
   * @param x input point
   * @param grad returns gradient evaluation
   * @param my_func_data pointer to the InnerOptimization object
   * 
   * @return function evaluation
   */  
  double evaluate_nlopt_grad (unsigned int n, const double *x,
			      double *grad, void *my_func_data);

  
  }
}

#endif
