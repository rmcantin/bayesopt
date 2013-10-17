/** \file inneroptimization.hpp 
    \brief C++ wrapper of the NLOPT library */
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


#ifndef __INNEROPTIMIZATION_HPP__
#define __INNEROPTIMIZATION_HPP__

#include "dll_stuff.h"
#include "optimization.hpp"

namespace bayesopt {

  // We plan to add more in the future since nlopt actually support many of them
  typedef enum {
    DIRECT,    ///< Global optimization
    LBFGS,     ///< Local, derivative based
    BOBYQA,    ///< Local, derivative free
    COMBINED   ///< Global exploration, local refinement (hand tuned)
  } innerOptAlgorithms;


  class InnerOptimization: public Optimization
  {
  public:
    InnerOptimization();
    virtual ~InnerOptimization(){};

    /** Sets the optimization algorithm  */
    void setAlgorithm(innerOptAlgorithms newAlg)
    { alg = newAlg; }

    /** Sets the optimization algorithm  */
    void setMaxEvals(size_t meval)
    { maxEvals = meval; }

    /** 
     * Limits of the hypercube. 
     * Currently, it assumes that all dimensions have the same limits.
     * @param down 
     * @param up 
     */
    void setLimits(double down, double up)
    {
      mDown = down;   mUp = up;
    }

    /** 
     * Compute the inner optimization algorithm
     * @param Xnext input: initial guess, output: result
     * @return error_code
     */
    int run(vectord &Xnext);
    

    /** 
     * Dummy function to be overriden by the actual function to be
     * evaluated.  
     * Note: it is not pure virtual because we might want
     * to use the other evaluate method
     * @param query input point
     * @return function value at query point
     */
    virtual double evaluate(const vectord& query) 
    {return 0.0;};


    /** 
     * Dummy function to be overriden by the actual function to be evaluated
     * Note: it is not pure virtual because we might want
     * to use the other evaluate method
     * @param query input point
     * @param grad output gradient at query point
     * @return function value at query point
     */
    virtual double evaluate(const vectord& query, 
				 vectord& grad) 
    {return 0.0;};


  private:

    int send_to_nlopt_optimize(double* x, int n, void* objPointer);	

    innerOptAlgorithms alg;
    double mDown, mUp;
    size_t maxEvals;
  };

}//namespace bayesopt

#endif
