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
#include "optimizable.hpp"
//#include "optimization.hpp"

namespace bayesopt {

  // We plan to add more in the future since nlopt actually support many of them
  typedef enum {
    DIRECT,    ///< Global optimization
    LBFGS,     ///< Local, derivative based
    BOBYQA,    ///< Local, derivative free
    COMBINED   ///< Global exploration, local refinement (hand tuned)
  } innerOptAlgorithms;


  class RBOptimizableWrapper
  {
  public:
    explicit RBOptimizableWrapper(RBOptimizable* rbo): rbo_(rbo){};
    virtual ~RBOptimizableWrapper(){};
    double evaluate(const vectord& query){return rbo_->evaluate(query);}
  private:
    RBOptimizable* rbo_;
  };

  class RGBOptimizableWrapper
  {
  public:
    explicit RGBOptimizableWrapper(RGBOptimizable* rgbo): rgbo_(rgbo){};
    virtual ~RGBOptimizableWrapper(){};
    double evaluate(const vectord& query, vectord& grad){return rgbo_->evaluate(query,grad);}
  private:
    RGBOptimizable* rgbo_;
  };


  class NLOPT_Optimization //: public Optimization
  {
  public:
    NLOPT_Optimization(RBOptimizable* rbo, size_t dim);
    NLOPT_Optimization(RGBOptimizable* rgbo, size_t dim);
    virtual ~NLOPT_Optimization();

    /** Sets the optimization algorithm  */
    void setAlgorithm(innerOptAlgorithms newAlg);

    /** Sets the optimization algorithm  */
    void setMaxEvals(size_t meval);

    /** Limits of the hypercube. */
    void setLimits(const vectord& down, const vectord& up);

    /** Limits of the hypercube assuming that all dimensions have the same limits. */
    void setLimits(double down, double up);

    /** Compute the inner optimization algorithm
     * @param Xnext input: initial guess, output: result
     * @return error_code
     */
    int run(vectord &Xnext);

    /** 
     * Wrapper of inner optimization to be evaluated by NLOPT
     * 
     * @param n # of dimensions
     * @param x input point
     * @param grad (NOT USED. Only for compatibily with NLOPT template, see evaluate_nlopt_grad)
     * @param my_func_data pointer to the NLOPT_Optimization object
     * 
     * @return function evaluation
     */  
    static double evaluate_nlopt (unsigned int n, const double *x,
				  double *grad, void *my_func_data);

    /** 
     * Wrapper of inner optimization to be evaluated by NLOPT
     * 
     * @param n # of dimensions
     * @param x input point
     * @param grad returns gradient evaluation
     * @param my_func_data pointer to the NLOPT_Optimization object
     * 
     * @return function evaluation
     */  
    static double evaluate_nlopt_grad (unsigned int n, const double *x,
				       double *grad, void *my_func_data);

  private:
    RBOptimizableWrapper *rbobj;
    RGBOptimizableWrapper *rgbobj;

    innerOptAlgorithms alg;
    vectord mDown, mUp;
    size_t maxEvals;
  };


  inline void NLOPT_Optimization::setAlgorithm(innerOptAlgorithms newAlg)
  { alg = newAlg; }

  inline void NLOPT_Optimization::setMaxEvals(size_t meval)
  { maxEvals = meval; }

  inline void NLOPT_Optimization::setLimits(const vectord& down, const vectord& up)
  { mDown = down;   mUp = up; }

  inline void NLOPT_Optimization::setLimits(double down, double up)
  { 
    for(size_t i = 0; i<mDown.size();++i) 
      {
	mDown(i) = down; mUp(i) = up;
      }
  };
}//namespace bayesopt

#endif
