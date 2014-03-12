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


#ifndef __NLOPT_OPTIMIZATION_HPP__
#define __NLOPT_OPTIMIZATION_HPP__

#include "optimizable.hpp"
#include "nlopt.hpp"


namespace bayesopt {

  // We plan to add more in the future since nlopt actually support many of them
  typedef enum {
    DIRECT,    ///< Global optimization
    LBFGS,     ///< Local, derivative based
    BOBYQA,    ///< Local, derivative free
    COMBINED   ///< Global exploration, local refinement (hand tuned)
  } OptAlgorithms;


  class NLOPT_Optimization: //: public Optimization
  {
  public:
    NLOPT_Optimization();
    virtual ~NLOPT_Optimization(){};

    /** Sets the optimization algorithm  */
    void setAlgorithm(innerOptAlgorithms newAlg);

    /** Sets the optimization algorithm  */
    void setMaxEvals(size_t meval);

    /** Limits of the hypercube. */
    void setLimits(const vectord& down, const vectord& up);

    /** Compute the inner optimization algorithm
     * @param Xnext input: initial guess, output: result
     * @return error_code
     */
    int run(vectord &Xnext);
    
    double evaluate(const std::vector<double> &x, std::vector<double> &grad, void *data);

    /** Dummy function to be overriden by the actual function to be
     * evaluated.  
     * Note: it is not pure virtual because we might want
     * to use the other evaluate method
     * @param query input point
     * @return function value at query point
     */
    //    virtual double evaluate(const vectord& query) = 0;

    /** Dummy function to be overriden by the actual function to be evaluated
     * Note: it is not pure virtual because we might want
     * to use the other evaluate method
     * @param query input point
     * @param grad output gradient at query point
     * @return function value at query point
     */
    // virtual double evaluate(const vectord& query, 
    // 				 vectord& grad)  {return 0.0;};


  private:

    //int send_to_nlopt_optimize(double* x, int n, void* objPointer);	

    RBOptimizable *rbobj;

    innerOptAlgorithms alg;
    vectord mDown;
    vectord mUp;
    size_t maxEvals;
  };


  inline void NLOPT_Optimization::setAlgorithm(innerOptAlgorithms newAlg)
  { alg = newAlg; }

  inline void NLOPT_Optimization::setMaxEvals(size_t meval)
  { maxEvals = meval; }

  inline void NLOPT_Optimization::setLimits(double down, double up)
  { mDown = down;   mUp = up; }

  inline double NLOPT_Optimization::evaluate(const std::vector<double> &x, std::vector<double> &grad, void *data)
  { 
    vectord query(x.size());
    std::copy(x.begin(),x.end(),query.begin);
    NLOPT_Optimization* ptr = static_cast<void *>(data)
    return ptr->rbobj->evaluate(query); 
  }

}//namespace bayesopt

#endif
