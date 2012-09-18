/** -*- c++ -*- \file bayesoptbase.hpp \brief Bayesian optimization module */
/*
-----------------------------------------------------------------------------
   This file is part of BayesOptimization, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOptimization is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOptimization is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with BayesOptimization.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/


#ifndef  _BAYESOPTBASE_HPP_
#define  _BAYESOPTBASE_HPP_

#include "ctypes.h"
#include "specialtypes.hpp"
#include "logger.hpp"

#include "nonparametricprocess.hpp"
#include "criteria.hpp"


/** \addtogroup BayesOptimization */
/*@{*/

/**
 * \brief Bayesian optimization using different non-parametric 
 * processes as surrogate functions (non-parametric processes). 
 */
class BayesOptBase : public Logger
{
 public:
  
  /** 
   * Constructor
   * 
   * @param validSet  Set of potential inputs
   * @param gp        Pointer to the surrogate model
   */
  BayesOptBase( bopt_params params,
		 bool uselogfile = false,
		 const char* logfilename = "bayesopt.log");

  /** 
   * Default destructor
   */
  virtual ~BayesOptBase();

  /** 
   * Execute the optimization process of the function defined in evaluateSample.
   * 
   * @see evaluateSample
   * @see checkReachability
   *
   * @param bestPoint returns point with the optimum value in a ublas::vector.
   * 
   * @return 0 if terminate successfully, any other value otherwise
   */
  virtual int optimize(vectord &bestPoint) = 0;
  


  /** 
   * Function that defines the actual mathematical function to be optimized.
   * Virtual function for polymorphism. This function must need to be modified 
   * according to the specific problem.
   *
   * @param query point to be evaluated. 
   * 
   * @return value of the function at the point evaluated.
   */
  virtual double evaluateSample( const vectord &query ) = 0;
  

  /** 
   * This function checks if the query is valid or not. It can be used 
   * to introduce arbitrary constrains. Since the Gaussian process 
   * assumes smoothness, constrains are managed by DIRECT, being highly
   * time consuming. If the constrain is very tricky, DIRECT will need
   * much more function evaluations.
   *
   * Note: This function is experimental. 
   * 
   * @param query point to be evaluated.
   * 
   * @return boolean value showing if the the function is valid at
   *         the query point or not.
   */ 
  virtual bool checkReachability( const vectord &query )
  { return true; };

protected:
  /** 
   * Set the surrogate function based on the current parameters
   * 
   * @return 0 if terminate successfully
   */
  int setSurrogateFunction();

  /** 
   * Sets the number of iterations
   */
  inline void setNumberIterations()
  {
    if ((mParameters.n_iterations <= 0) || (mParameters.n_iterations > MAX_ITERATIONS))
      mParameters.n_iterations = MAX_ITERATIONS;
  };

  /** 
   * Sets the number of initial samples.
   * If the number is not set properly n<=0, then the number of initial samples
   * is fixed 10% of the total budget
   */
  inline size_t setInitSet()
  { 
    if (mParameters.n_init_samples <= 0)
      return static_cast<size_t>(ceil(0.1*mParameters.n_iterations));
    else
      return mParameters.n_init_samples;
  };

  inline double evaluateCriteria( const vectord &query )
  {
    bool reachable = checkReachability(query);
    if (!reachable)  return 0.0;
    return crit.evaluate(mGP,query);       
  }  // evaluateCriteria

  virtual int findOptimal(vectord &xOpt) = 0;

  int nextPoint( vectord &Xnext );  

protected:

  NonParametricProcess* mGP;        ///< Pointer to surrogate model
  Criteria crit;                    ///< Criteria model
  bopt_params mParameters;
};

/**@}*/
// end namespaces


#endif
