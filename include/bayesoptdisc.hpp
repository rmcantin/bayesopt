/**
 * @file   bayesopt.hpp
 * @author Ruben Martinez-Cantin <rmcantin@ist.isr.utl.pt>
 * @date   Thu Mar 26 02:12:36 2009
 * 
 * @brief  Sequential Krigging Optimization (SKO) 
 * 
 * This file implements Sequential Krigging Optimization using different 
 * non-parametric processes as surrogate (krigging) functions.
 *
 * 
 * Copyright: See COPYING file that comes with this distribution
 */


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


#ifndef  _BAYESOPTDISC_HPP_
#define  _BAYESOPTDISC_HPP_

#include "specialtypes.hpp"
#include "elementwiseUblas.hpp"

#include "inneroptimization.hpp"
#include "nonparametricprocess.hpp"
#include "logger.hpp"

#include "ctypes.h"
#include "criteria.hpp"

// Included here to simplify the C++-API
#include "gaussprocess.hpp"
#include "basicgaussprocess.hpp"
#include "studenttprocess.hpp"


/** \addtogroup BayesOptimization */
/*@{*/

/**
 * \brief Sequential Kriging Optimization using different non-parametric 
 * processes as surrogate (kriging) functions. 
 */
class SKO_DISC : public Logger
{
 public:
  
  /** 
   * Constructor
   * 
   * @param validSet  Set of potential inputs
   * @param gp        Pointer to the surrogate model
   */
  SKO_DISC( vecOfvec &validSet,
	    sko_params params,
	    bool uselogfile = false,
	    const char* logfilename = "bayesopt.log");

  /** 
   * Default destructor
   * 
   * @return 
   */
  virtual ~SKO_DISC();

  /** 
   * Execute the optimization process of the function defined in evaluateSample.
   * We assume that the function is defined in the [0,1] hypercube, as a 
   * normalized representation of the bound constrains.
   * 
   * @see scaleInput
   * @see evaluateSample
   *
   * @param bestPoint returns the optimum value in a ublas::vector defined in 
   * the hypercube [0,1], it might also be used as an initial point
   * 
   * @return 1 if terminate successfully, 0 otherwise
   */
  inline int optimize( vectord &bestPoint);
  


  /** 
   * Function that defines the actual mathematical function to be optimized.
   * Virtual function for polymorphism. This function must need to be modified 
   * according to the specific problem.
   *
   * @param query point to be evaluated. 
   * 
   * @return value of the function at the point evaluated.
   */
  virtual double evaluateSample( const vectord &query ) 
  { return 0.0; }
  
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
   * Chooses which criterium to optimize in the inner loop.
   * 
   * @param c criterium name
   */
  inline void setNumberIterations()
  {
    if ((mParameters.n_iterations <= 0) || (mParameters.n_iterations > MAX_ITERATIONS))
      mParameters.n_iterations = MAX_ITERATIONS;
  };

  inline size_t setInitSet()
  {  
    // Configuration simplified.
    // The number of initial samples is fixed 10% of the total budget
    if (mParameters.n_init_samples <= 0)
      return static_cast<size_t>(ceil(0.1*mParameters.n_iterations));
    else
      return mParameters.n_init_samples;
  };

  inline double evaluateCriteria( const vectord &query )
  {
    bool reachable = checkReachability(query);
    if (!reachable)  return 0.0;
    return crit.evaluate(*mGP,query);       
  }  // evaluateCriteria

  int sampleRandomPoints( size_t nSamples );
  int findOptimal(vectord &xOpt);
  int nextPoint( vectord &Xnext );  

protected:

  vecOfvec mInputSet;               ///< List of input points
  NonParametricProcess* mGP;        ///< Pointer to surrogate model
  Criteria crit;                    ///< Criteria model
  sko_params mParameters;
};

/**@}*/
// end namespaces


#endif
