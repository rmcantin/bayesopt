/**
 * @file   bayesoptcont.hpp
 * @brief  Sequential Krigging Optimization (SKO)
 * 
 * This file implements Sequential Krigging Optimization using different 
 * non-parametric processes as surrogate (krigging) functions.
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


#ifndef  _BAYESOPTCONT_HPP_
#define  _BAYESOPTCONT_HPP_

//TODO: Check if everything is needed
#include "specialtypes.hpp"
#include "elementwiseUblas.hpp"

#include "inneroptimization.hpp"
#include "nonparametricprocess.hpp"

#include "ctypes.h"
#include "criteria.hpp"

// Included here to simplify the C++-API
#include "gaussprocess.hpp"
#include "basicgaussprocess.hpp"
#include "studenttprocess.hpp"


/** \addtogroup BayesOptimization */
/**@{*/

/**
 * \brief Sequential Kriging Optimization using different non-parametric 
 * processes as surrogate (kriging) functions. 
 */
class SKO: public InnerOptimization
{
 public:
  
  /** 
   * Constructor
   * 
   * @param gp        Pointer to the surrogate model
   */
  SKO( NonParametricProcess* gp = NULL ); 

  /** 
   * Default destructor
   * 
   * @return 
   */
  virtual ~SKO();

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
   * @param nIterations number of iterations (budget)
   * 
   * @return 1 if terminate successfully, 0 otherwise
   */
  inline int optimize( vectord &bestPoint, 
		       size_t nIterations )
  {
    size_t dim = bestPoint.size();
    vectord lowerBound = zvectord(dim);
    vectord upperBound = svectord(dim,1.0);
  
    return optimize(bestPoint,lowerBound,upperBound,nIterations);
  }



  /** 
   * Execute the optimization process of the function defined in evaluateSample.
   * We assume that the function is defined in the hypercube defined by 
   * the lower and upper vectors. 
   *
   * @param bestPoint returns the optimum value in a ublas::vector x, 
   * it might also be used as an initial point
   * @param lowerBound vector with the lower bounds of the hypercube 
   * @param upperBound vector with the upper bounds of the hypercube 
   * @param nIterations number of iterations (budget)
   * 
   * @return 1 if terminate successfully, 0 otherwise
   */
  int optimize( vectord &bestPoint,
		vectord &lowerBound,
		vectord &upperBound,
		size_t nIterations );


  /** 
   * Chooses which criterium to optimize in the inner loop.
   * 
   * @param c criterium name
   */
  void setCriteria (criterium_name c)
  {crit_name = c;}

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
  { return 0.0; };
  
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

  /** 
   * Function that returns the corresponding criteria of a series 
   * of queries in the hypercube [0,1] in order to choose the best point to try
   * the next iteration.
   * 
   * @param query point in the hypercube [0,1] to evaluate the Gaussian process
   * 
   * @return negative criteria (Expected Improvement, LCB, A-optimality, etc.).
   */	

  virtual double innerEvaluate( const vectord &query, 
				vectord &grad )
  {   return evaluateCriteria(query); }

protected:

  inline double evaluateCriteria( const vectord &query )
  {
    bool reachable = checkReachability(query);
    if (!reachable)  return 0.0;
    return crit.evaluate(*mGP,query);       
  }  // evaluateCriteria

  inline vectord unnormalizeVector( const vectord &vin)
  {
    vectord vout = ublas_elementwise_prod(vin,mRangeBound);
    return ublas_elementwise_add(vout, mLowerBound);
  }  // unnormalizeVector
    

  inline double evaluateNormalizedSample( const vectord &query)
  { 
    vectord unnormalizedQuery = unnormalizeVector(query);
    return evaluateSample(unnormalizedQuery);
  } // evaluateNormalizedSample

  int sampleInitialPoints( size_t nSamples, size_t nDims, bool useLatinBox );
  int nextPoint( vectord &Xnext );  

protected:

  NonParametricProcess* mGP;        ///< Pointer to surrogate model
  Criteria crit;                    ///< Criteria model
  criterium_name crit_name;         ///< Name of the criteria
  size_t mMaxIterations;            ///< Maximum SKO evaluations (budget)
  vectord mLowerBound, mRangeBound; ///< Lower bound and range of the input space
  int mVerbose;                     ///< Verbose level

};

/**@}*/


#endif
