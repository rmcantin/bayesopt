/** -*- c++ -*- \file bayesoptcont.hpp \brief Continuous Bayesian optimization */
/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2012 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
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

#ifndef  _BAYESOPTCONT_HPP_
#define  _BAYESOPTCONT_HPP_

#include "boundingbox.hpp"
#include "inneroptimization.hpp"
#include "bayesoptbase.hpp"


/** \addtogroup BayesOpt */
/**@{*/

/**
 * \brief Bayesian optimization using different non-parametric
 * processes as distributions over surrogate functions. The
 * exploration spaces is assumed to be continous and bounded.
 */
class BAYESOPT_API BayesOptContinuous: public InnerOptimization, 
			public BayesOptBase
{
 public:
   
  /** 
   * Default constructor
   */
  BayesOptContinuous();

  /** 
   * Constructor
   * @param params set of parameters (see parameters.h)
   */
  BayesOptContinuous( bopt_params params);

  /** 
   * Default destructor
   * 
   * @return 
   */
  virtual ~BayesOptContinuous();

  /** 
   * \brief Execute the optimization process of the function defined in
   * evaluateSample.  If no bounding box is defined, we assume that
   * the function is defined in the [0,1] hypercube, as a normalized
   * representation of the bound constrains.
   * 
   * @see scaleInput
   * @see evaluateSample
   *
   * @param bestPoint returns the optimum value in a ublas::vector.
   * @return 0 if terminate successfully, nonzero otherwise
   */
  int optimize(vectord &bestPoint);



  /** 
   * \brief Sets the bounding box. 
   *
   * @param lowerBound vector with the lower bounds of the hypercube 
   * @param upperBound vector with the upper bounds of the hypercube 
   * 
   * @return 0 if terminate successfully, nonzero otherwise
   */
  inline int setBoundingBox( const vectord &lowerBound,
			     const vectord &upperBound)
  {
    if (mBB != NULL)
      delete mBB;

    mBB = new BoundingBox<vectord>(lowerBound,upperBound);

    FILE_LOG(logINFO) << "Bounds: ";
    FILE_LOG(logINFO) << lowerBound;
    FILE_LOG(logINFO) << upperBound;

    return 0;
  };


protected:


  /** 
   * \brief Returns the corresponding criteria of a series of queries
   * in the hypercube [0,1] in order to choose the best point to try
   * the next iteration.
   * 
   * @param query point in the hypercube [0,1] to evaluate the
   * Gaussian process
   * 
   * @return negative criteria (Expected Improvement, LCB,
   * A-optimality, etc.).
   */	
  double innerEvaluate( const vectord &query )
  {
    return evaluateCriteria(query);       
  };  // evaluateCriteria

    
  /** 
   * \brief Wrapper for the target function normalize in the hypercube
   * [0,1]
   * @param query point to evaluate in [0,1] hypercube
   * @return actual return value of the target function
   */
  inline double evaluateNormalizedSample( const vectord &query )
  { 
    vectord unnormalizedQuery = mBB->unnormalizeVector(query);
    return evaluateSample(unnormalizedQuery);
  }; // evaluateNormalizedSample


  /** 
   * \brief Print data for every step according to the verbose level
   * 
   * @param iteration iteration number 
   * @param xNext next point
   * @param yNext function value at next point
   * 
   * @return error code
   */
  int plotStepData(size_t iteration, const vectord& xNext,
		   double yNext);

  /** \brief Sample a set of points to initialize the surrogate function.
   * It uses pure random sampling or uniform Latin Hypercube sampling.
   * @return error code
   */
  int sampleInitialPoints();

  /** 
   * \brief Wrapper of the innerOptimization class to find the optimal
   * point acording to the criteria.
   * 
   * @param xOpt optimal point
   * @return error code
   */
  inline int findOptimal(vectord &xOpt)
  { return innerOptimize(xOpt);};

protected:

  BoundingBox<vectord> *mBB;      ///< Bounding Box (input space limits)

};

/**@}*/


#endif
