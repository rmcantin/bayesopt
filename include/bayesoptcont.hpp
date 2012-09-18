/** -*- c++ -*- \file bayesoptcont.hpp \brief Continuous Bayesian optimization */
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

//#include <iostream>
//#include <fstream>

//TODO: Check if everything is needed
#include "boundingbox.hpp"
#include "inneroptimization.hpp"
#include "bayesoptbase.hpp"


/** \addtogroup BayesOptimization */
/**@{*/

/**
 * \brief Sequential Kriging Optimization using different non-parametric 
 * processes as surrogate (kriging) functions. 
 */
class BayesOptContinuous: public InnerOptimization, BayesOptBase
{
 public:
  
  /** 
   * Constructor
   * 
   * @param gp        Pointer to the surrogate model
   */
  BayesOptContinuous( bopt_params parameters,
       bool uselogfile = false,
       const char* logfilename = "bayesopt.log"); 

  /** 
   * Default destructor
   * 
   * @return 
   */
  virtual ~BayesOptContinuous();

  /** 
   * Execute the optimization process of the function defined in evaluateSample.
   * If no bounding box is defined, we assume that the function is defined in the 
   * [0,1] hypercube, as a normalized representation of the bound constrains.
   * 
   * @see scaleInput
   * @see evaluateSample
   *
   * @param bestPoint returns the optimum value in a ublas::vector, it might also
   * be used as an initial point
   * 
   * @return 0 if terminate successfully, nonzero otherwise
   */
  int optimize(vectord &bestPoint);



  /** 
   * Sets the bounding box. 
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

    if (mParameters.verbose_level > 1)
      {
	mOutput << "Bounds: "<< std::endl;
	mOutput << lowerBound << std::endl;
	mOutput << upperBound << std::endl;
      }
    return 0;
  };


protected:


  /** 
   * Function that returns the corresponding criteria of a series 
   * of queries in the hypercube [0,1] in order to choose the best point to try
   * the next iteration.
   * 
   * @param query point in the hypercube [0,1] to evaluate the Gaussian process
   * 
   * @return negative criteria (Expected Improvement, LCB, A-optimality, etc.).
   */	
  double innerEvaluate( const vectord &query )
  {
    return evaluateCriteria(query);       
  };  // evaluateCriteria

    

  inline double evaluateNormalizedSample( const vectord &query)
  { 
    vectord unnormalizedQuery = mBB->unnormalizeVector(query);
    return evaluateSample(unnormalizedQuery);
  }; // evaluateNormalizedSample

  int sampleInitialPoints( size_t nSamples, size_t nDims, bool useLatinBox );

  inline int findOptimal(vectord &xOpt)
  { return innerOptimize(xOpt);};

protected:

  BoundingBox<vectord> *mBB;

};

/**@}*/


#endif
