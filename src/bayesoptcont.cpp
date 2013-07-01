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

#include <limits>
#include "lhs.hpp"
#include "randgen.hpp"
#include "log.hpp"
#include "bayesoptcont.hpp"

namespace bayesopt  {
  
  ContinuousModel::ContinuousModel():
    BayesOptBase(), mBB(NULL)
  { 
    setAlgorithm(DIRECT);
  } // Def Constructor

  ContinuousModel::ContinuousModel(size_t dim, bopt_params parameters):
    BayesOptBase(dim,parameters), mBB(NULL)
  { 
    setAlgorithm(DIRECT);
  } // Constructor

  ContinuousModel::~ContinuousModel()
  {
    if (mBB != NULL)
      delete mBB;
  } // Default destructor

  int ContinuousModel::initializeOptimization()
  {
    if (mBB == NULL)
      {
	vectord lowerBound = zvectord(mDims);
	vectord upperBound = svectord(mDims,1.0);
	mBB = new utils::BoundingBox<vectord>(lowerBound,upperBound);
      }
    sampleInitialPoints();
    return 0;
  }

  vectord ContinuousModel::getFinalResult()
  {
    return mBB->unnormalizeVector(mGP->getPointAtMinimum());
  }


  int ContinuousModel::setBoundingBox(const vectord &lowerBound,
			      const vectord &upperBound)
  {
    if (mBB != NULL)
      delete mBB;
    
    mBB = new utils::BoundingBox<vectord>(lowerBound,upperBound);
    
    FILE_LOG(logINFO) << "Bounds: ";
    FILE_LOG(logINFO) << lowerBound;
    FILE_LOG(logINFO) << upperBound;
    
    return 0;
  } //setBoundingBox



  int ContinuousModel::plotStepData(size_t iteration, const vectord& xNext,
			    double yNext)
  {
    FILE_LOG(logINFO) << "Iteration: " << iteration+1 << " of " 
		      << mParameters.n_iterations << " | Total samples: " 
		      << iteration+1+mParameters.n_init_samples ;
    FILE_LOG(logINFO) << "Query: " << mBB->unnormalizeVector(xNext); ;
    FILE_LOG(logINFO) << "Query outcome: " << yNext ;
    FILE_LOG(logINFO) << "Best query: " 
		      << mBB->unnormalizeVector(mGP->getPointAtMinimum()); 
    FILE_LOG(logINFO) << "Best outcome: " <<  mGP->getValueAtMinimum();
    
    return 0;
  } //plotStepData

  int ContinuousModel::sampleInitialPoints()
  {
    
    size_t nSamples = mParameters.n_init_samples;
    int useLatinBox = 2;
    
    matrixd xPoints(nSamples,mDims);
    vectord yPoints(nSamples);
    vectord sample(mDims);
    randEngine mtRandom;
    
    if (useLatinBox == 1)           utils::lhs(xPoints, mtRandom);
    else if (useLatinBox == 2)           utils::sobol(xPoints, 0);
    else                utils::uniformSampling(xPoints, mtRandom);
    
    for(size_t i = 0; i < nSamples; i++)
      {
	sample = row(xPoints,i);
	yPoints(i) = evaluateSampleInternal(sample);
      }
    
    mGP->setSamples(xPoints,yPoints);
    mGP->fitInitialSurrogate();
    
    // For logging purpose
    if(mParameters.verbose_level > 0)
      {
	FILE_LOG(logDEBUG) << "Initial points:" ;
	double ymin = (std::numeric_limits<double>::max)();
	for(size_t i = 0; i < nSamples; i++)
	  {
	    sample = row(xPoints,i);
	    FILE_LOG(logDEBUG) << "Normalized X:" << sample ;
	    
	    if (yPoints(i)<ymin)  ymin = yPoints(i);
	      
	    FILE_LOG(logDEBUG) << "X:" << mBB->unnormalizeVector(sample) 
			       << "|Y:" << yPoints(i) << "|Min:" << ymin ;
	  }  
      }
    return 0;
  } // sampleInitialPoints

}  //namespace bayesopt
