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

#include "lhs.hpp"
#include "randgen.hpp"
#include "bayesoptcont.hpp"



BayesOptContinuous::BayesOptContinuous( bopt_params parameters,
       bool uselogfile,
       const char* logfilename):
  BayesOptBase(parameters,uselogfile,logfilename), 
  mBB(NULL)
{ 
  setAlgorithm(DIRECT);
} // Constructor

BayesOptContinuous::~BayesOptContinuous()
{
  if (mBB != NULL)
    delete mBB;

} // Default destructor




int BayesOptContinuous::optimize(vectord &bestPoint)
{
  mDims = bestPoint.size();

  if (mBB == NULL)
    {
      vectord lowerBound = zvectord(mDims);
      vectord upperBound = svectord(mDims,1.0);
      mBB = new BoundingBox<vectord>(lowerBound,upperBound);
    }

  sampleInitialPoints();

  vectord xNext(mDims);
  for (size_t ii = 0; ii < mParameters.n_iterations; ii++)
    {      
      // Find what is the next point.
      nextPoint(xNext);
      double yNext = evaluateNormalizedSample(xNext);
      mGP->updateSurrogateModel(xNext,yNext); 
      plotStepData(ii,xNext,yNext);
    }

  bestPoint = mBB->unnormalizeVector(mGP->getPointAtMinimum());

  return 1;
} // optimize

int BayesOptContinuous::plotStepData(size_t iteration, const vectord& xNext,
				     double yNext)
{
  vectord xScaled = mBB->unnormalizeVector(xNext);
  vectord xOpt = mBB->unnormalizeVector(mGP->getPointAtMinimum());
  FILE_LOG(logINFO) << "Iteration: " << iteration+1 << " of " 
		       << mParameters.n_iterations << " | Total samples: " 
		       << iteration+1+mParameters.n_init_samples ;
  FILE_LOG(logINFO) << "Trying point at: " << xScaled ;
  FILE_LOG(logINFO) << "Current outcome: " << yNext ;
  FILE_LOG(logINFO) << "Best found at: " << xOpt ; 
  FILE_LOG(logINFO) << "Best outcome: " <<  mGP->getValueAtMinimum();

  return 1;
}

int BayesOptContinuous::sampleInitialPoints()
{
   
  size_t nSamples = mParameters.n_init_samples;
  bool useLatinBox = true;

  matrixd xPoints(nSamples,mDims);
  vectord yPoints(nSamples);
  vectord sample(mDims);
  randEngine mtRandom;

  if (useLatinBox)
      lhs(xPoints, mtRandom);
  else
      uniformSampling(xPoints, mtRandom);

  for(size_t i = 0; i < nSamples; i++)
    {
      sample = row(xPoints,i);
      yPoints(i) = evaluateNormalizedSample(sample);
    }

  mGP->setSamples(xPoints,yPoints);
  mGP->fitInitialSurrogate();

  // For logging purpose
  if(mParameters.verbose_level > 0)
    {
      FILE_LOG(logDEBUG) << "Initial points:" ;
      double ymin = std::numeric_limits<double>::max();
      for(size_t i = 0; i < nSamples; i++)
	{
	  sample = row(xPoints,i);
	  FILE_LOG(logDEBUG) << sample ;
	  
	  if (mParameters.verbose_level > 1)
	    { 
	      if(yPoints(i)<ymin) 
		ymin = yPoints(i);
	      
	      FILE_LOG(logDEBUG) << ymin << "|" << yPoints(i) << "|" 
				 << mBB->unnormalizeVector(sample) ;
	    }
	}  
    }
  return 1;
} // sampleInitialPoints




