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

#include "lhs.hpp"
#include "randgen.hpp"
#include "bayesoptcont.hpp"



SKO_CONT::SKO_CONT( sko_params parameters,
       bool uselogfile,
       const char* logfilename):
  SKO_BASE(parameters,uselogfile,logfilename), 
  mBB(NULL)
{ 
  setNumberIterations();
  setAlgorithm(direct);
  setSurrogateFunction();
} // Constructor

SKO_CONT::~SKO_CONT()
{
  if (mBB != NULL)
    delete mBB;
} // Default destructor


int SKO_CONT::optimize( vectord &bestPoint)
{
  crit.resetHedgeValues();

  size_t nDims = bestPoint.size();

  if (mBB == NULL)
    {
      vectord lowerBound = zvectord(nDims);
      vectord upperBound = svectord(nDims,1.0);
      mBB = new BoundingBox<vectord>(lowerBound,upperBound);
    }
  
  vectord xNext(nDims);
  double yNext;

  size_t nLHSSamples = setInitSet();

  if (mParameters.verbose_level > 0) 
    mOutput << "Sampling initial points..." << std::endl;

  clock_t start_time = clock();
  sampleInitialPoints(nLHSSamples, nDims, true);

  if (mParameters.verbose_level > 0) 
    mOutput << "DONE" << std::endl;

  for (size_t ii = 0; ii < mParameters.n_iterations; ii++)
    {      
      // Find what is the next point.
      nextPoint(xNext);
      
      if(mParameters.verbose_level >0)
	{ 
	  vectord xScaled = mBB->unnormalizeVector(xNext);
	  vectord xOpt = mBB->unnormalizeVector(mGP->getPointAtMinimum());
	  mOutput << "Iteration: " << ii+1 << " of " << mParameters.n_iterations;
	  mOutput << " | Total samples: " << ii+1+nLHSSamples << std::endl;
	  mOutput << "Trying point at: " << xScaled << std::endl;
	  mOutput << "Best found at: " << xOpt << std::endl; 
	  mOutput << "Best outcome: " <<  mGP->getValueAtMinimum() <<  std::endl; 
	}

      yNext = evaluateNormalizedSample(xNext);
      mGP->addNewPointToGP(xNext,yNext); 

      if(mParameters.verbose_level>1)
	{
	  clock_t lap = clock();
	  double t_lap = static_cast<double>(lap-start_time)/CLOCKS_PER_SEC;
	  mOutput << t_lap << "|" << mGP->getValueAtMinimum() << "|";
	  mOutput << yNext << "|" << mBB->unnormalizeVector(xNext) << std::endl;
	}
         
    }

  bestPoint = mBB->unnormalizeVector(mGP->getPointAtMinimum());

  return 1;
} // optimize

int SKO_CONT::sampleInitialPoints( size_t nSamples, size_t nDims,
			      bool useLatinBox)
{
  /** \brief Sample a set of points to initialize GP fit
   * Use pure random sampling or uniform Latin Hypercube sampling
   * as appeared in Jones EGO
   */
   

  matrixd xPoints(nSamples,nDims);
  vectord yPoints(nSamples);
  vectord sample(nDims);
  randEngine mtRandom;
  clock_t start_time = clock();
 
  if (useLatinBox)
      lhs(xPoints, mtRandom);
  else
      uniformSampling(xPoints, mtRandom);

  double ymin = std::numeric_limits<double>::max();

  for(size_t i = 0; i < nSamples; i++)
    {
      sample = row(xPoints,i);
      yPoints(i) = evaluateNormalizedSample(sample);

      if(mParameters.verbose_level > 0)
	{
	  mOutput << sample << std::endl;

	  if (mParameters.verbose_level > 1)
	    { 
	      if(yPoints(i)<ymin) 
		ymin = yPoints(i);

	      clock_t lap = clock();
	      double t_lap = static_cast<double>(lap-start_time)/CLOCKS_PER_SEC;
	      mOutput << t_lap << "|" << ymin << "|";
	      mOutput << yPoints(i) << "|" 
		       << mBB->unnormalizeVector(sample) << std::endl;
	    }
	}
    }

  mGP->setSamples(xPoints,yPoints);
  mGP->fitGP();

  return 1;
} // sampleInitialPoints




