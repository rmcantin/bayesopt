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
#include <ctime>
#include <limits>
#include <algorithm>
#include "bayesoptcont.hpp"
#include "lhs.hpp"
#include "randgen.hpp"
#include "basicgaussprocess.hpp"


SKO::SKO( NonParametricProcess* gp):
  InnerOptimization(), mVerbose(0), mLogFile()
{ 
  setAlgorithm(direct);
  crit_name = c_gp_hedge;
  if (gp == NULL)
    mGP = new BasicGaussianProcess(KERNEL_THETA,DEF_REGULARIZER);
  else
    mGP = gp;
} // Constructor

SKO::~SKO()
{} // Default destructor

int SKO::optimize( vectord &bestPoint,
		   vectord &lowerBound,
		   vectord &upperBound,
		   size_t nIterations )
{
  mVerbose = 1;

  if (mVerbose < 0)
    {
      mLogFile.open("log_bopt.out");
      if ( !mLogFile.is_open() )  mVerbose = 1;
    }

  crit.resetHedgeValues();
 
  mLowerBound = lowerBound;
  mRangeBound = upperBound - lowerBound;

  if (mVerbose > 0)
    {
      std::cout << "Bounds: "<< std::endl;
      std::cout << lowerBound << std::endl;
      std::cout << upperBound << std::endl;
    }


  size_t nDims = bestPoint.size();
 
  vectord xNext(nDims);
  double yNext;

  // Configuration simplified.
  // The number of initial samples is fixed 10% of the total budget
  if (nIterations <= 0) 
    nIterations = MAX_ITERATIONS;

  size_t nLHSSamples;

  if (nInitSet <= 0)
    nLHSSamples = static_cast<size_t>(ceil(0.1*nIterations));
  else
    nLHSSamples = nInitSet;

  if (mVerbose > 0) std::cout << "Sampling initial points..." << std::endl;

  clock_t start_time = clock();
  sampleInitialPoints(nLHSSamples, nDims, true);

  if (mVerbose > 0) std::cout << "DONE" << std::endl;

  for (size_t ii = 0; ii < nIterations; ii++)
    {      
      // Find what is the next point.
      nextPoint(xNext);
      
      if(mVerbose >0)
	{ 
	  vectord xScaled = unnormalizeVector(xNext);
	  vectord xOpt = unnormalizeVector(mGP->getPointAtMinimum());
	  std::cout << "Iteration: " << ii+1 << " of " << nIterations;
	  std::cout << " | Total samples: " << ii+1+nLHSSamples << std::endl;
	  std::cout << "Trying point at: " << xScaled << std::endl;
	  std::cout << "Best found at: " << xOpt << std::endl; 
	  std::cout << "Best outcome: " <<  mGP->getValueAtMinimum() <<  std::endl; 
	}
      yNext = evaluateNormalizedSample(xNext);
      mGP->addNewPointToGP(xNext,yNext); 

      if(mVerbose<0)
	{
	  clock_t lap = clock();
	  double t_lap = static_cast<double>(lap-start_time)/static_cast<double>(CLOCKS_PER_SEC);
	  mLogFile << t_lap << "|" << mGP->getValueAtMinimum() << "|";
	  mLogFile << yNext << "|" << unnormalizeVector(xNext) << std::endl;
	}
         
    }

  bestPoint = unnormalizeVector(mGP->getPointAtMinimum());

  if(mVerbose < 0)
    mLogFile.close();

  return 1;
} // optimize

int SKO::sampleInitialPoints( size_t nSamples, size_t nDims,
			      bool useLatinBox)
{
  /** \brief Sample a set of points to initialize GP fit
   * Use pure random sampling or uniform Latin Hypercube sampling
   * as appeared in Jones EGO
   */
   

  matrixd xPoints(nSamples,nDims);
  vectord yPoints(nSamples);
  vectord sample(nDims);
  randEngine mtRandom(100u);
  clock_t start_time = clock();
 
  if (useLatinBox)
      lhs(xPoints, mtRandom);
  else
      uniformSampling(xPoints, mtRandom);

  double ymin = std::numeric_limits<double>::max();

  for(size_t i = 0; i < nSamples; i++)
    {
      sample = row(xPoints,i);
      if(mVerbose >0)
	std::cout << sample << std::endl;

      yPoints(i) = evaluateNormalizedSample(sample);

      if (mVerbose < 0)
	{ 
	  if(yPoints(i)<ymin) 
	    ymin = yPoints(i);
	  clock_t lap = clock();
	  double t_lap = static_cast<double>(lap-start_time)/static_cast<double>(CLOCKS_PER_SEC);
	  mLogFile << t_lap << "|" << ymin << "|";
	  mLogFile << yPoints(i) << "|" << unnormalizeVector(sample) << std::endl;
	}
    }

  mGP->setSamples(xPoints,yPoints);
  mGP->fitGP();

  return 1;
} // sampleInitialPoints


int SKO::nextPoint(vectord &Xnext)
{
  crit.resetAnnealValues();
  if (crit_name == c_gp_hedge)
    {
      vectord best_ei(Xnext);
      vectord best_lcb(Xnext);
      vectord best_poi(Xnext);
      double l_ei,l_lcb,l_poi,foo;

      crit.setCriterium(c_ei);
      innerOptimize(best_ei);
      mGP->prediction(best_ei,l_ei,foo);
      
      crit.setCriterium(c_lcb);
      innerOptimize(best_lcb);
      mGP->prediction(best_lcb,l_lcb,foo);

      crit.setCriterium(c_poi);
      innerOptimize(best_poi);
      mGP->prediction(best_poi,l_poi,foo);

      // Since we want to find the minimum, the predicted value is loss value, not a
      // reward value.
      criterium_name better = crit.update_hedge(l_ei,l_lcb,l_poi);
      switch(better)
	{
	case c_ei: 
	  Xnext = best_ei;
	  if (mVerbose > 0) std::cout << "EI used." << std::endl;
	  break;
	case c_lcb: 
	  Xnext = best_lcb; 
	  if (mVerbose > 0) std::cout << "LCB used." << std::endl;
	  break;
	case c_poi: 
	  Xnext = best_poi; 
	  if (mVerbose > 0) std::cout << "POI used." << std::endl;
	  break;
	default: return -1;
	}
      return 1;
    }
  else
    {
      crit.setCriterium(crit_name);
      return innerOptimize(Xnext);
    }
}



