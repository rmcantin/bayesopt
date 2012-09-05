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


SKO::SKO( sko_params parameters,
       bool uselogfile,
       const char* logfilename):
  InnerOptimization(), 
  Logger(uselogfile,logfilename),
  mGP(NULL)
{ 
  mParameters = parameters;
  setNumberIterations();
  setAlgorithm(direct);
  setSurrogateFunction();
} // Constructor

SKO::~SKO()
{
  if (mGP != NULL)
    delete mGP;
} // Default destructor

int SKO::setSurrogateFunction()
{
  if (mGP != NULL)
    delete mGP;
 
  switch(mParameters.s_name)
    {
    case s_gaussianProcess: 
      mGP = new BasicGaussianProcess(mParameters.noise); break;

    case s_gaussianProcessHyperPriors: 
      mGP = new GaussianProcess(mParameters.noise, mParameters.alpha,
				mParameters.beta,mParameters.delta);  break;

    case s_studentTProcess:
      mGP = new StudentTProcess(mParameters.noise); break;

    default:
      std::cout << "Error: surrogate function not supported." << std::endl;
      return -1;
    }
  
  mGP->setKernel(mParameters.theta,mParameters.k_name);
  return 0;
}

int SKO::optimize( vectord &bestPoint,
		   vectord &lowerBound,
		   vectord &upperBound)
{
  crit.resetHedgeValues();
 
  mLowerBound = lowerBound;
  mRangeBound = upperBound - lowerBound;

  if (mParameters.verbose_level > 1)
    {
      mOutput << "Bounds: "<< std::endl;
      mOutput << lowerBound << std::endl;
      mOutput << upperBound << std::endl;
    }


  size_t nDims = bestPoint.size();
 
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
	  vectord xScaled = unnormalizeVector(xNext);
	  vectord xOpt = unnormalizeVector(mGP->getPointAtMinimum());
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
	  mOutput << yNext << "|" << unnormalizeVector(xNext) << std::endl;
	}
         
    }

  bestPoint = unnormalizeVector(mGP->getPointAtMinimum());

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
		       << unnormalizeVector(sample) << std::endl;
	    }
	}
    }

  mGP->setSamples(xPoints,yPoints);
  mGP->fitGP();

  return 1;
} // sampleInitialPoints


int SKO::nextPoint(vectord &Xnext)
{
  crit.resetAnnealValues();
  if (mParameters.c_name == c_gp_hedge)
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
	  if (mParameters.verbose_level > 0) mOutput << "EI used." << std::endl;
	  break;
	case c_lcb: 
	  Xnext = best_lcb; 
	  if (mParameters.verbose_level > 0) mOutput << "LCB used." << std::endl;
	  break;
	case c_poi: 
	  Xnext = best_poi; 
	  if (mParameters.verbose_level > 0) mOutput << "POI used." << std::endl;
	  break;
	default: return -1;
	}
      return 1;
    }
  else
    {
      crit.setCriterium(mParameters.c_name);
      return innerOptimize(Xnext);
    }
}



