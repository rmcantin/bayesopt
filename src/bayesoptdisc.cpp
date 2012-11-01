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

#include "bayesoptdisc.hpp"
#include "randgen.hpp"
#include "lhs.hpp"


BayesOptDiscrete::BayesOptDiscrete( const vecOfvec &validSet):
  BayesOptBase(), mInputSet(validSet)
{} // Constructor


BayesOptDiscrete::BayesOptDiscrete( const vecOfvec &validSet, 
				    bopt_params parameters):
  BayesOptBase(parameters), mInputSet(validSet)
{} // Constructor


BayesOptDiscrete::~BayesOptDiscrete()
{} // Default destructor


int BayesOptDiscrete::optimize( vectord &bestPoint )
{
  mDims = mInputSet[0].size();
  sampleInitialPoints();

  vectord xNext(mDims);
  for (size_t ii = 0; ii < mParameters.n_iterations; ++ii)
    {      
      nextPoint(xNext);
      double yNext = evaluateSample(xNext);
      mGP->updateSurrogateModel(xNext,yNext);   
      plotStepData(ii,xNext,yNext);
    }

  bestPoint = mGP->getPointAtMinimum();

  return 1;
} // optimize


int BayesOptDiscrete::plotStepData(size_t iteration, const vectord& xNext,
				   double yNext)
{
  if(mParameters.verbose_level >0)
    { 
      FILE_LOG(logINFO) << "Iteration: " << iteration+1 << " of " 
			<< mParameters.n_iterations << " | Total samples: " 
			<< iteration+1+mParameters.n_init_samples ;
      FILE_LOG(logINFO) << "Trying point at: " << xNext ;
      FILE_LOG(logINFO) << "Current outcome: " << yNext ;
      FILE_LOG(logINFO) << "Best found at: " << mGP->getPointAtMinimum() ; 
      FILE_LOG(logINFO) << "Best outcome: " <<  mGP->getValueAtMinimum() ;    
    }

  return 1;
}


int BayesOptDiscrete::sampleInitialPoints()
{
  size_t nSamples = mParameters.n_init_samples;
  double yPoint;
  randEngine rng;
  vecOfvec perms = mInputSet;

  // By using random permutations, we guarantee that 
  // the same point is not selected twice
  randomPerms(perms,rng);

  vectord xPoint(mInputSet[0].size());
  for(size_t i = 0; i < nSamples; i++)
    {
      xPoint = perms[i];
      yPoint = evaluateSample(xPoint);
      mGP->addSample(xPoint,yPoint);
    }

  mGP->fitInitialSurrogate();

  // For logging purpose
  if(mParameters.verbose_level > 0)
    {
      FILE_LOG(logDEBUG) << "Initial points:" ;
      double ymin = std::numeric_limits<double>::max();
      for(size_t i = 0; i < nSamples; i++)
	{
	  yPoint = mGP->getSample(i,xPoint);
	  FILE_LOG(logDEBUG) << xPoint ;
	  
	  if (mParameters.verbose_level > 1)
	    { 
	      if(yPoint<ymin) 
		ymin = yPoint;
	      
	      FILE_LOG(logDEBUG) << ymin << "|" << yPoint ;
	    }
	}  
    }
  return 1;
} // sampleInitialPoints


int BayesOptDiscrete::findOptimal(vectord &xOpt)
{
  double current, min;
  
  xOpt = *mInputSet.begin();
  min = evaluateCriteria(xOpt);
  
  for(vecOfvecIterator it = mInputSet.begin();
      it != mInputSet.end(); ++it)
    {
      current = evaluateCriteria(*it);
      if (current < min)
	{
	  xOpt = *it;  
	  min = current;
	}
    }
  return 1;
}



