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

#include "bayesoptdisc.hpp"
#include "randgen.hpp"
#include "lhs.hpp"
#include "basicgaussprocess.hpp"


SKO_DISC::SKO_DISC( vecOfvec &validSet, NonParametricProcess* gp ):
  mInputSet(validSet), mVerbose(0), mLogFile()
{
  //mObservedNodes = zvectord(validSet);
  crit_name = c_gp_hedge;

  if (gp == NULL) 
    mGP = new BasicGaussianProcess(KERNEL_THETA,DEF_REGULARIZER);
  else            
    mGP = gp;
} // Constructor


SKO_DISC::~SKO_DISC()
{} // Default destructor


int SKO_DISC::optimize( vectord &bestPoint, 
			size_t nIterations )
{
  mVerbose = 1;

  if (mVerbose < 0)
    {
      mLogFile.open("log_bopt.out");
      if ( !mLogFile.is_open() )  mVerbose = 1;
    }

  crit.resetHedgeValues();

  size_t nDims = mInputSet[0].size();

  vectord xNext(nDims);
  double yNext;

  // Configuration simplified.
  // The number of initial samples is fixed 10% of the total budget
  if (nIterations <= 0) 
    nIterations = MAX_ITERATIONS;

  size_t nRandomSamples;

  if (nInitSet <= 0)
    nRandomSamples = static_cast<size_t>(ceil(0.1*nIterations));
  else
    nRandomSamples = nInitSet;

  if (mVerbose > 0) std::cout << "Sampling initial points..." << std::endl;

  sampleRandomPoints(nRandomSamples);

  if (mVerbose > 0) std::cout << "DONE" << std::endl;

  for (size_t ii = 0; ii < nIterations; ii++)
    {      
      // FIXME: Find what is the next point.
      nextPoint(xNext);
    
      if(mVerbose >0)
	{ 
	  std::cout << "Iteration: " << ii+1 << " of " << nIterations;
	  std::cout << " | Total samples: " << ii+1+nRandomSamples << std::endl;
	  std::cout << "Trying point at: " << xNext << std::endl;
	  std::cout << "Best found at: " << mGP->getPointAtMinimum() << std::endl; 
	  std::cout << "Best outcome: " <<  mGP->getValueAtMinimum() <<  std::endl; 
	}
     
      yNext = evaluateSample(xNext);
      mGP->addNewPointToGP(xNext,yNext);     
    }

  bestPoint = mGP->getPointAtMinimum();

  return 1;
} // optimize

int SKO_DISC::sampleRandomPoints( size_t nSamples )
{
  /** \brief Sample a set of points to initialize GP fit
   */
  double yPoint;
  randEngine rng;
  vecOfvec perms = mInputSet;

  // By using random permutations, we guarantee that 
  // the same point is not selected twice
  randomPerms(perms,rng);
  
  for(size_t i = 0; i < nSamples; i++)
    {
      yPoint = evaluateSample(perms[i]);
      mGP->addSample(perms[i],yPoint);
    }
  /*randInt sample(mtRandom, intUniformDist(0,mInputSet.size()-1));

  for(size_t i = 0; i < nSamples; i++)
    {
      size_t index = sample();
      vectord xPoint = mInputSet[index];
      if(mVerbose >0)
	std::cout << xPoint << std::endl;
      yPoint = evaluateSample(xPoint);
      mGP->addSample(xPoint,yPoint);
    }
  */
  mGP->fitGP();

  return 1;
} // sampleInitialPoints


int SKO_DISC::findOptimal(vectord &xOpt)
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

int SKO_DISC::nextPoint(vectord &Xnext)
{
  crit.resetAnnealValues();
  if (crit_name == c_gp_hedge)
    {
      vectord best_ei(Xnext);
      vectord best_lcb(Xnext);
      vectord best_poi(Xnext);
      double r_ei,r_lcb,r_poi,foo;

      crit.setCriterium(c_ei);
      findOptimal(best_ei);
      mGP->prediction(best_ei,r_ei,foo);
      
      crit.setCriterium(c_lcb);
      findOptimal(best_lcb);
      mGP->prediction(best_lcb,r_lcb,foo);

      crit.setCriterium(c_poi);
      findOptimal(best_poi);
      mGP->prediction(best_poi,r_poi,foo);

      criterium_name better = crit.update_hedge(r_ei,r_lcb,r_poi);
      switch(better)
	{
	case c_ei: Xnext = best_ei; break;
	case c_lcb: Xnext = best_lcb; break;
	case c_poi: Xnext = best_poi; break;
	default: return -1;
	}
      return 1;
    }
  else
    {
      crit.setCriterium(crit_name);
      return findOptimal(Xnext);
    }
}



