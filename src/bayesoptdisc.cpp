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


SKO_DISC::SKO_DISC( vecOfvec &validSet, sko_params parameters,
       bool uselogfile,
       const char* logfilename):
  Logger(uselogfile,logfilename),
  mInputSet(validSet), mGP(NULL)
{ 
  mParameters = parameters;
  setNumberIterations();
} // Constructor


SKO_DISC::~SKO_DISC()
{
   delete mGP;
} // Default destructor


int SKO_DISC::setSurrogateFunction()
{
  if (mGP != NULL)
    delete mGP;
 
  switch(mParameters.s_name)
    {
    case s_gaussianProcess: 
      mGP = new BasicGaussianProcess(mParameters.theta,mParameters.noise); break;

    case s_gaussianProcessHyperPriors: 
      mGP = new GaussianProcess(mParameters.theta,mParameters.noise,
			       mParameters.alpha,mParameters.beta,
			       mParameters.delta);  break;

    case s_studentTProcess:
      mGP = new StudentTProcess(mParameters.theta,mParameters.noise); break;

    default:
      std::cout << "Error: surrogate function not supported." << std::endl;
      return -1;
    }

  mGP->setKernel(mParameters.k_name);
  return 0;
}

int SKO_DISC::optimize( vectord &bestPoint )
{

  crit.resetHedgeValues();

  size_t nDims = mInputSet[0].size();

  vectord xNext(nDims);
  double yNext;

  size_t nRandomSamples = setInitSet();

  if (mParameters.verbose_level > 0) 
    mOutput << "Sampling initial points..." << std::endl;

  sampleRandomPoints(nRandomSamples);

  if (mParameters.verbose_level > 0) 
    mOutput << "DONE" << std::endl;

  for (size_t ii = 0; ii < mParameters.n_iterations; ii++)
    {      
      // FIXME: Find what is the next point.
      nextPoint(xNext);
    
      if(mParameters.verbose_level >0)
	{ 
	  mOutput << "Iteration: " << ii+1 << " of " << mParameters.n_iterations;
	  mOutput << " | Total samples: " << ii+1+nRandomSamples << std::endl;
	  mOutput << "Trying point at: " << xNext << std::endl;
	  mOutput << "Best found at: " << mGP->getPointAtMinimum() << std::endl; 
	  mOutput << "Best outcome: " <<  mGP->getValueAtMinimum() <<  std::endl; 
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
      vectord xPoint = perms[i];
      if(mParameters.verbose_level >0)
	mOutput << xPoint << std::endl;
      yPoint = evaluateSample(xPoint);
      mGP->addSample(xPoint,yPoint);
    }
  /*randInt sample(mtRandom, intUniformDist(0,mInputSet.size()-1));

  for(size_t i = 0; i < nSamples; i++)
    {
      size_t index = sample();
      vectord xPoint = mInputSet[index];
      if(mParameters.verbose_level >0)
	mOutput << xPoint << std::endl;
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
  if (mParameters.c_name == c_gp_hedge)
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
      crit.setCriterium(mParameters.c_name);
      return findOptimal(Xnext);
    }
}



