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

#include "bayesopt.hpp"
#include "lhs.hpp"
#include "randgen.hpp"
#include "basicgaussprocess.hpp"


SKO_DISC::SKO_DISC( vecOfvec &validSet, NonParametricProcess* gp):
  mInputSet(validSet),
  mVerbose(0)
{
  mObservedNodes 
  crit_name = c_gp_hedge;
  if (gp == NULL)
    mGP = new BasicGaussianProcess(KERNEL_THETA,DEF_REGULARIZER);
  else
    mGP = gp;
} // Constructor

SKO_DISC::~SKO_DISC()
{} // Default destructor

int SKO_DISC::optimize( size_t &bestPointIndex, 
			size_t nIterations )
{
  mVerbose = 1;

  crit.resetHedgeValues();
  vectord xNext(nDims);
  double yNext;

  size_t nLHSSamples = std::min(N_LHS_EVALS_PER_DIM*nDims,
				MAX_LHS_EVALUATIONS);
  
  if (  ( nIterations > (MAX_ITERATIONS - nLHSSamples) )  
	|| ( nIterations <= 0) )
    nIterations = MAX_ITERATIONS - nLHSSamples;

  if (mVerbose > 0) std::cout << "Sampling initial points..." << std::endl;

  sampleInitialPoints(nLHSSamples);

  if (mVerbose > 0) std::cout << "DONE" << std::endl;

  for (size_t ii = 0; ii < nIterations; ii++)
    {      
      // Find what is the next point.
      size_t nextIndex = findNextIndex;
    
      if(mVerbose >0)
	{ 
	  std::cout << "Iteration " << ii+1 << "  |  ";
	  std::cout << "# of samples " << ii+1+nLHSSamples << std::endl;
	  std::cout << "Trying: " << xNext << std::endl;
	  std::cout << "Best: " << mGP->getPointAtMinimum() << std::endl; 
	  std::cout << "Best outcome: " <<  mGP->getValueAtMinimum() <<  std::endl; 
	}
     
      yNext = evaluateNormalizedSample(xNext);
      mGP->addNewPointToGP(xNext,yNext);     
    }

  bestPoint = mGP->getPointAtMinimum();

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
 
  if (useLatinBox)
      lhs(xPoints, mtRandom);
  else
      uniformSampling(xPoints, mtRandom);

  for(size_t i = 0; i < nSamples; i++)
    {
      sample = row(xPoints,i);
      if(mVerbose >0)
	std::cout << sample << std::endl;
      yPoints(i) = evaluateNormalizedSample(sample);
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
      double r_ei,r_lcb,r_poi,foo;

      crit.setCriterium(c_ei);
      innerOptimize(best_ei);
      mGP->prediction(best_ei,r_ei,foo);
      
      crit.setCriterium(c_lcb);
      innerOptimize(best_lcb);
      mGP->prediction(best_lcb,r_lcb,foo);

      crit.setCriterium(c_poi);
      innerOptimize(best_poi);
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
      return innerOptimize(Xnext);
    }
}



