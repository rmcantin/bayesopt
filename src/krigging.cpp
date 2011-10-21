//
// C++ Implementation: krigging
//
// Description: 
//
//
// Author: Ruben Martinez-Cantin  <rmcantin@unizar.es>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "krigging.hpp"
#include "lhs.hpp"


SKO::SKO():
  mGP(),
  mMaxIterations(MAX_ITERATIONS), mMaxDim(MAX_DIM), 
  mVerbose(0)
{} // Default constructor


SKO::SKO( double theta, double p,
	  double alpha, double beta, 
	  double delta, double noise,
	  size_t nIter, bool useCool):
  mGP(theta,noise),
  mMaxIterations(nIter), mMaxDim(MAX_DIM),
  mVerbose(0)
{} // Constructor

SKO::SKO( gp_params params,
	  size_t nIter, bool useCool):
  mGP(params.theta,params.noise),
  mMaxIterations(nIter), mMaxDim(MAX_DIM),
  mVerbose(0)
{ } // Constructor


SKO::~SKO()
{} // Default destructor


int SKO::optimize( vectord &bestPoint,
		   randEngine& mtRandom)
{
  size_t dim = bestPoint.size();
  vectord lowerBound = zvectord(dim);
  vectord upperBound = svectord(dim,1.0);
  
  return optimize(bestPoint,lowerBound,upperBound,mtRandom);
}

int SKO::optimize( vectord &bestPoint,
		   vectord &lowerBound,
		   vectord &upperBound,
		   randEngine& mtRandom)
{
  mVerbose = 1;
 
  mLowerBound = lowerBound;
  mRangeBound = upperBound - lowerBound;

  size_t nDims = bestPoint.size();
  
  if ( nDims > mMaxDim )
    { 
      std::cout << "Warning!! Too many dimensions!! " << std::endl;
      std::cout << "This algorithm is efficient up to " << mMaxDim 
		<< " dimensions." << std::endl;
    }

  vectord xNext(nDims);
  double yNext;
  size_t nLHSSamples = N_LHS_EVALS_PER_DIM * nDims;

  if (nLHSSamples > MAX_LHS_EVALUATIONS)
    nLHSSamples = MAX_LHS_EVALUATIONS;
  
  if (  ( mMaxIterations > (MAX_ITERATIONS - nLHSSamples) )  
	|| ( mMaxIterations <= 0) )
    mMaxIterations = MAX_ITERATIONS - nLHSSamples;

  std::cout << "Sampling initial points..." << std::endl;
  sampleInitialPoints(nLHSSamples, nDims, true, mtRandom);
  std::cout << "DONE" << std::endl;

  for (size_t ii = 0; ii < mMaxIterations; ii++)
    {      
      nextPoint(xNext);
    
      if(mVerbose >0)
	{ 
	  std::cout << "Iteration " << ii+1 << std::endl;
	  std::cout << "Trying: " << xNext << std::endl;
	  std::cout << "Best: " << mGP.getPointAtMinimum() << std::endl; 
	  std::cout << "Best outcome: " <<  mGP.getValueAtMinimum() <<  std::endl; 
	}
     
      yNext = evaluateNormalizedSample(xNext);
      mGP.addNewPointToGP(xNext,yNext);     
    }

  bestPoint = mGP.getPointAtMinimum();

  return 1;
} // optimize

int SKO::sampleInitialPoints( size_t nSamples, size_t nDims,
			      bool useLatinBox, randEngine& mtRandom )
{
  /** \brief Sample a set of points to initialize GP fit
   * Use pure random sampling or uniform Latin Hypercube sampling
   * as appeared in Jones EGO
   */
   
  matrixd xPoints(nSamples,nDims);
  vectord yPoints(nSamples);
  vectord sample(nDims);

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

  mGP.setSamples(xPoints,yPoints);
  mGP.fitGP();

  return 1;
} // sampleInitialPoints


double SKO::evaluateCriteria(const vectord &query)
{
  bool reachable = checkReachability(query);
  if (!reachable)
    return 0.0;

  return crit.evaluate(mGP,query);
       
}  // evaluateCriteria


