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

#include "krig_config.h"

#ifdef USE_DIRECT_FORTRAN
  #include "direct.hpp"
#else
  // NLOPT
  #include <nlopt.h>
  #include "nloptwpr.hpp"
#endif

#include "krigging.hpp"
#include "lhs.hpp"
#include "criteria.hpp"

SKO::SKO():
  mGP(),
  mMaxIterations(MAX_ITERATIONS), mMaxDim(MAX_DIM), 
  mUseCool(false), mG(1), 
  mLCBparam(1.0), mUseEI(true), mVerbose(0)
{} // Default constructor


SKO::SKO( double theta, double p,
	  double alpha, double beta, 
	  double delta, double noise,
	  size_t nIter, bool useCool):
  mGP(theta,p,alpha,beta,delta,noise),
  mMaxIterations(nIter), mMaxDim(MAX_DIM),
  mUseCool(useCool), mG(1), 
  mLCBparam(1.0),mUseEI(true), mVerbose(0)
{} // Constructor

SKO::SKO( gp_params params,
	  size_t nIter, bool useCool):
  mGP(params),
  mMaxIterations(nIter), mMaxDim(MAX_DIM),
  mUseCool(useCool), mG(1), 
  mLCBparam(1.0),mUseEI(true), mVerbose(0)
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
      if (mUseCool)
	updateCoolingScheme(mMaxIterations,ii);
 
      if(mVerbose >0)
	std::cout << "Iteration " << ii+1 << std::endl;
      
      nextPoint(xNext);
    
      if(mVerbose >0)
	{ 
	  std::cout << "Trying: " << xNext << std::endl;
	  std::cout << "Best: " << mGP.getPointAtMinimum() << 
	    mGP.getValueAtMinimum() <<  std::endl; 
	}
      yNext = evaluateNormalizedSample(xNext);
      mGP.addNewPointToGP(xNext,yNext);     
    }

  bestPoint = mGP.getPointAtMinimum();

  return 1;
} // optimize


int SKO::updateCoolingScheme( size_t nTotalIterations, 
			      size_t nCurrentIteration)
{

  double iterPercentaje = static_cast<double>(nTotalIterations) 
                        / static_cast<double>(nCurrentIteration);

  if (iterPercentaje < 0.2)
    mG = 20;
  else if (iterPercentaje < 0.3)
    mG = 10; 
  else if (iterPercentaje < 0.4)
    mG = 5;     
  else if (iterPercentaje < 0.5)
    mG = 2;    
  else if (iterPercentaje < 0.7)
    mG = 1;       

  return 1;
}
	
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

  double yPred, sPred;

  mGP.prediction(query,yPred,sPred);
  double yMin = mGP.getValueAtMinimum(); 

  if (mUseEI)
    return criteria::negativeExpectedImprovement(yPred,sPred,yMin,mG);
  else
    return criteria::lowerConfidenceBound(yPred,sPred,mLCBparam);
       
}  // evaluateCriteria

int SKO::nextPoint(vectord &Xnext)
{   
    double x[128];
    void *objPointer = dynamic_cast<void *>(this);
    int n = static_cast<int>(Xnext.size());
    int error;

    if (objPointer == 0)
      std::cout << "Error casting the current object!" << std::endl;

    error = nextPoint(x, n, objPointer);

    // There should be a clever way to do this.
    array_adaptor<double> shared(n, x);
    vector<double, array_adaptor<double> > Xshared(n, shared); 

    Xnext = Xshared;
    
    return error;
} // nextPoint (uBlas)

int SKO::nextPoint(double* x, int n, void* objPointer)
{
    double u[128], l[128];
    double fmin = 1;
    int maxf = MAX_DIRECT_EVALUATIONS;    
    int ierror;

    for (int i = 0; i < n; ++i) {
	l[i] = 0.;
	u[i] = 1.;
    }
 
#ifdef USE_DIRECT_FORTRAN

    int (*fpointer)(int *, double *, double *, 
		    int *, int *,int *, double *,
		    int *, char *, int *, int);
    fpointer = &(DIRECT::criteriawrap_);

    int maxT = MAX_DIRECT_ITERATIONS;
    DIRECT::direct(fpointer, x, &n, &fmin, l, u, 
		   &ierror, &maxf, &maxT, objPointer);	


#else /* USE_DIRECT_FORTRAN */
    double (*fpointer)(unsigned int, const double *, double *, void *);
    fpointer = &(NLOPT_WPR::evaluate_criteria_nlopt);

    double coef = 0.8;
    nlopt_opt opt;
    opt = nlopt_create(NLOPT_GN_ORIG_DIRECT_L, n); /* algorithm and dims */
    nlopt_set_lower_bounds(opt, l);
    nlopt_set_upper_bounds(opt, u);
    nlopt_set_min_objective(opt, fpointer, objPointer);
    nlopt_set_maxeval(opt, round(maxf*coef) ) ;

    nlopt_result errortype = nlopt_optimize(opt, x, &fmin);

    if (coef < 1) 
      {
	//opt = nlopt_create(NLOPT_LN_BOBYQA, n); /* algorithm and dims */
	opt = nlopt_create(NLOPT_LN_NELDERMEAD, n); /* algorithm and dims */
	nlopt_set_lower_bounds(opt, l);
	nlopt_set_upper_bounds(opt, u);
	nlopt_set_min_objective(opt, fpointer, objPointer);
	nlopt_set_maxeval(opt, maxf-round(maxf*coef));
	
	errortype = nlopt_optimize(opt, x, &fmin);
      }
      

    if(mVerbose)
      {std::cout << "Error:" << errortype << std::endl;}

    ierror = static_cast<int>(errortype);
#endif /* USE_DIRECT_FORTRAN */

    return ierror;

} // nextPoint (C array)


double SKO::evaluateNormalizedSample( const vectord &query)
{ 
  vectord unnormalizedQuery = ublas_elementwise_prod(query,
							    mRangeBound);
  
  unnormalizedQuery = ublas_elementwise_add(unnormalizedQuery,
					    mLowerBound);
    
  return evaluateSample(unnormalizedQuery);
} // evaluateNormalizedSample
  


