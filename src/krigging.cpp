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


#include <cmath>
#include <ctime>

//#include <unistd.h>
#include "krigging.hpp"
#include "ublasinv.hpp"
#include "lhs.hpp"

#include "kernels.hpp"
#include "meanfuncs.hpp"

#include "krig_config.h"

#ifdef USE_DIRECT_FORTRAN
  #include "direct.hpp"
#else
  // NLOPT
  #include <nlopt.h>
  #include "nloptwpr.hpp"
#endif
  

GaussianProcess::GaussianProcess():
  mTheta(KERNEL_THETA), mP(KERNEL_P),
  mAlpha(PRIOR_ALPHA), mBeta (PRIOR_BETA),
  mDelta2(PRIOR_DELTA_SQ), mRegularizer(DEF_REGULARIZER),
  mMinIndex(0), mMaxIndex(0), mVerbose(0)
{} // Default constructor

GaussianProcess::GaussianProcess( double theta, double p,
				  double alpha, double beta, 
				  double delta, double noise):
  mTheta(theta), mP(p),
  mAlpha(alpha), mBeta (beta),
  mDelta2(delta), mRegularizer(noise),
  mMinIndex(0), mMaxIndex(0), mVerbose(0)
{} // Constructor

GaussianProcess::GaussianProcess( gp_params params ):
  mTheta(params.theta), mP(params.p),
  mAlpha(params.alpha), mBeta(params.beta),
  mDelta2(params.delta), mRegularizer(params.noise),
  mMinIndex(0), mMaxIndex(0), mVerbose(0)
{} // Constructor


GaussianProcess::~GaussianProcess()
{} // Default destructor


int GaussianProcess::computeCorrMatrix()
{
  size_t nSamples = mGPX.size1();
  bool inversionFlag;
  matrix<double> corrMatrix(nSamples,nSamples);
  
  if ( (nSamples != mInvR.size1()) || (nSamples != mInvR.size2()) )
    mInvR.resize(nSamples,nSamples);
  
  for (size_t ii=0; ii< nSamples; ii++)
    {
      for (size_t jj=0; jj < ii; jj++)
	{
	  corrMatrix(ii,jj) = correlationFunction(row(mGPX,ii), row(mGPX,jj));
	  corrMatrix(jj,ii) = corrMatrix(ii,jj);
	}
      corrMatrix(ii,ii) = correlationFunction(row(mGPX,ii),
					      row(mGPX,ii)) + mRegularizer;
    }
  
  inversionFlag = InvertMatrix(corrMatrix,mInvR);
  if (inversionFlag == false)    return -2;

  return 1;
}

vector<double> GaussianProcess::computeCrossCorrelation(const vector<double> &query)
{
  vector<double> knx(mGPX.size1());

  for (size_t ii=0; ii<mGPX.size1(); ++ii)
    {
      knx(ii) = correlationFunction(row(mGPX,ii), query);
    }
  return knx;
}

int GaussianProcess::fitGP()
{
  /** Computes the GP based on mGPX
   *  This function is hightly inefficient O(N^3). Use it only at 
   *  the begining
   */
  size_t nSamples = mGPX.size1();
  for (size_t ii=0; ii<nSamples; ii++)
    checkBoundsY(ii);

  //  normalizeData();

  int error = computeCorrMatrix();

  if (error < 0)
    return error;

  return precomputeGPParams();
} // fitGP


int GaussianProcess::precomputeGPParams()
{
  size_t nSamples = mGPX.size1();
  vector<double> colU(nSamples);

  for (size_t ii=0; ii< nSamples; ii++) 
    colU(ii) = meanFunction(row(mGPX,ii));

  mUInvR = prod(colU,mInvR);
  mUInvRUDelta = inner_prod(mUInvR,colU) + 1/mDelta2;
  
  vector<double> YInvR(nSamples);
  double YInvRY;
  
  //Test: Replace mYNorm for mGPY

  mMu =  inner_prod(mUInvR,mGPY) / mUInvRUDelta;
  noalias(YInvR) = prod(mGPY,mInvR);
  YInvRY = inner_prod(YInvR,mGPY);
  
  mSig = (mBeta + YInvRY - mMu*mMu/mUInvRUDelta) / (mAlpha + (nSamples+1) + 2);
  
  scalar_vector<double> colMu(nSamples,mMu);
  mYUmu = mGPY - colMu;
  
  return 1;
}


int GaussianProcess::addNewPointToGP(const vector<double> &Xnew, 
				     double Ynew)
{
  /** Add new point efficiently using Matrix Decomposition Lemma
   *  for the inversion of the correlation matrix. Maybe it is faster
   *  to just construct and invert a new matrix each time.
   */   

  size_t nSamples = mGPX.size1();
  size_t XDim = mGPX.size2();
  size_t NewDim = Xnew.size();
  
  scalar_vector<double> colU(nSamples+1,1.0);
  vector<double> Li(nSamples);
  vector<double> wInvR(nSamples);
  double wInvRw;
  double selfCorrelation, Ni;
  
  if (XDim != NewDim)
    {
      std::cout << "Dimensional Error" << std::endl;
      return -1;
    }

  addSample(Xnew,Ynew);
  checkBoundsY(nSamples);
  //  normalizeData();

  if(mVerbose>1)
    std::cout << mGPY(nSamples) << " vs. " << mGPY(mMinIndex) << std::endl;
    
  vector<double> correlationNewValue = computeCrossCorrelation(Xnew);
  
  selfCorrelation = correlationFunction(Xnew, Xnew) + mRegularizer;
  
  noalias(wInvR) = prod(correlationNewValue,mInvR);
  wInvRw = inner_prod(wInvR,correlationNewValue);
  Ni = 1/(selfCorrelation + wInvRw);
  noalias(Li) = -Ni * wInvR;
  mInvR += outer_prod(Li,Li) / Ni;
  
  mInvR.resize(nSamples+1,nSamples+1);
  
  Li.resize(nSamples+1);
  Li(nSamples) = Ni;
  
  row(mInvR,nSamples) = Li;
  column(mInvR,nSamples) = Li;
  
  return precomputeGPParams();
} // addNewPointToGP



double GaussianProcess::correlationFunction( const vector<double> &x1, 
					     const vector<double> &x2 )
{
    
  /* double diff, sum = 0.0, prod = 1.0;

    vector<double> xdiff = x1 - x2;
	
    for (size_t ii = 0; ii<x1.size(); ++ii) 
      {
	diff  = fabs(xdiff(ii)) / mTheta;
	sum  += diff;
	prod *= 1 + diff;
      }

      return(prod * exp(-sum));*/
       double grad;
       double theta = 0.05;
       return kernels::SEIso(x1,x2,grad,theta);
}  // correlationFunction


double GaussianProcess::meanFunction( const vector<double> &x)
{
  return means::One(x);
}

int GaussianProcess::prediction( const vector<double> &query,
				 double& yPred, double& sPred)
{
  vector<double> rInvR(mGPX.size1());
  double kn;
  double uInvRr, rInvRr;

  vector<double> colR = computeCrossCorrelation(query);
  kn = correlationFunction(query, query);
  
  noalias(rInvR) = prod(colR,mInvR);	
  rInvRr = inner_prod(rInvR,colR);
  uInvRr = inner_prod(mUInvR,colR);
  
  yPred = mMu + inner_prod( rInvR, mYUmu );
  sPred = sqrt( mSig * (kn - rInvRr + (1.0 - uInvRr) * (1.0 - uInvRr) 
			/ mUInvRUDelta ) );

  return 1;
}
	






/////////////////////////////////////////////////////////////////////



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


int SKO::optimize( vector<double> &bestPoint,
		   randEngine& mtRandom)
{
  size_t dim = bestPoint.size();
  vector<double> lowerBound = zero_vector<double>(dim);
  vector<double> upperBound = scalar_vector<double>(dim,1.0);
  
  return optimize(bestPoint,lowerBound,upperBound,mtRandom);
}

int SKO::optimize( vector<double> &bestPoint,
		   vector<double> &lowerBound,
		   vector<double> &upperBound,
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

  vector<double> xNext(nDims);
  double yNext;
  size_t nLHSSamples = N_LHS_EVALS_PER_DIM * nDims;

  if (nLHSSamples > MAX_LHS_EVALUATIONS)
    nLHSSamples = MAX_LHS_EVALUATIONS;
  
  if (  ( mMaxIterations > (MAX_ITERATIONS - nLHSSamples) )  
	|| ( mMaxIterations <= 0) )
    mMaxIterations = MAX_ITERATIONS - nLHSSamples;

  std::cout << "Sampling initial points...";
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
/*
int SKO::allocateMatrices(size_t nSamples, size_t nDims)
{
  if ( ( nSamples != mGPX.size1() ) || ( nDims != mGPX.size2() ) )
    mGPX.resize(nSamples,nDims,false);
  
  if ( ( nSamples != mInvR.size1() ) || ( nSamples != mInvR.size2() ) )
    mInvR.resize(nSamples,nSamples,false);
  
  if ( nSamples != mGPY.size() )
    mGPY.resize(nSamples,false);
  
  return 1;
}
*/
	
int SKO::sampleInitialPoints( size_t nSamples, size_t nDims,
			      bool useLatinBox, randEngine& mtRandom )
{
  /** \brief Sample a set of points to initialize GP fit
   * Use pure random sampling or uniform Latin Hypercube sampling
   * as appeared in Jones EGO
   */
   
  matrix<double> xPoints(nSamples,nDims);
  vector<double> yPoints(nSamples);
  vector<double> sample(nDims);

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




double SKO::lowerConfidenceBound(const vector<double> &query)
{    
  bool reachable = checkReachability(query);
  if (!reachable)
    return 0.0;
  
  double yPred, sPred, result;

  mGP.prediction(query,yPred,sPred);
  result = yPred -  mLCBparam*sPred;

  return result;

}

double SKO::negativeExpectedImprovement(const vector<double> &query)
{
  bool reachable = checkReachability(query);
  if (!reachable)
    return 0.0;

  double yPred, yDiff, yNorm, sPred;
  double result;

  mGP.prediction(query,yPred,sPred);
  yDiff = mGP.getValueAtMinimum() - yPred; 
// Because data is normalized, therefore Y minimum is 0

  yNorm = yDiff / sPred;
  
  if (mG == 1)
    {
      result = -1.0 * ( yDiff * cdf(yNorm) + sPred * pdf(yNorm) );
    }
  else
    {
      double pdfD = pdf(yNorm);
      double Tm2 = cdf(yNorm);
      double Tm1 = pdfD;
      double fg = factorial(mG);
      double Tact;
      
      double sumEI = pow(yNorm,mG)*Tm2 - mG*pow(yNorm,mG-1)*Tm1;

      for (unsigned int ii = 2; ii < mG; ii++) 
	{
	  Tact = (ii-1)*Tm2 - pdfD*pow(yNorm,ii-1);
	  sumEI += pow(-1.0,ii)*(fg/(factorial(ii)*factorial(mG-ii)))*
	    pow(yNorm,mG-ii)*Tact;
	  
	  //roll-up
	  Tm2 = Tm1;   Tm1 = Tact;
	}
      result = -1.0 * pow(sPred,mG) * sumEI;
    }
  
  return result;
       
}  // negativeExpectedImprovement

int SKO::nextPoint(vector<double> &Xnext)
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


    int maxT = MAX_DIRECT_ITERATIONS;

    if (mUseEI)
      fpointer = &(DIRECT::negeiwrap_);
    else
      fpointer = &(DIRECT::lcbwrap_);

    DIRECT::direct(fpointer, x, &n, &fmin, l, u, 
		   &ierror, &maxf, &maxT, objPointer);	


#else /* USE_DIRECT_FORTRAN */
    double (*fpointer)(unsigned int, const double *, double *, void *);

    if (mUseEI)
      fpointer = &(NLOPT_WPR::negeiwrap_nlopt);
    else 
      fpointer = &(NLOPT_WPR::lcbwrap_nlopt);

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


double SKO::evaluateNormalizedSample( const vector<double> &query)
{ 
  vector<double> unnormalizedQuery = ublas_elementwise_prod(query,
							    mRangeBound);
  
  unnormalizedQuery = ublas_elementwise_add(unnormalizedQuery,
					    mLowerBound);
    
  return evaluateSample(unnormalizedQuery);
} // evaluateNormalizedSample
  


unsigned int SKO::factorial(unsigned int no, unsigned int a)
{
  // termination condition
  if (0 == no || 1 == no)
    return a;
  
  // Tail recursive call
  return factorial(no - 1, no * a);
} //factorial

double SKO::pdf(double x)
{
  return (1 / (sqrt(2 * M_PI)) * exp(-(x*x)/2));
} //pdf
  
double SKO::cdf(double x)
{
  /** \brief Abromowitz and Stegun approximation of Normal CDF
   * 
   * Extracted from 
   * http://www.sitmo.com/doc/Calculating_the_Cumulative_Normal_Distribution
   * The most used algorithm is algorithm 26.2.17 from Abromowitz and Stegun, 
   * Handbook of Mathematical Functions. It has a maximum absolute error of 7.5e^-8.
   * 
   */
	
  static const double b1 =  0.319381530;
  static const double b2 = -0.356563782;
  static const double b3 =  1.781477937;
  static const double b4 = -1.821255978;
  static const double b5 =  1.330274429;
  static const double p  =  0.2316419;
  static const double c  =  0.39894228;
  
  if(x >= 0.0) {
    double t = 1.0 / ( 1.0 + p * x );
    return (1.0 - c * exp( -x * x / 2.0 ) * t *
	    ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
  }
  else {
    double t = 1.0 / ( 1.0 - p * x );
    return ( c * exp( -x * x / 2.0 ) * t *
	     ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
  }
} // cdf
