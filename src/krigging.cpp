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


#include "krig_config.h"

#ifdef USE_DIRECT_FORTRAN
  #include "direct.hpp"
#else
  // NLOPT
  #include <nlopt.h>
  #include "nloptwpr.hpp"
#endif
  


Krigging::Krigging():
mTheta(KERNEL_THETA), mP(KERNEL_P),
mAlpha(PRIOR_ALPHA), mBeta (PRIOR_BETA),
mDelta2(PRIOR_DELTA_SQ), mRegularizer(DEF_REGULARIZER),
mMaxIterations(MAX_ITERATIONS), mMaxDim(MAX_DIM), 
mUseCool(false), mG(1), mMinY(0.0), mMaxY(0.0),
mVerbose(0),mUseEI(true)
{} // Default constructor


Krigging::Krigging( double theta, double p,
		    double alpha, double beta, 
		    double delta, double noise,
		    size_t nIter, bool useCool):
mTheta(theta), mP(p),
mAlpha(alpha), mBeta (beta),
mDelta2(delta), mRegularizer(noise),
mMaxIterations(nIter), mMaxDim(MAX_DIM),
mUseCool(useCool), mG(1), 
mMinY(0.0), mMaxY(0.0),
mVerbose(0), mUseEI(true)
{} // Constructor

Krigging::Krigging( gp_params params,
		    size_t nIter, bool useCool):
mTheta(params.theta), mP(params.p),
mAlpha(params.alpha), mBeta(params.beta),
mDelta2(params.delta), mRegularizer(params.noise),
mMaxIterations(nIter), mMaxDim(MAX_DIM),
mUseCool(useCool), mG(1), 
mMinY(0.0), mMaxY(0.0),
mVerbose(0), mUseEI(true)
{ } // Constructor


Krigging::~Krigging()
{} // Default destructor


int Krigging::optimize( vector<double> &bestPoint,
			randEngine& mtRandom)
{
  size_t dim = bestPoint.size();
  vector<double> lowerBound = zero_vector<double>(dim);
  vector<double> upperBound = scalar_vector<double>(dim,1.0);
  
  return optimize(bestPoint,lowerBound,upperBound,mtRandom);
}

int Krigging::optimize( vector<double> &bestPoint,
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
      std::cout << "This algorithm is efficient up to " << mMaxDim << " dimensions." << std::endl;
    }

  vector<double> xNext(nDims);
  size_t nLHSSamples = N_LHS_EVALS_PER_DIM * nDims;

  if (nLHSSamples > MAX_LHS_EVALUATIONS)
    nLHSSamples = MAX_LHS_EVALUATIONS;
  
  if (  ( mMaxIterations > (MAX_ITERATIONS - nLHSSamples) )  || ( mMaxIterations <= 0) )
    mMaxIterations = MAX_ITERATIONS - nLHSSamples;

  std::cout << "Sampling initial points...";
  sampleInitialPoints(nLHSSamples, nDims, true, mtRandom);
  std::cout << "DONE" << std::endl;

  for (size_t ii = 0; ii < mMaxIterations; ii++)
    {
      if (mUseCool)
	{ updateCoolingScheme(mMaxIterations,ii); }
 
      if(mVerbose >0)
	{  std::cout << "Iteration " << ii+1 << std::endl;}
      
      nextPoint(xNext);
    
      if(mVerbose >0)
	{ std::cout << "Trying: " << xNext << std::endl;
	  std::cout << "Best: " << mMinX << std::endl;}			    

      addNewPointToGP(xNext);     
    }

  bestPoint = mMinX;

  return 1;
} // optimize


int Krigging::updateCoolingScheme( size_t nTotalIterations, size_t nCurrentIteration)
{

  double iterPercentaje = static_cast<double>(nTotalIterations)/static_cast<double>(nCurrentIteration);

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

int Krigging::allocateMatrices(size_t nSamples, size_t nDims)
{
  if ( ( nSamples != mGPX.size1() ) || ( nDims != mGPX.size2() ) )
    mGPX.resize(nSamples,nDims,false);
  
  if ( ( nSamples != mInvR.size1() ) || ( nSamples != mInvR.size2() ) )
    mInvR.resize(nSamples,nSamples,false);
  
  if ( nSamples != mGPY.size() )
    mGPY.resize(nSamples,false);
  
  if ( nDims != mMinX.size() )
    mMinX.resize(nDims,false);

  return 1;
}

	
int Krigging::sampleInitialPoints( size_t nSamples, size_t nDims,
				   bool useLatinBox, randEngine& mtRandom )
{
  /** \brief Sample a set of points to initialize GP fit
   * Use pure random sampling or uniform Latin Hypercube sampling
   * as appeared in Jones EGO
   */
   
  vector<double> Xsample(nDims);
  
  allocateMatrices(nSamples, nDims);
  
  if (useLatinBox)
      lhs(mGPX, mtRandom);
  else
      uniformSampling(mGPX, mtRandom);

  for(size_t i = 0; i < nSamples; i++)
    {
      Xsample = row(mGPX,i);
      if(mVerbose >0)
	std::cout << Xsample << std::endl;
      mGPY(i) = evaluateNormalizedSample(Xsample);
      checkBoundsY(i,Xsample);
    }

  fitGP();
  return 1;
} // sampleInitialPoints

int Krigging::checkBoundsY( size_t i, 
			    const vector<double>& Xsample )
{
  if (i == 0)
    {
      mMinY = mGPY(i);
      mMinX = Xsample;
      mMaxY = mGPY(i);
    }
  else
    {
      if ( mMinY > mGPY(i) )
	{
	  mMinY = mGPY(i);
	  mMinX = Xsample;
	}
      if ( mMaxY < mGPY(i) )
	{
	  mMaxY = mGPY(i);
	}
    }
  
  return 1;
}

int Krigging::fitGP()
{
  /** Computes the GP based on mGPX and mGPY
   *  This function is hightly inefficient O(N^3). Use it only at 
   *  the beggining
   */
  size_t Xsamples = mGPX.size1();
  bool inversionFlag;
  scalar_vector<double> colU(Xsamples, 1.0);
  matrix<double> corrMatrix(Xsamples,Xsamples);
  
  if ( (Xsamples != mInvR.size1()) || (Xsamples != mInvR.size2()) )
    mInvR.resize(Xsamples,Xsamples);

  normalizeData();
  
  for (size_t ii=0; ii< Xsamples; ii++)
    {
      for (size_t jj=0; jj < ii; jj++)
	{
	  corrMatrix(ii,jj) = correlationFunction(row(mGPX,ii), row(mGPX,jj));
	  corrMatrix(jj,ii) = corrMatrix(ii,jj);
	}
      corrMatrix(ii,ii) = correlationFunction(row(mGPX,ii),row(mGPX,ii)) + mRegularizer;
    }
  
  inversionFlag = InvertMatrix(corrMatrix,mInvR);
  if (inversionFlag == false)    return -2;
  
  mUInvR = prod(colU,mInvR);
  mUInvRUDelta = inner_prod(mUInvR,colU) + 1/mDelta2;
  
  vector<double> YInvR(mGPX.size1());
  double YInvRY;
  
  mMu =  inner_prod(mUInvR,mYNorm) / mUInvRUDelta;
  noalias(YInvR) = prod(mYNorm,mInvR);
  YInvRY = inner_prod(YInvR,mYNorm);
  
  mSig = (mBeta + YInvRY - mMu*mMu/mUInvRUDelta) / (mAlpha + (Xsamples+1) + 2);
  
  scalar_vector<double> colMu(mGPX.size1(),mMu);
  mYUmu = mYNorm - colMu;
  
  return 1;
}



int Krigging::addNewPointToGP(const vector<double> &Xnew)
{
  /** Add new point efficiently using Matrix Decomposition Lemma
   *  for the inversion of the correlation matrix. Maybe it is faster
   *  to just construct and invert a new matrix each time.
   */   

  size_t Xsamples = mGPX.size1();
  size_t Xdim = mGPX.size2();
  size_t Vdim = Xnew.size();
  
  scalar_vector<double> colU(Xsamples+1,1.0);
  vector<double> correlationNewValue(Xsamples);
  vector<double> Li(Xsamples);
  vector<double> wInvR(Xsamples);
  double wInvRw;
  double selfCorrelation, Ni;
  
  if (Xdim != Vdim)
    {
      std::cout << "Dimensional Error" << std::endl;
      return -1;
    }
  else
    {
      mGPX.resize(Xsamples+1,Xdim);
      mGPY.resize(Xsamples+1);
      
      row(mGPX,Xsamples) = Xnew;
      std::cout << "Print In: " << Xnew << std::endl;
      mGPY(Xsamples) = evaluateNormalizedSample(Xnew);
      checkBoundsY(Xsamples,Xnew);
      if(mVerbose>1)
	{  std::cout << mGPY(Xsamples) << " vs. " << mMinY << std::endl; }
    }
  //fitGP();
  
  
  normalizeData();
  
  for (size_t ii=0; ii< Xsamples; ii++)
    {
      correlationNewValue(ii) = correlationFunction(row(mGPX,ii), Xnew);
    }
  
  selfCorrelation = correlationFunction(Xnew, Xnew) + mRegularizer;
  
  noalias(wInvR) = prod(correlationNewValue,mInvR);
  wInvRw = inner_prod(wInvR,correlationNewValue);
  Ni = 1/(selfCorrelation + wInvRw);
  noalias(Li) = -Ni * wInvR;
  mInvR += outer_prod(Li,Li) / Ni;
  
  mInvR.resize(Xsamples+1,Xsamples+1);
  
  Li.resize(Xsamples+1);
  Li(Xsamples) = Ni;
  
  row(mInvR,Xsamples) = Li;
  column(mInvR,Xsamples) = Li;
  
  mUInvR = prod(colU,mInvR);
  mUInvRUDelta = inner_prod(mUInvR,colU) + 1/mDelta2;
  
  vector<double> YInvR(mGPX.size1());
  double YInvRY;
  
  mMu =  inner_prod(mUInvR,mYNorm) / mUInvRUDelta;
  noalias(YInvR) = prod(mYNorm,mInvR);
  YInvRY = inner_prod(YInvR,mYNorm);
  
  mSig = (mBeta + YInvRY - mMu*mMu/mUInvRUDelta) / (mAlpha + (Xsamples+1) + 2);
  
  scalar_vector<double> colMu(mGPX.size1(),mMu);
  mYUmu = mYNorm - colMu;
  
  return 1;
} // addNewPoint

double Krigging::correlationFunction(const vector<double> &x1, const vector<double> &x2)
{
  /** \brief GP Kernel computation
   * Kernel correlation based on
   * Matern linear function
   */
  
  double diff, sum = 0.0, prod = 1.0;
  vector<double> xdiff = x1 - x2;
	
  for (size_t ii = 0; ii<x1.size(); ++ii) 
    {
      diff  = fabs(xdiff(ii)) / mTheta;
      sum  += diff;
      prod *= 1 + diff;
    }

  return(prod * exp(-sum));
}  // correlationFunction
	
double Krigging::lowerConfidenceBound(const vector<double> &query)
{  
  vector<double> colR(mGPX.size1());
  vector<double> rInvR(mGPX.size1());
  double kn, yPred, sPred;
  double uInvRr, rInvRr, result;
  
  bool reachable = checkReachability(query);

  if (!reachable)
    return 0.0;
  
  for (size_t ii=0; ii< mGPX.size1(); ++ii)
    {
      colR(ii) = correlationFunction(row(mGPX,ii), query);
    }
  
  kn = correlationFunction(query, query) + mRegularizer;
  
  noalias(rInvR) = prod(colR,mInvR);	
  rInvRr = inner_prod(rInvR,colR);
  uInvRr = inner_prod(mUInvR,colR);
  
  yPred = mMu + inner_prod( rInvR, mYUmu );
  sPred = sqrt( mSig * (kn - rInvRr + (1.0 - uInvRr) * (1.0 - uInvRr) / mUInvRUDelta ) );

  double alpha = 1.0;

  result = yPred - alpha*sPred;

  return result;

}

double Krigging::negativeExpectedImprovement(const vector<double> &query)
{
  vector<double> colR(mGPX.size1());
  vector<double> rInvR(mGPX.size1());
  double kn, yPred, yDiff, yNorm, sPred;
  double uInvRr, rInvRr, result;
  
  bool reachable = checkReachability(query);
  
  if (!reachable)
    return 0.0;
  
  for (size_t ii=0; ii< mGPX.size1(); ++ii)
    {
      colR(ii) = correlationFunction(row(mGPX,ii), query);
    }
  
  kn = correlationFunction(query, query) + mRegularizer;
  
  noalias(rInvR) = prod(colR,mInvR);	
  rInvRr = inner_prod(rInvR,colR);
  uInvRr = inner_prod(mUInvR,colR);
  
  yPred = mMu + inner_prod( rInvR, mYUmu );
  sPred = sqrt( mSig * (kn - rInvRr + (1.0 - uInvRr) * (1.0 - uInvRr) / mUInvRUDelta ) );
  
  yDiff = - yPred; // Because data is normalized, therefore Y minimum is 0
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
	  sumEI += pow(-1.0,ii)*(fg/(factorial(ii)*factorial(mG-ii)))*pow(yNorm,mG-ii)*Tact;
	  
	  //roll-up
	  Tm2 = Tm1;   Tm1 = Tact;
	}
      result = -1.0 * pow(sPred,mG) * sumEI;
    }
  
  return result;
       
}  // negativeExpectedImprovement

int Krigging::nextPoint(vector<double> &Xnext)
{   
    double u[128], x[128], l[128];
    void *objPointer = dynamic_cast<void *>(this);
    double fmin = 1;
    int maxf = MAX_DIRECT_EVALUATIONS;
    int n = static_cast<int>(mGPX.size2());
    int ierror;

    if (objPointer == 0)
      std::cout << "Error casting the current object!" << std::endl;

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

    DIRECT::direct(fpointer, x, &n, &fmin, l, u, &ierror, &maxf, &maxT, objPointer);	


#else
    double (*fpointer)(unsigned int, const double *, double *, void *);

    if (mUseEI)
      fpointer = &(NLOPT_WPR::negeiwrap_nlopt);
    else 
      fpointer = &(NLOPT_WPR::lcbwrap_nlopt);

    nlopt_opt opt;
    opt = nlopt_create(NLOPT_GN_ORIG_DIRECT_L, n); /* algorithm and dims */
    nlopt_set_lower_bounds(opt, l);
    nlopt_set_upper_bounds(opt, u);
    nlopt_set_min_objective(opt, fpointer, objPointer);
    nlopt_set_maxeval(opt, maxf);

    nlopt_result errortype = nlopt_optimize(opt, x, &fmin);
      

    if(mVerbose)
      {std::cout << "Error:" << errortype << std::endl;}
      /*nlopt_minimize(NLOPT_GN_ORIG_DIRECT_L,
					    n, fpointer, objPointer,
					    l, u, x, &fmin,
					    -HUGE_VAL, 0.0, 0.0, 0.0, 
					    NULL, maxf, 1000000.0);*/

    ierror = static_cast<int>(errortype);
#endif

    // There should be a clever way to do this.
    array_adaptor<double> shared(n, x);
    vector<double, array_adaptor<double> > Xshared(n, shared); 

    Xnext = Xshared;
    
    return ierror;

} // nextPoint


double Krigging::evaluateNormalizedSample( const vector<double> &query)
{ 
  vector<double> unnormalizedQuery = ublas_elementwise_prod(query,
							    mRangeBound);
  
  unnormalizedQuery = ublas_elementwise_add(unnormalizedQuery,
					    mLowerBound);
    
  return evaluateSample(unnormalizedQuery);
} // evaluateNormalizedSample
  

void Krigging::normalizeData()
{
  scalar_vector<double> MinYVec(mGPY.size(), mMinY);
  mYNorm = (mGPY - MinYVec) * ( 1/(mMaxY-mMinY) );
} //normalizeData



unsigned int Krigging::factorial(unsigned int no, unsigned int a)
{
  // termination condition
  if (0 == no || 1 == no)
    return a;
  
  // Tail recursive call
  return factorial(no - 1, no * a);
} //factorial

double Krigging::pdf(double x)
{
  return (1 / (sqrt(2 * M_PI)) * exp(-(x*x)/2));
} //pdf
  
double Krigging::cdf(double x)
{
  /** \brief Abromowitz and Stegun approximation of Normal CDF
   * 
   * Extracted from 
   * http://www.sitmo.com/doc/Calculating_the_Cumulative_Normal_Distribution
   * The most used algorithm is algorithm 26.2.17 from Abromowitz and Stegun, 
   * Handbook of Mathematical Functions. It has a maximum absolute error of 7.5e^-8.
   * 
   */
	
  const double b1 =  0.319381530;
  const double b2 = -0.356563782;
  const double b3 =  1.781477937;
  const double b4 = -1.821255978;
  const double b5 =  1.330274429;
  const double p  =  0.2316419;
  const double c  =  0.39894228;
  
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

