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


//#include <cmath>
//#include <ctime>

#include "gaussprocess.hpp"
#include "ublasinv.hpp"
#include "kernels.hpp"
#include "meanfuncs.hpp"

  

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
  size_t nSamples = mGPXX.size();
  bool inversionFlag;
  matrix<double> corrMatrix(nSamples,nSamples);
  
  if ( (nSamples != mInvR.size1()) || (nSamples != mInvR.size2()) )
    mInvR.resize(nSamples,nSamples);
  
  for (size_t ii=0; ii< nSamples; ii++)
    {
      for (size_t jj=0; jj < ii; jj++)
	{
	  corrMatrix(ii,jj) = correlationFunction(mGPXX[ii], mGPXX[jj]);
	  corrMatrix(jj,ii) = corrMatrix(ii,jj);
	}
      corrMatrix(ii,ii) = correlationFunction(mGPXX[ii],
					      mGPXX[ii]) + mRegularizer;
    }
  
  inversionFlag = InvertMatrix(corrMatrix,mInvR);
  if (inversionFlag == false)    return -2;

  return 1;
}

vector<double> GaussianProcess::computeCrossCorrelation(
				     const vector<double> &query)
{
  vector<double> knx(mGPXX.size());

  //TODO: Replace by transform
  for (size_t ii=0; ii<mGPXX.size(); ++ii)
    {
      knx(ii) = correlationFunction(mGPXX[ii], query);
    }
  return knx;
}
/*
double GaussianProcess::negativeLogLikelihood(double param,
					      vector<double> &grad)
{
  vector<double> alpha = dot(mInvR,mGPY);
    
}
*/
int GaussianProcess::fitGP()
{
  /** Computes the GP based on mGPXX
   *  This function is hightly inefficient O(N^3). Use it only at 
   *  the begining
   */
  size_t nSamples = mGPXX.size();
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
  size_t nSamples = mGPXX.size();
  vector<double> colU(nSamples);

  //TODO: Replace by transform
  for (size_t ii=0; ii< nSamples; ii++) 
    colU(ii) = meanFunction(mGPXX[ii]);

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

  size_t nSamples = mGPXX.size();
  size_t XDim = mGPXX[1].size();
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
  
  //TODO: There must be a better way to do this.
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
  double grad;
  return kernels::SEIso(x1,x2,grad,mTheta);
}  // correlationFunction


double GaussianProcess::correlationFunction( const vector<double> &x1, 
					     const vector<double> &x2,
					     double param, double &grad)
{
  return kernels::SEIso(x1,x2,grad,param);
}  // correlationFunction


double GaussianProcess::meanFunction( const vector<double> &x)
{
  return means::One(x);
}

int GaussianProcess::prediction( const vector<double> &query,
				 double& yPred, double& sPred)
{
  vector<double> rInvR(mGPXX.size());
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
	

