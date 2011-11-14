#include "gaussprocess.hpp"
#include "cholesky.hpp"
#include "trace.hpp"


  
GaussianProcess::GaussianProcess( double theta, double noise,
				  double alpha, double beta, 
				  double delta):
  BasicGaussianProcess(theta,noise),
  mAlpha(alpha), mBeta (beta),
  mDelta2(delta)
{
  setAlgorithm(bobyqa);
  setLimits(0.,100.);
}  // Constructor



GaussianProcess::~GaussianProcess()
{} // Default destructor




double GaussianProcess::negativeLogLikelihood(double& grad,
					      size_t index)
{
  matrixd K = computeCorrMatrix(mRegularizer,0);
  size_t n = K.size1();
  matrixd L(n,n);
  cholesky_decompose(K,L);

  vectord colU(n);

  //TODO: Replace by transform
  for (size_t ii=0; ii< n; ii++) 
    colU(ii) = meanFunction(mGPXX[ii]);

  vectord alphU(colU);
  boost::numeric::ublas::inplace_solve(L,alphU,boost::numeric::ublas::lower_tag());
  double eta = inner_prod(colU,alphU) + 1/mDelta2;
  
  vectord alphY(mGPY);
  boost::numeric::ublas::inplace_solve(L,alphY,boost::numeric::ublas::lower_tag());
  double mu     = inner_prod(colU,alphY) / eta;
  double YInvRY = inner_prod(mGPY,alphY);
    
  double sigma = (mBeta + YInvRY - mu*mu*eta) / (mAlpha + (n+1) + 2);
  
  svectord colMu(n,mu);
  vectord yumu = mGPY - colMu;
  
  alphY = yumu;
  boost::numeric::ublas::inplace_solve(L,alphY,boost::numeric::ublas::lower_tag());

  double lik1 = inner_prod(yumu,alphY) / (2*sigma); 
  double lik2 = trace(L) + 0.5*n*log(sigma) + n*0.91893853320467; //log(2*pi)/2

  return lik1 + lik2 + mBeta/2 * mTheta - (mAlpha+1) * log(mTheta);
}


int GaussianProcess::prediction( const vectord &query,
				 double& yPred, double& sPred)
{
  vectord rInvR(mGPXX.size());
  double kn;
  double uInvRr, rInvRr;
  double meanf = meanFunction(query);

  vectord colR = computeCrossCorrelation(query);
  kn = correlationFunction(query, query);
  
  noalias(rInvR) = prod(colR,mInvR);	
  rInvRr = inner_prod(rInvR,colR);
  uInvRr = inner_prod(mUInvR,colR);
  
  yPred = meanf*mMu + inner_prod( rInvR, mYUmu );
  sPred = sqrt( mSig * (kn - rInvRr + (meanf - uInvRr) * (meanf - uInvRr) 
			/ mUInvRUDelta ) );

  return 1;
}
	

int GaussianProcess::precomputeGPParams()
{
  size_t nSamples = mGPXX.size();
  vectord colU(nSamples);

  //TODO: Replace by transform
  for (size_t ii=0; ii< nSamples; ii++) 
    colU(ii) = meanFunction(mGPXX[ii]);

  mUInvR = prod(colU,mInvR);
  mUInvRUDelta = inner_prod(mUInvR,colU) + 1/mDelta2;
  
  vectord YInvR(nSamples);
  double YInvRY;
  
  mMu =  inner_prod(mUInvR,mGPY) / mUInvRUDelta;
  
  noalias(YInvR) = prod(mGPY,mInvR);
  YInvRY = inner_prod(YInvR,mGPY);
  
  mSig = (mBeta + YInvRY - mMu*mMu*mUInvRUDelta) / (mAlpha + (nSamples+1) + 2);
  
  svectord colMu(nSamples,mMu);
  mYUmu = mGPY - colMu;

  return 1;
}


