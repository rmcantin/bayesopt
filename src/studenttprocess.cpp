#include "studenttprocess.hpp"
#include "cholesky.hpp"
#include "trace.hpp"

  
StudentTProcess::StudentTProcess( double theta, double noise):
  NonParametricProcess(theta,noise)
{
  setAlgorithm(bobyqa);
  setLimits(0.,100.);
}  // Constructor



StudentTProcess::~StudentTProcess()
{} // Default destructor




double StudentTProcess::negativeLogLikelihood(double &grad,
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
  double eta = inner_prod(colU,alphU);
  
  vectord alphY(mGPY);
  boost::numeric::ublas::inplace_solve(L,alphY,boost::numeric::ublas::lower_tag());
  double mu     = inner_prod(colU,alphY) / eta;
  double YInvRY = inner_prod(mGPY,alphY);
    
  double sigma = (YInvRY - mu*mu*eta) / (n-1);

  double negloglik = 0.5*( (n-1)*log(sigma) + trace(L) + log(eta) );

  return negloglik;
}


int StudentTProcess::prediction( const vectord &query,
				 double& yPred, double& sPred)
{
  size_t n = mGPXX.size();
  vectord rInvR(n);
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

  return n-1;
}
	

int StudentTProcess::precomputeGPParams()
{
  size_t nSamples = mGPXX.size();
  vectord colU(nSamples);

  //TODO: Replace by transform
  for (size_t ii=0; ii< nSamples; ii++) 
    colU(ii) = meanFunction(mGPXX[ii]);

  mUInvR = prod(colU,mInvR);
  mUInvRUDelta = inner_prod(mUInvR,colU);
  
  vectord YInvR(nSamples);
  double YInvRY;
  
  mMu =  inner_prod(mUInvR,mGPY) / mUInvRUDelta;
  
  noalias(YInvR) = prod(mGPY,mInvR);
  YInvRY = inner_prod(YInvR,mGPY);
  
  mSig = (YInvRY - mMu*mMu*mUInvRUDelta) / (nSamples-1);
  
  svectord colMu(nSamples,mMu);
  mYUmu = mGPY - colMu;

  std::cout << "Mu, Sigma" << mMu <<", "<< mSig << std::endl;
  std::cout << "Eta" << mUInvRUDelta << std::endl;

  return 1;
}


