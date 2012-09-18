
#include "gaussian_process_ign.hpp"
#include "cholesky.hpp"
#include "trace_ublas.hpp"


  
GaussianProcessIGN::GaussianProcessIGN( double noise, double alpha, 
					double beta, double delta):
  GaussianProcess(noise),
  mAlpha(alpha), mBeta (beta), mDelta2(delta)
{}  // Constructor



GaussianProcessIGN::~GaussianProcessIGN()
{} // Default destructor




double GaussianProcessIGN::negativeLogLikelihood(size_t index)
{
  matrixd K = computeCorrMatrix();
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

  double th = mKernel->getScale(index);

  return lik1 + lik2 + mBeta/2 * th -
    (mAlpha+1) * log(th);
}


int GaussianProcessIGN::prediction( const vectord &query,
				 double& yPred, double& sPred)
{
  size_t n = mGPXX.size();
  vectord rInvR(n);
  double kn;
  double uInvRr, rInvRr;
  double meanf = meanFunction(query);

  vectord colR = computeCrossCorrelation(query);
  kn = (*mKernel)(query, query);
  
  svectord colMu(n,mMu);
  vectord yumu = mGPY - meanf*colMu;
  
  noalias(rInvR) = prod(colR,mInvR);	
  rInvRr = inner_prod(rInvR,colR);
  uInvRr = inner_prod(mUInvR,colR);
  
  yPred = meanf*mMu + inner_prod( rInvR, yumu );
  sPred = sqrt( mSig * (kn - rInvRr + (meanf - uInvRr) * (meanf - uInvRr) 
			/ mUInvRUDelta ) );

  return 1;
}
	

int GaussianProcessIGN::precomputeGPParams()
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
  
  return 1;
}


