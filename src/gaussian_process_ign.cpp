
#include "gaussian_process_ign.hpp"
#include "cholesky.hpp"
#include "trace_ublas.hpp"

using boost::numeric::ublas::inplace_solve;
using boost::numeric::ublas::lower_tag;
using boost::numeric::ublas::lower;
  
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
  inplace_solve(L,alphU,lower_tag());
  double eta = inner_prod(colU,alphU) + 1/mDelta2;
  
  vectord alphY(mGPY);
  inplace_solve(L,alphY,lower_tag());
  double mu     = inner_prod(colU,alphY) / eta;
  double YInvRY = inner_prod(mGPY,alphY);
    
  double sigma = (mBeta + YInvRY - mu*mu*eta) / (mAlpha + (n+1) + 2);
  
  svectord colMu(n,mu);
  vectord yumu = mGPY - colMu;
  
  alphY = yumu;
  inplace_solve(L,alphY,lower_tag());

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
  //vectord rInvR(n);
  double kn;
  double uInvRr, rInvRr;
  double meanf = meanFunction(query);

  vectord colR = computeCrossCorrelation(query);
  kn = (*mKernel)(query, query);
  
  svectord colMean(n, meanf * mMu);
  vectord invRy = mGPY - colMean;
  cholesky_solve(mL,invRy,lower());

  // vectord ymu = mGPY - colMean;
  // vectord rInvR = prod(colR,mInvR);

  vectord invRr(colR);
  cholesky_solve(mL,invRr,lower());

  //  noalias(rInvR) = prod(colR,mInvR);	
  rInvRr = inner_prod(colR,invRr);//inner_prod(rInvR,colR);//
  //TODO:
  uInvRr = inner_prod(mMeanV,invRr);
  
  yPred = meanf*mMu + inner_prod( colR, invRy );
  sPred = sqrt( mSig * (kn - rInvRr + (meanf - uInvRr) * (meanf - uInvRr) 
			/ mUInvRUDelta ) );

  return 1;
}
	

int GaussianProcessIGN::precomputeGPParams()
{
  size_t nSamples = mGPXX.size();
  mMeanV.resize(nSamples);
  //mMeanV(nSamples);
  //std::transform(mGPXX.begin(),mGPXX.end(),mMeanV.begin(),NonParametricProcess::meanFunction);

  //TODO: Replace by transform
  for (size_t ii=0; ii< nSamples; ++ii) 
    mMeanV(ii) = meanFunction(mGPXX[ii]);

  mAlphaV.resize(nSamples,false);
  mAlphaV = mGPY;
  cholesky_solve(mL,mAlphaV,lower());

  vectord alphaMean(mMeanV);
  inplace_solve(mL,alphaMean,lower_tag());
  mUInvRUDelta = inner_prod(alphaMean,alphaMean) + 1/mDelta2;

  mMu = inner_prod(mMeanV,mAlphaV) / mUInvRUDelta;
  double YInvRY = inner_prod(mGPY,mAlphaV);

  // mInvR.resize(nSamples,nSamples);
  // mInvR.assign(boost::numeric::ublas::identity_matrix<value_type>(nSamples));
  // cholesky_solve(mL,mInvR,lower);

  //  mUInvR = prod(colU,mInvR);
  // mUInvRUDelta = inner_prod(mUInvR,colU) + 1/mDelta2;
  
  /*  vectord YInvR(nSamples);
  double YInvRY;
  
  mMu =  inner_prod(mUInvR,mGPY) / mUInvRUDelta;
  
  noalias(YInvR) = prod(mGPY,mInvR);
  YInvRY = inner_prod(YInvR,mGPY);*/
  
  mSig = (mBeta + YInvRY - mMu*mMu*mUInvRUDelta) / (mAlpha + (nSamples+1) + 2);
  
  return 1;
}


