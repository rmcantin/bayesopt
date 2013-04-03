/*
-----------------------------------------------------------------------------
   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/

#include "gaussian_process_ml.hpp"
#include "cholesky.hpp"
#include "trace_ublas.hpp"
#include "gauss_distribution.hpp"

using boost::numeric::ublas::inplace_solve;
using boost::numeric::ublas::lower_tag;
using boost::numeric::ublas::lower;
  

GaussianProcessML::GaussianProcessIGN(double noise):
  NonParametricProcess(noise)
{
  d_.reset(new GaussianDistribution());
}  // Constructor



GaussianProcessIGN::~GaussianProcessIGN()
{} // Default destructor




double GaussianProcessIGN::negativeLogLikelihood()
{
  matrixd K = computeCorrMatrix();
  size_t n = K.size1();
  matrixd L(n,n);
  utils::cholesky_decompose(K,L);

  matrixd featAll = mMean->getFeatures(mGPXX);

  vectord alphU(mMeanV);
  inplace_solve(L,alphU,lower_tag());
  double eta = inner_prod(alphU,alphU) + 1/mDelta2;
  
  vectord alphY(mGPY);
  inplace_solve(L,alphY,lower_tag());
  double mu     = inner_prod(alphU,alphY) / eta;
  double YInvRY = inner_prod(alphY,alphY);
    
  double sigma = (mBeta + YInvRY - mu*mu*eta) / (mAlpha + (n+1) + 2);
  
  vectord yumu = mGPY - mMeanV*mu;
  alphY = yumu;
  inplace_solve(L,alphY,lower_tag());

  double lik1 = inner_prod(alphY,alphY) / (2*sigma); 
  double lik2 = log_trace(L) + 0.5*n*log(sigma) + n*0.91893853320467; //log(2*pi)/2
  
  //TODO: This must be wrong.
  size_t index = 1;
  vectord th = mKernel->getScale(index);

  return lik1 + lik2 + mBeta/2 * th(0) -
    (mAlpha+1) * log(th(0));
}


ProbabilityDistribution* GaussianProcessIGN::prediction( const vectord &query )
{
  double kn;
  double uInvRr, rInvRr, rInvRy;
  double meanf = mMean->getMean(query);

  vectord colR = computeCrossCorrelation(query);
  kn = (*mKernel)(query, query);

#if USE_CHOL
  vectord invRr(colR);
  inplace_solve(mL,invRr,lower_tag());
  rInvRr = inner_prod(invRr,invRr);
  uInvRr = inner_prod(mUInvR, invRr);  
  rInvRy = inner_prod(invRr, mInvRy );
#else
  vectord rInvR = prod(colR,mInvR);
	
  rInvRr = inner_prod(rInvR,colR);
  uInvRr = inner_prod(mUInvR,colR);
  rInvRy = inner_prod(colR,mInvRy);
#endif
  
  double yPred = meanf*mMu + rInvRy;
  double sPred = sqrt( mSig * (kn - rInvRr + (meanf - uInvRr) * (meanf - uInvRr) 
			       / mUInvRUDelta ) );

  d_->setMeanAndStd(yPred,sPred);
  return d_.get();
}

int GaussianProcessIGN::precomputePrediction()
{
  size_t ns = mGPXX.size();
  size_t nf = mMean->nFeatures();

  mFKF =

  mUInvR.resize(nSamples,false);
  mInvRy.resize(nSamples,false);

#if USE_CHOL
  mAlphaV.resize(nSamples,false);
  mAlphaV = mGPY;
  utils::cholesky_solve(mL,mAlphaV,lower());

  //  vectord alphaMean(mMeanV);
  mUInvR = mMeanV;
  inplace_solve(mL,mUInvR,lower_tag());
  mUInvRUDelta = inner_prod(mUInvR,mUInvR) + 1/mDelta2;

  mMu = inner_prod(mMeanV,mAlphaV) / mUInvRUDelta;
  double YInvRY = inner_prod(mGPY,mAlphaV);

  mInvRy = mGPY - mMeanV*mMu;
  inplace_solve(mL,mInvRy,lower_tag());
#else
  // mInvR.resize(nSamples,nSamples);
  // mInvR.assign(boost::numeric::ublas::identity_matrix<double>(nSamples));
  // cholesky_solve(mL,mInvR,lower());
  mUInvR = prod(mMeanV,mInvR);
  mUInvRUDelta = inner_prod(mUInvR,mMeanV) + 1/mDelta2;
  
  vectord YInvR(nSamples);
  double YInvRY;
  
  mMu =  inner_prod(mUInvR,mGPY) / mUInvRUDelta;
  
  noalias(YInvR) = prod(mGPY,mInvR);
  YInvRY = inner_prod(YInvR,mGPY);
  vectord ymu = mGPY - mMeanV*mMu;
  mInvRy = prod(mInvR,ymu);
#endif  

  mSig = (mBeta + YInvRY - mMu*mMu*mUInvRUDelta) / (mAlpha + (nSamples+1) + 2);
  
  return 1;
}
