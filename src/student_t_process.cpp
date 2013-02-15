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

#include "student_t_process.hpp"
#include "cholesky.hpp"
#include "trace_ublas.hpp"
#include "student_t_distribution.hpp"

using boost::numeric::ublas::inplace_solve;
using boost::numeric::ublas::lower_tag;
using boost::numeric::ublas::lower;

  
StudentTProcess::StudentTProcess(size_t dim, double noise):
  NonParametricProcess(dim, noise)
{}  // Constructor



StudentTProcess::~StudentTProcess()
{} // Default destructor


double StudentTProcess::negativeLogLikelihood()
{
  matrixd K = computeCorrMatrix();
  size_t n = K.size1();
  matrixd L(n,n);
  cholesky_decompose(K,L);

  vectord alphU(mMeanV);
  inplace_solve(L,alphU,lower_tag());
  double eta = inner_prod(alphU,alphU);
  
  vectord alphY(mGPY);
  inplace_solve(L,alphY,lower_tag());
  double mu     = inner_prod(alphU,alphY) / eta;
  double YInvRY = inner_prod(alphY,alphY);
    
  double sigma = (YInvRY - mu*mu*eta) / (n-1);

  // TODO: Revise equation.
  // Obtained from SEQUENTIAL DESIGN OF COMPUTER EXPERIMENTS
  // TO MINIMIZE INTEGRATED RESPONSE FUNCTIONS
  // http://www3.stat.sinica.edu.tw/statistica/oldpdf/A10n46.pdf
  double negloglik = log_trace(L) + 0.5*( (n-1)*log(sigma) + log(eta) );

  return negloglik;
}


ProbabilityDistribution* StudentTProcess::prediction(const vectord &query)
{
  double yPred, sPred, rInvRu, rInvRr, rInvRy;
  double meanf = mMean->getMean(query);
  
  vectord colR = computeCrossCorrelation(query);
  double kn = (*mKernel)(query, query);
  
#if USE_CHOL
  vectord invRr(colR);
  inplace_solve(mL,invRr,lower_tag());
  rInvRr = inner_prod(invRr,invRr);
  rInvRu = inner_prod(invRr,mUInvR);  
  rInvRy = inner_prod(invRr, mInvRy );
#else
  vectord rInvR(mGPXX.size());
  noalias(rInvR) = prod(colR,mInvR);	
  rInvRr = inner_prod(rInvR,colR);
  rInvRu = inner_prod(rInvR,mMeanV);
  rInvRy = inner_prod(colR,mInvRy);
#endif

  yPred = meanf*mMu + rInvRy;
  sPred = sqrt( mSig * (kn - rInvRr + (meanf - rInvRu) * (meanf - rInvRu) 
			/ mUInvRUDelta ) );

  d_->setMeanAndStd(yPred,sPred);
  return d_.get();
}


int StudentTProcess::precomputePrediction()
{
  size_t nSamples = mGPXX.size();
  double YInvRY;
  d_.reset(new StudentTDistribution(nSamples-1));

  mInvRy.resize(nSamples,false);

#if USE_CHOL
  mUInvR = mMeanV;
  inplace_solve(mL,mUInvR,lower_tag());
  mUInvRUDelta = inner_prod(mUInvR,mUInvR);

  mInvRy = mGPY;
  inplace_solve(mL,mInvRy,lower_tag());
  mMu =  inner_prod(mUInvR,mInvRy) / mUInvRUDelta;
  YInvRY = inner_prod(mInvRy,mInvRy);

  mInvRy = mGPY - mMeanV*mMu;
  inplace_solve(mL,mInvRy,lower_tag()); 
#else
  mUInvR = prod(mMeanV,mInvR);
  mUInvRUDelta = inner_prod(mUInvR,mMeanV);
  mMu =  inner_prod(mUInvR,mGPY) / mUInvRUDelta;
  noalias(mInvRy) = prod(mInvR,mGPY); 
  YInvRY = inner_prod(mGPY,mInvRy);

  vectord yumu = mGPY - mMeanV*mMu;
  mInvRy = prod(mInvR,yumu);
#endif
  
  mSig = (YInvRY - mMu*mMu*mUInvRUDelta) / (nSamples-1);

  return 1;
}


	
