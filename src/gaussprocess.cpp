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


#include "gaussprocess.hpp"

  

GaussianProcess::GaussianProcess():
  NonParametricProcess(), mTheta(KERNEL_THETA), mP(KERNEL_P),
  mAlpha(PRIOR_ALPHA), mBeta (PRIOR_BETA),
  mDelta2(PRIOR_DELTA_SQ), mRegularizer(DEF_REGULARIZER)
{} // Default constructor

GaussianProcess::GaussianProcess( double theta, double p,
				  double alpha, double beta, 
				  double delta, double noise):
  NonParametricProcess(), mTheta(theta), mP(p),
  mAlpha(alpha), mBeta (beta),
  mDelta2(delta), mRegularizer(noise)
{}  // Constructor

GaussianProcess::GaussianProcess( gp_params params ):
  NonParametricProcess(), mTheta(params.theta), mP(params.p),
  mAlpha(params.alpha), mBeta(params.beta),
  mDelta2(params.delta), mRegularizer(params.noise)
{}  // Constructor



GaussianProcess::~GaussianProcess()
{} // Default destructor



/*
double GaussianProcess::negativeLogLikelihood(double param,
					      vectord &grad)
{
  vectord alpha = dot(mInvR,mGPY);
    
}
*/

int GaussianProcess::prediction( const vectord &query,
				 double& yPred, double& sPred)
{
  vectord rInvR(mGPXX.size());
  double kn;
  double uInvRr, rInvRr;

  vectord colR = computeCrossCorrelation(query);
  kn = correlationFunction(query, query);
  
  noalias(rInvR) = prod(colR,mInvR);	
  rInvRr = inner_prod(rInvR,colR);
  uInvRr = inner_prod(mUInvR,colR);
  
  yPred = mMu + inner_prod( rInvR, mYUmu );
  sPred = sqrt( mSig * (kn - rInvRr + (1.0 - uInvRr) * (1.0 - uInvRr) 
			/ mUInvRUDelta ) );

  return 1;
}
	

int GaussianProcess::fitGP()
{
  size_t nSamples = mGPXX.size();
  for (size_t ii=0; ii<nSamples; ii++)
    checkBoundsY(ii);

  int error = computeInverseCorrMatrix(mRegularizer);

  if (error < 0)
    return error;

  return precomputeGPParams();
} // fitGP


int GaussianProcess::addNewPointToGP(const vectord &Xnew, 
				     double Ynew)
{
  size_t nSamples = mGPXX.size();
  size_t XDim = mGPXX[1].size();
  size_t NewDim = Xnew.size();
  
  svectord colU(nSamples+1,1.0);
  vectord Li(nSamples);
  vectord wInvR(nSamples);
  double wInvRw;
  double selfCorrelation, Ni;
  
  if (XDim != NewDim)
    {
      std::cout << "Dimensional Error" << std::endl;
      return -1;
    }
    
  vectord correlationNewValue = computeCrossCorrelation(Xnew);
  
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

  addSample(Xnew,Ynew);
  checkBoundsY(nSamples);
  
  return precomputeGPParams();
} // addNewPointToGP


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
  
  mSig = (mBeta + YInvRY - mMu*mMu/mUInvRUDelta) / (mAlpha + (nSamples+1) + 2);
  
  svectord colMu(nSamples,mMu);
  mYUmu = mGPY - colMu;
  
  return 1;
}


