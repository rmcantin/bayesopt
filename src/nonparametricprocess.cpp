/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2012 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/


#include <cstdio>
#include "nonparametricprocess.hpp"
#include "cholesky.hpp"
#include "ublas_extra.hpp"

NonParametricProcess::NonParametricProcess(double noise):
  InnerOptimization(), mRegularizer(noise)
{ 
  kernelSet = false; meanSet = false;
  mMinIndex = 0;     mMaxIndex = 0;   

  setAlgorithm(BOBYQA);
  setLimits(0.,100.);
}


NonParametricProcess::~NonParametricProcess()
{
  if (kernelSet) delete mKernel;
  if (meanSet) delete mMean;
}

int NonParametricProcess::setKernel (const vectord &thetav,
				     kernel_name k_name)
{
  if (kernelSet)   delete mKernel;

  kernelSet = false;

  mKernel = Kernel::create(k_name, thetav);
  if (mKernel == NULL)   return -1;

  kernelSet = true;
  return 0;
}


int NonParametricProcess::setMean (const vectord &muv,
				   mean_name m_name)
{
  if (meanSet)  delete mMean;
  meanSet = false;

  mMean = ParametricFunction::create(m_name,muv);
  if (mMean == NULL) return -1;

  meanSet = true;
  return 0;
}


int NonParametricProcess::fitGP()
{
  vectord optimalTheta = mKernel->getScale();
  int error = -1;

  std::cout << "Initial theta: " << optimalTheta << std::endl;
  innerOptimize(optimalTheta);
  mKernel->setScale(optimalTheta);
  std::cout << "Final theta: " << optimalTheta << std::endl;

  //#ifdef USE_CHOLESKY
  error = computeCholeskyCorrelation();
  //#else
  error = computeInverseCorrelation();
  //#end

  if (error < 0)   
    return error;

  return precomputePrediction();
} // fitGP

int NonParametricProcess::addNewPointToGP( const vectord &Xnew,
					   double Ynew)
{
  assert( mGPXX[1].size() == Xnew.size() );

  vectord newK = computeCrossCorrelation(Xnew);
  double selfCorrelation = (*mKernel)(Xnew, Xnew) + mRegularizer;
  
  addSample(Xnew,Ynew);

  //TODO: Choose one!
#if USE_CHOL
  addNewPointToCholesky(newK,selfCorrelation);
#else
  addNewPointToInverse(newK,selfCorrelation);
#endif
  
  return precomputePrediction();
} // addNewPointToGP


int NonParametricProcess::addNewPointToCholesky(const vectord& correlation,
						double selfcorrelation)
{
  //TODO: Optimize
  vectord newK = correlation;
  append(newK, selfcorrelation);
  cholesky_add_row(mL,newK);

  return 1;
}


int NonParametricProcess::computeCholeskyCorrelation()
{
  size_t nSamples = mGPXX.size();
  mL.resize(nSamples,nSamples);
  
  const matrixd K = computeCorrMatrix();
  return cholesky_decompose(K,mL);
}

int NonParametricProcess::addNewPointToInverse(const vectord& correlation,
					       double selfcorrelation)
{
  size_t nSamples = correlation.size();
  
  vectord Li(nSamples);
  vectord wInvR(nSamples);
  double wInvRw;
  double Ni;

  noalias(wInvR) = prod(correlation,mInvR);
  wInvRw = inner_prod(wInvR,correlation);

  Ni = 1/(selfcorrelation - wInvRw);

  noalias(Li) = -Ni * wInvR;
  mInvR += outer_prod(Li,Li) / Ni;
  
  //TODO: There must be a better way to do this.
  mInvR.resize(nSamples+1,nSamples+1,true);
  
  Li.resize(nSamples+1);
  Li(nSamples) = Ni;
  
  row(mInvR,nSamples) = Li;
  column(mInvR,nSamples) = Li;

  return 1;

}


int NonParametricProcess::computeInverseCorrelation()
{
  size_t nSamples = mGPXX.size();
  if ( (nSamples != mInvR.size1()) || (nSamples != mInvR.size2()) )
    mInvR.resize(nSamples,nSamples);
    
  matrixd corrMatrix = computeCorrMatrix();

  //return InvertMatrix(corrMatrix,mInvR);
  return inverse_cholesky(corrMatrix,mInvR);
}




matrixd NonParametricProcess::computeCorrMatrix(int dth_index)
{
  size_t nSamples = mGPXX.size();
  matrixd corrMatrix(nSamples,nSamples);
  
  for (size_t ii=0; ii< nSamples; ++ii)
    {
      for (size_t jj=0; jj < ii; ++jj)
	{
	  corrMatrix(ii,jj) = (*mKernel)(mGPXX[ii], mGPXX[jj], 
						  dth_index);
	  corrMatrix(jj,ii) = corrMatrix(ii,jj);
	}
      corrMatrix(ii,ii) = (*mKernel)(mGPXX[ii],mGPXX[ii], dth_index);
      if (dth_index < 0) 
	corrMatrix(ii,ii) += mRegularizer;
    }
  return corrMatrix;
}



vectord NonParametricProcess::computeCrossCorrelation(const vectord &query)
{
  vectord knx(mGPXX.size());

  //TODO: Replace by transform
  for (size_t ii=0; ii<mGPXX.size(); ++ii)
    knx(ii) = (*mKernel)(mGPXX[ii], query);

  return knx;
}
