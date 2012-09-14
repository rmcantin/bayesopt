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

#include "nonparametricprocess.hpp"
#include "cholesky.hpp"

NonParametricProcess::NonParametricProcess(double noise):
  InnerOptimization(),  mKernel(NULL), mRegularizer(noise)
{ 
  mMinIndex = 0; 
  mMaxIndex = 0;   
  setAlgorithm(BOBYQA);
  setLimits(0.,100.);
};


NonParametricProcess::~NonParametricProcess()
{
  if (mKernel != NULL)
    delete mKernel;
};

int NonParametricProcess::setKernel (const vectord &thetav,
				     kernel_name k_name)
{
  if (mKernel != NULL)
    delete mKernel;

  switch(k_name)
    {
    case K_MATERN_ISO1: mKernel = new MaternIso1(); break;
    case K_MATERN_ISO3: mKernel = new MaternIso3(); break;
    case K_MATERN_ISO5: mKernel = new MaternIso5(); break;
    case K_SE_ISO: mKernel = new SEIso(); break;
    case K_SE_ARD: mKernel = new SEArd(); break;
    default:
      std::cout << "Error: kernel function not supported." << std::endl;
      return -1;
    }

  mKernel->setScale(thetav);

  return 0;
};

int NonParametricProcess::fitGP()
{
  size_t nSamples = mGPXX.size();
  for (size_t ii=0; ii<nSamples; ++ii)
    checkBoundsY(ii);

  vectord optimalTheta = mKernel->getScale();
  
  std::cout << "Initial theta: " << optimalTheta << std::endl;
  innerOptimize(optimalTheta);
  mKernel->setScale(optimalTheta);
  std::cout << "Final theta: " << optimalTheta << std::endl;

  int error = computeInverseCorrMatrix();

  if (error < 0)
    return error;

  return precomputeGPParams();
} // fitGP


int NonParametricProcess::addNewPointToGP( const vectord &Xnew,
					   double Ynew)
{
  size_t nSamples = mGPXX.size();
  size_t XDim = mGPXX[1].size();
  
  vectord Li(nSamples);
  vectord wInvR(nSamples);
  double wInvRw;
  double selfCorrelation, Ni;
  
  if (XDim != Xnew.size())
    {
      std::cout << "Dimensional Error" << std::endl;
      return -1;
    }
    
  vectord correlationNewValue = computeCrossCorrelation(Xnew);
  
  selfCorrelation = (*mKernel)(Xnew, Xnew) + mRegularizer;
  
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



int NonParametricProcess::computeInverseCorrMatrix()
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

  /*  std::transform(mGPXX.begin(),mGPXX.end(),knx.begin(),
      std::bind2nd(mKernel(),query))*/
}
