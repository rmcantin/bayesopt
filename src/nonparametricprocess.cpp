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
#include "log.hpp"
#include "cholesky.hpp"
#include "ublas_extra.hpp"

#include "gaussian_process.hpp"
#include "gaussian_process_ign.hpp"
#include "student_t_process.hpp"


NonParametricProcess::NonParametricProcess(double noise):
  InnerOptimization(), mRegularizer(noise)
{ 
  mMinIndex = 0;     mMaxIndex = 0;   
  setAlgorithm(BOBYQA);
  setLimits(0.,100.);
}

NonParametricProcess::~NonParametricProcess(){}


NonParametricProcess* NonParametricProcess::create(bopt_params parameters)
{
  NonParametricProcess* s_ptr;

  switch(parameters.s_name)
    {
    case S_GAUSSIAN_PROCESS: 
      s_ptr = new GaussianProcess(parameters.noise); 
      break;

    case S_GAUSSIAN_PROCESS_INV_GAMMA_NORMAL:
      s_ptr = new GaussianProcessIGN(parameters.noise, parameters.alpha,
				     parameters.beta,parameters.delta);  
      break;

    case S_STUDENT_T_PROCESS_JEFFREYS: 
      if (parameters.m_name == M_ZERO)
	{
	  FILE_LOG(logWARNING) << "Zero mean incompatible with Student's t process,"
			       << "using one-mean instead.";
	  parameters.m_name = M_ONE;
	}
      s_ptr = new StudentTProcess(parameters.noise); 
      break;

    default:
      FILE_LOG(logERROR) << "Error: surrogate function not supported.";
      return NULL;
    }
  
  s_ptr->setKernel(parameters.theta,parameters.n_theta,parameters.k_name);
  s_ptr->setMean(parameters.mu,parameters.n_mu,parameters.m_name);
  return s_ptr;
};


int NonParametricProcess::fitInitialSurrogate()
{
  vectord optimalTheta = mKernel->getScale();
  int error = -1;

  FILE_LOG(logDEBUG) << "Initial theta: " << optimalTheta;
  innerOptimize(optimalTheta);
  mKernel->setScale(optimalTheta);
  FILE_LOG(logDEBUG) << "Final theta: " << optimalTheta;

  //TODO: Choose one!
#if USE_CHOL
  error = computeCholeskyCorrelation();
#else
  error = computeInverseCorrelation();
#endif

  if (error < 0)
    {
      FILE_LOG(logERROR) << "Error computing the correlation matrix";
      return error;
    }   

  error = precomputePrediction(); 

  if (error < 0)
    {
      FILE_LOG(logERROR) << "Error pre-computing the prediction distribution";
      return error;
    }   

  return 0; 
} // fitInitialSurrogate


int NonParametricProcess::updateSurrogateModel( const vectord &Xnew,
						double Ynew)
{
  assert( mGPXX[1].size() == Xnew.size() );

  const vectord newK = computeCrossCorrelation(Xnew);
  double selfCorrelation = (*mKernel)(Xnew, Xnew) + mRegularizer;
  
  addSample(Xnew,Ynew);

  //TODO: Choose one!
#if USE_CHOL
  addNewPointToCholesky(newK,selfCorrelation);
#else
  addNewPointToInverse(newK,selfCorrelation);
#endif

  int error = precomputePrediction(); 
  if (error < 0)
    {
      FILE_LOG(logERROR) << "Error pre-computing the prediction distribution";
      return error;
    }   

  return 0; 
} // updateSurrogateModel


//////////////////////////////////////////////////////////////////////////////
//// Getters and Setters
void NonParametricProcess::setSamples(const matrixd &x, const vectord &y)
{
  mGPY = y;
  for (size_t i=0; i<x.size1(); ++i)
    {
      mGPXX.push_back(row(x,i));
      checkBoundsY(i);
    } 
  mMeanV = (*mMean)(mGPXX);
}

void NonParametricProcess::addSample(const vectord &x, double y)
{
  mGPXX.push_back(x);
  mGPY.resize(mGPY.size()+1);  mGPY(mGPY.size()-1) = y;
  checkBoundsY(mGPY.size()-1);
  mMeanV.resize(mMeanV.size()+1);  
  mMeanV(mMeanV.size()-1) = mMean->getMean(x);
};

double NonParametricProcess::getSample(size_t index, vectord &x)
{
  x = mGPXX[index];
  return mGPY(index);
}

vectord NonParametricProcess::getPointAtMinimum()
{ return mGPXX[mMinIndex]; };

double NonParametricProcess::getValueAtMinimum()
{ return mGPY(mMinIndex); };

int NonParametricProcess::setKernel (const vectord &thetav,
				     kernel_name k_name)
{
  mKernel.reset(Kernel::create(k_name, thetav));
  if (mKernel == NULL)   return -1;
  else  return 0;
}


int NonParametricProcess::setMean (const vectord &muv,
				   mean_name m_name)
{
  mMean.reset(ParametricFunction::create(m_name,muv));
  if (mMean == NULL) return -1;
  else  return 0;
}


int NonParametricProcess::addNewPointToCholesky(const vectord& correlation,
						double selfcorrelation)
{
  vectord newK(correlation);
  append(newK, selfcorrelation);
  cholesky_add_row(mL,newK);
  return 1;
}


int NonParametricProcess::computeCholeskyCorrelation()
{
  size_t nSamples = mGPXX.size();
  mL.resize(nSamples,nSamples);
  
  //  const matrixd K = computeCorrMatrix();
  matrixd K(nSamples,nSamples);
  computeCorrMatrix(K);
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
  const size_t nSamples = mGPXX.size();
  if ( (nSamples != mInvR.size1()) || (nSamples != mInvR.size2()) )
    mInvR.resize(nSamples,nSamples);
    
  const matrixd corrMatrix = computeCorrMatrix();
  return inverse_cholesky(corrMatrix,mInvR);
}


int NonParametricProcess::computeCorrMatrix(matrixd& corrMatrix)
{
  assert(corrMatrix.size1() == mGPXX.size());
  assert(corrMatrix.size2() == mGPXX.size());
  const size_t nSamples = mGPXX.size();
  
  for (size_t ii=0; ii< nSamples; ++ii)
    {
      for (size_t jj=0; jj < ii; ++jj)
	{
	  corrMatrix(ii,jj) = (*mKernel)(mGPXX[ii], mGPXX[jj]);
	  corrMatrix(jj,ii) = corrMatrix(ii,jj);
	}
      corrMatrix(ii,ii) = (*mKernel)(mGPXX[ii],mGPXX[ii]) + mRegularizer;
    }
  return 1;
}



matrixd NonParametricProcess::computeCorrMatrix()
{
  const size_t nSamples = mGPXX.size();
  matrixd corrMatrix(nSamples,nSamples);
  
  for (size_t ii=0; ii< nSamples; ++ii)
    {
      for (size_t jj=0; jj < ii; ++jj)
	{
	  corrMatrix(ii,jj) = (*mKernel)(mGPXX[ii], mGPXX[jj]);
	  corrMatrix(jj,ii) = corrMatrix(ii,jj);
	}
      corrMatrix(ii,ii) = (*mKernel)(mGPXX[ii],mGPXX[ii]) + mRegularizer;
    }
  return corrMatrix;
}

matrixd NonParametricProcess::computeDerivativeCorrMatrix(int dth_index)
{
  const size_t nSamples = mGPXX.size();
  matrixd corrMatrix(nSamples,nSamples);
  
  for (size_t ii=0; ii< nSamples; ++ii)
    {
      for (size_t jj=0; jj < ii; ++jj)
	{
	  corrMatrix(ii,jj) = mKernel->getGradient(mGPXX[ii],mGPXX[jj], 
						   dth_index);
	  corrMatrix(jj,ii) = corrMatrix(ii,jj);
	}
      corrMatrix(ii,ii) = mKernel->getGradient(mGPXX[ii],mGPXX[ii],dth_index);
    }
  return corrMatrix;
}



vectord NonParametricProcess::computeCrossCorrelation(const vectord &query)
{
  vectord knx(mGPXX.size());

  std::vector<vectord>::const_iterator x_it  = mGPXX.begin();
  std::vector<vectord>::const_iterator x_end = mGPXX.end();
  vectord::iterator k_it = knx.begin();
  while(x_it != x_end)
    {
      *k_it++ = (*mKernel)(*x_it++, query);
    }
    
  //TODO: Replace by transform
  // for (size_t ii=0; ii<mGPXX.size(); ++ii)
  //   knx(ii) = (*mKernel)(mGPXX[ii], query);

  return knx;
}

