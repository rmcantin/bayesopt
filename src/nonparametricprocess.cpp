/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2013 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
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


namespace bayesopt
{
  
  NonParametricProcess::NonParametricProcess(size_t dim, double noise):
    InnerOptimization(), mRegularizer(noise), dim_(dim)
  { 
    mMinIndex = 0;     mMaxIndex = 0;   
    setAlgorithm(BOBYQA);
    setLimits(0.,100.);
  }

  NonParametricProcess::~NonParametricProcess(){}


  NonParametricProcess* NonParametricProcess::create(size_t dim, 
						     bopt_params parameters)
  {
    NonParametricProcess* s_ptr;

    switch(parameters.s_name)
      {
      case S_GAUSSIAN_PROCESS: 
	s_ptr = new GaussianProcess(dim,parameters.noise); 
	break;

      case S_GAUSSIAN_PROCESS_INV_GAMMA_NORMAL:
	s_ptr = new GaussianProcessIGN(dim,parameters.noise, parameters.alpha,
				       parameters.beta,parameters.delta);  
	break;

      case S_STUDENT_T_PROCESS_JEFFREYS: 
	if (parameters.m_name == M_ZERO)
	  {
	    FILE_LOG(logWARNING) << "Zero mean incompatible with Student's t "
				 << "process, using one-mean instead.";
	    parameters.m_name = M_ONE;
	  }
	s_ptr = new StudentTProcess(dim,parameters.noise); 
	break;

      default:
	FILE_LOG(logERROR) << "Error: surrogate function not supported.";
	return NULL;
      }
  
    s_ptr->setKernel(parameters.theta,parameters.n_theta,
		     parameters.k_s_name,dim);
    s_ptr->setKernelPrior(parameters.theta,
			  parameters.s_theta,parameters.n_theta);
    s_ptr->setMean(parameters.mu,parameters.n_mu,parameters.m_s_name,dim);
    return s_ptr;
  };


  int NonParametricProcess::fitInitialSurrogate()
  {
    vectord optimalTheta = mKernel->getHyperParameters();
    int error = -1;

    FILE_LOG(logDEBUG) << "Computing kernel parameters. Seed: " << optimalTheta;
    innerOptimize(optimalTheta);
    mKernel->setHyperParameters(optimalTheta);
    FILE_LOG(logDEBUG) << "Final kernel parameters: " << optimalTheta;

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


  int NonParametricProcess::fullUpdateSurrogateModel( const vectord &Xnew,
						      double Ynew)
  {
    assert( mGPXX[1].size() == Xnew.size() );
    addSample(Xnew,Ynew);
    return fitInitialSurrogate();
  } // fullUpdateSurrogateModel


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
    mFeatM = mMean->getAllFeatures(mGPXX);
  }

  void NonParametricProcess::addSample(const vectord &x, double y)
  {
    using boost::numeric::ublas::column;

    mGPXX.push_back(x);
    mGPY.resize(mGPY.size()+1);  mGPY(mGPY.size()-1) = y;
    checkBoundsY(mGPY.size()-1);

    mMeanV.resize(mMeanV.size()+1);  
    mMeanV(mMeanV.size()-1) = mMean->getMean(x);

    mFeatM.resize(mFeatM.size1(),mFeatM.size2()+1);  
    column(mFeatM,mFeatM.size2()-1) = mMean->getFeatures(x);

  };

  double NonParametricProcess::getSample(size_t index, vectord &x)
  {
    x = mGPXX[index];
    return mGPY(index);
  }

  int NonParametricProcess::setKernel (const vectord &thetav, 
				       kernel_name k_name, 
				       size_t dim)
  {
    mKernel.reset(mKFactory.create(k_name, dim));
    if (mKernel == NULL)   
      {
	return -1;
      }
    else
      {
	mKernel->setHyperParameters(thetav);
	return 0;
      }
  }


  int NonParametricProcess::setKernel (const vectord &thetav, 
				       std::string k_name, 
				       size_t dim)
  {
    mKernel.reset(mKFactory.create(k_name, dim));
    if (mKernel == NULL)   
      {
	return -1;
      }
    else
      {
	mKernel->setHyperParameters(thetav);
	return 0;
      }
  }


  int NonParametricProcess::setMean (const vectord &muv,
				     mean_name m_name,
				     size_t dim)
  {
    mMean.reset(mPFactory.create(m_name,dim));
    if (mMean == NULL) 
      {
	return -1; 
      }
    else 
      {
	mMean->setParameters(muv);
	return 0;
      } 
  }

  int NonParametricProcess::setMean (const vectord &muv,
				     std::string m_name,
				     size_t dim)
  {
    mMean.reset(mPFactory.create(m_name,dim));
    if (mMean == NULL) 
      {
	return -1; 
      }
    else 
      {
	mMean->setParameters(muv);
	return 0;
      } 
  }



  int NonParametricProcess::addNewPointToCholesky(const vectord& correlation,
						  double selfcorrelation)
  {
    vectord newK(correlation);
    append(newK, selfcorrelation);
    utils::cholesky_add_row(mL,newK);
    return 1;
  }


  int NonParametricProcess::computeCholeskyCorrelation()
  {
    size_t nSamples = mGPXX.size();
    mL.resize(nSamples,nSamples);
  
    //  const matrixd K = computeCorrMatrix();
    matrixd K(nSamples,nSamples);
    computeCorrMatrix(K);
    return utils::cholesky_decompose(K,mL);
  }

  int NonParametricProcess::addNewPointToInverse(const vectord& correlation,
						 double selfcorrelation)
  {
    size_t nSamples = correlation.size();
  
    vectord wInvR = prod(correlation,mInvR);
    double wInvRw = inner_prod(wInvR,correlation);
    double Ni = 1/(selfcorrelation - wInvRw);
    vectord Li = -Ni * wInvR;
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
    return utils::inverse_cholesky(corrMatrix,mInvR);
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
	    corrMatrix(ii,jj) = mKernel->gradient(mGPXX[ii],mGPXX[jj], 
						  dth_index);
	    corrMatrix(jj,ii) = corrMatrix(ii,jj);
	  }
	corrMatrix(ii,ii) = mKernel->gradient(mGPXX[ii],mGPXX[ii],dth_index);
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
    
    return knx;
  }


} //namespace bayesopt
