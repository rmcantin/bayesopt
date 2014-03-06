
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
#include <cstdlib>
#include "nonparametricprocess.hpp"
#include "log.hpp"
#include "cholesky.hpp"
#include "ublas_extra.hpp"

#include "gaussian_process.hpp"
#include "gaussian_process_ml.hpp"
#include "gaussian_process_normal.hpp"
#include "student_t_process_jef.hpp"
#include "student_t_process_nig.hpp"


namespace bayesopt
{
  
  NonParametricProcess::NonParametricProcess(size_t dim, bopt_params parameters):
    dim_(dim), mRegularizer(parameters.noise),
    mKernel(dim, parameters), mMean(dim, parameters)
  { 
    setLearnType(parameters.l_type);
  }

  NonParametricProcess::~NonParametricProcess()
  {}


  NonParametricProcess* NonParametricProcess::create(size_t dim, 
						     bopt_params parameters)
  {
    NonParametricProcess* s_ptr;

    std::string name = parameters.surr_name;

    if (!name.compare("sGaussianProcess"))
      s_ptr = new GaussianProcess(dim,parameters);
    else  if(!name.compare("sGaussianProcessML"))
      s_ptr = new GaussianProcessML(dim,parameters);
    else  if(!name.compare("sGaussianProcessNormal"))
      s_ptr = new GaussianProcessNormal(dim,parameters);
    else if (!name.compare("sStudentTProcessJef"))
      s_ptr = new StudentTProcessNIG(dim,parameters); 
    else if (!name.compare("sStudentTProcessNIG"))
      s_ptr = new StudentTProcessNIG(dim,parameters); 
    else
      {
	FILE_LOG(logERROR) << "Error: surrogate function not supported.";
	return NULL;
      }
    return s_ptr;
  };


  int NonParametricProcess::fitSurrogateModel()
  {
    int error = updateKernelParameters();
    error += precomputeSurrogate();
    return error;
  };


  int NonParametricProcess::precomputeSurrogate()
  {
    int error = -1;
    error = computeCholeskyCorrelation();
    if (error)
      {
	FILE_LOG(logERROR) << "Error computing the correlation matrix";
	exit(EXIT_FAILURE);
      }   

    error = precomputePrediction(); 
    if (error)
      {
	FILE_LOG(logERROR) << "Error pre-computing the prediction distribution";
	exit(EXIT_FAILURE);
      }   

    return error; 
  } // fitSurrogateModel


  int NonParametricProcess::updateSurrogateModel( const vectord &Xnew,
						  double Ynew, bool retrain)
  {
    assert( Dataset.mX[1].size() == Xnew.size() );

    if (retrain)
      {
	addSample(Xnew,Ynew);
	FILE_LOG(logDEBUG) << "Retraining model parameters";
	return fitSurrogateModel();	
      }
    else
      {
	const vectord newK = computeCrossCorrelation(Xnew);
	double selfCorrelation = computeSelfCorrelation(Xnew) + mRegularizer;
	addSample(Xnew,Ynew);
	addNewPointToCholesky(newK,selfCorrelation);

	int error = precomputePrediction(); 
	if (error < 0)
	  {
	    FILE_LOG(logERROR) << "Error pre-computing the prediction distribution";
	    exit(EXIT_FAILURE);
	  }   
	return error; 
      }
    return 0; //JIC
  } // updateSurrogateModel


  //////////////////////////////////////////////////////////////////////////////
  //// Getters and Setters
  void NonParametricProcess::setSamples(const matrixd &x, const vectord &y)
  {
    mData.setSamples(x,y);
    mMean.setPoints(mData.mX);  //Because it expects a vecOfvec instead of a matrixd
  }

  void NonParametricProcess::addSample(const vectord &x, double y)
  {
    mData.addSample(x,y);
    mMean.addNewPoint(x);
  };

  
  int NonParametricProcess::addNewPointToCholesky(const vectord& correlation,
						  double selfcorrelation)
  {
    vectord newK(correlation);
    utils::append(newK, selfcorrelation);
    utils::cholesky_add_row(mL,newK);
    return 0;
  }


  int NonParametricProcess::computeCholeskyCorrelation()
  {
    size_t nSamples = mData.getNSamples();
    mL.resize(nSamples,nSamples);
  
    //  const matrixd K = computeCorrMatrix();
    matrixd K(nSamples,nSamples);
    computeCorrMatrix(K);
    return utils::cholesky_decompose(K,mL);
  }


  int NonParametricProcess::computeCorrMatrix(matrixd& corrMatrix)
  {
    return mKernel.computeCorrMatrix(mData.mX,corrMatrix,mRegularizer);
  }



  matrixd NonParametricProcess::computeCorrMatrix()
  {
    const size_t nSamples = mData.getNSamples();
    matrixd corrMatrix(nSamples,nSamples);
    int error = mKernel.computeCorrMatrix(mData.mX,corrMatrix,mRegularizer);
    return corrMatrix;
  }

  matrixd NonParametricProcess::computeDerivativeCorrMatrix(int dth_index)
  {
    const size_t nSamples = mData.getNSamples();
    matrixd corrMatrix(nSamples,nSamples);
    int error = mKernel.computeDerivativeCorrMatrix(mData.mX,corrMatrix,dth_index);
    return corrMatrix;
  }


  vectord NonParametricProcess::computeCrossCorrelation(const vectord &query)
  {
    return mKernel.computeCrossCorrelation(mData.mX,query);
  }

  double NonParametricProcess::computeSelfCorrelation(const vectord& query)
  {
    return mKernel.computeSelfCorrelation(query);
  }

} //namespace bayesopt
