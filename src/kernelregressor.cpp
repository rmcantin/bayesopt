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
#include <stdexcept>
#include <boost/lexical_cast.hpp>

#include "kernelregressor.hpp"

#include "log.hpp"
#include "ublas_extra.hpp"


namespace bayesopt
{
  KernelRegressor::KernelRegressor(size_t dim, bopt_params parameters,
				   Dataset& data):
    NonParametricProcess(dim,parameters,data), mRegularizer(parameters.noise),
    mKernel(dim, parameters)
  { 
    setLearnType(parameters.l_type);
  }

  KernelRegressor::~KernelRegressor(){}



  void KernelRegressor::updateSurrogateModel( const vectord &Xnew,
					      double Ynew, bool retrain)
  {
    assert( mData.mX[1].size() == Xnew.size() );

    if (retrain)
      {
	FILE_LOG(logDEBUG) << "Retraining model parameters";
	addSample(Xnew,Ynew);
	fitSurrogateModel();	
      }
    else
      {
	addSample(Xnew,Ynew);
	vectord newK = computeCrossCorrelation(Xnew);
	newK(newK.size()) += mRegularizer;
	utils::cholesky_add_row(mL,newK);
	//double selfCorrelation = computeSelfCorrelation(Xnew) + mRegularizer;
	//addNewPointToCholesky(newK,selfCorrelation);
	precomputePrediction(); 
      }
  } // updateSurrogateModel


  void KernelRegressor::computeCholeskyCorrelation()
  {
    size_t nSamples = mData.getNSamples();
    mL.resize(nSamples,nSamples);
  
    //  const matrixd K = computeCorrMatrix();
    matrixd K(nSamples,nSamples);
    computeCorrMatrix(K);
    size_t line_error = utils::cholesky_decompose(K,mL);
    if (line_error) 
      throw std::runtime_error("Cholesky decomposition error at line " + 
			       boost::lexical_cast<std::string>(line_error));
  }

  matrixd KernelRegressor::computeDerivativeCorrMatrix(int dth_index)
  {
    const size_t nSamples = mData.getNSamples();
    matrixd corrMatrix(nSamples,nSamples);
    mKernel.computeDerivativeCorrMatrix(mData.mX,corrMatrix,dth_index);
    return corrMatrix;
  }

} //namespace bayesopt
