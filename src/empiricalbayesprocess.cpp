
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


#include "empiricalbayesprocess.hpp"
#include "log.hpp"
#include "optimizekernel.hpp"	


namespace bayesopt
{
  EmpiricalBayesProcess::EmpiricalBayesProcess(size_t dim, bopt_params parameters):
    NonParametricProcess(dim,parameters)
  { 
    if (mLearnType == L_BAYES)
      {
	FILE_LOG(logERROR) << "Empirical Bayes model and full Bayes learning are incompatible.";
	throw 1;
      }

    kOptimizer = new OptimizeKernel(this);

    //TODO: Generalize
    if (parameters.l_type == L_ML)
      {
	kOptimizer->setAlgorithm(BOBYQA);    // local search to avoid underfitting
      }
    else
      {
	kOptimizer->setAlgorithm(COMBINED);
      }
    kOptimizer->setLimits(1e-10,100.);
  }

  EmpiricalBayesProcess::~EmpiricalBayesProcess()
  {
    delete kOptimizer;
  }


  int EmpiricalBayesProcess::updateKernelParameters()
  {
    if (mLearnType == L_FIXED)
      {
	FILE_LOG(logDEBUG) << "Fixed hyperparameters. Not learning";
	return 0;
      }
    else
      {
	int error = -1;
	vectord optimalTheta = mKernel.getHyperParameters();
	
	FILE_LOG(logDEBUG) << "Computing kernel parameters. Initial: " 
			   << optimalTheta;

	kOptimizer->run(optimalTheta);
	error = mKernel.setHyperParameters(optimalTheta);

	if (error)
	  {
	    FILE_LOG(logERROR) << "Error updating kernel parameters.";
	    exit(EXIT_FAILURE);
	  }   

	FILE_LOG(logDEBUG) << "Final kernel parameters: " << optimalTheta;	
	return error;
      }
  };

  double EmpiricalBayesProcess::evaluateKernelParams(const vectord& query)
  { 
    if (mLearnType == L_FIXED)
      {
	FILE_LOG(logERROR) << "Fixed hyperparameters should not change.";
	return -1;
      }
    else
      {
	int error = mKernel.setHyperParameters(query);
	if (error) 
	  {
	    FILE_LOG(logERROR) << "Problem optimizing kernel parameters."; 
	    exit(EXIT_FAILURE);	
	  }
	return evaluateKernelParams();
      }
  };
  
  double EmpiricalBayesProcess::evaluateKernelParams()
  { 
    double result;
    switch(mLearnType)
      {
      case L_ML:
	result = negativeTotalLogLikelihood(); break;
      case L_FIXED:
	result = negativeLogLikelihood(); break;
      case L_MAP:
	// It is a minus because the prior is the positive and we want
	// the negative.
	result = negativeLogLikelihood()-mKernel.kernelLogPrior(); break;
      case L_LOO:
	result = negativeCrossValidation(); break;
      default:
	FILE_LOG(logERROR) << "Learning type not supported";
      }	  
    return result;
  }


  double EmpiricalBayesProcess::negativeCrossValidation()
  {
    // This is highly ineffient implementation for comparison purposes.
    size_t n = mGPXX.size();
    size_t last = n-1;
    int error = 0;
    double sum = 0.0;
    vecOfvec tempXX(mGPXX);
    vectord tempY(mGPY);
    vectord tempM(mMeanV);
    matrixd tempF(mFeatM);
    for(size_t i = 0; i<n; ++i)
      {
	vectord x = mGPXX[0];  double y = mGPY(0);
	double m = mMeanV(0);

	mGPXX.erase(mGPXX.begin()); 
	utils::erase(mGPY,mGPY.begin());
	utils::erase(mMeanV,mMeanV.begin());
	utils::erase_column(mFeatM,0);

	precomputeSurrogate();
	ProbabilityDistribution* pd = prediction(x);
	sum += log(pd->pdf(y));
	mGPXX.push_back(x);     
	mGPY.resize(mGPY.size()+1);  mGPY(mGPY.size()-1) = y;
	mMeanV.resize(mGPY.size());  mMeanV(mGPY.size()-1) = m;
	mFeatM.resize(mFeatM.size1(),mFeatM.size2()+1);  
	mFeatM = tempF;
      }
      std::cout << "End" << mGPY.size();
    return -sum;
  }

} // namespace bayesopt
