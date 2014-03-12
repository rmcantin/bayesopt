
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
//#include "optimizekernel.hpp"	


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

    size_t nhp = mKernel.nHyperParameters();
    kOptimizer = new NLOPT_Optimization(this,nhp);

    //TODO: Generalize
    if (parameters.l_type == L_ML)
      {
	kOptimizer->setAlgorithm(BOBYQA);    // local search to avoid underfitting
      }
    else
      {
	kOptimizer->setAlgorithm(COMBINED);
      }
    kOptimizer->setLimits(svectord(nhp,1e-10),svectord(nhp,100.));
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
    size_t n = mData.getNSamples();
    size_t last = n-1;
    int error = 0;
    double sum = 0.0;

    matrixd tempF(mMean.mFeatM);

    // Data point used for cross validation
    vectord x(n);
    double y;

    // We take the first element, use it for validation and then paste
    // it at the end. Thus, after every iteration, the first element
    // is different and, at the end, all the elements should have
    // rotated.
    for(size_t i = 0; i<n; ++i)
      {
	// Take the first element
	y = mData.getSample(0,x);

	// Remove it for cross validation
	mData.mX.erase(mData.mX.begin()); 
	utils::erase(mData.mY,mData.mY.begin());
	utils::erase_column(mMean.mFeatM,0);

	// Compute the cross validation
	precomputeSurrogate();
	ProbabilityDistribution* pd = prediction(x);
	sum += log(pd->pdf(y));

	//Paste it back at the end
	mData.addSample(x,y);
	mMean.mFeatM.resize(mMean.mFeatM.size1(),mMean.mFeatM.size2()+1);  
	mMean.mFeatM = tempF;
      }
    std::cout << "End" << mData.getNSamples();
    return -sum;   //Because we are minimizing.
  }

} // namespace bayesopt
