
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
  EmpiricalBayesProcess::EmpiricalBayesProcess(size_t dim, bopt_params parameters, 
					       const Dataset& data):
    KernelRegressor(dim,parameters,data)
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


  void EmpiricalBayesProcess::updateKernelParameters()
  {
    if (mLearnType == L_FIXED)
      {
	FILE_LOG(logDEBUG) << "Fixed hyperparameters. Not learning";
      }
    else
      {
	vectord optimalTheta = mKernel.getHyperParameters();
	
	FILE_LOG(logDEBUG) << "Initial kernel parameters: " << optimalTheta;
	kOptimizer->run(optimalTheta);
	mKernel.setHyperParameters(optimalTheta);
	FILE_LOG(logDEBUG) << "Final kernel parameters: " << optimalTheta;	
      }
  };

  double EmpiricalBayesProcess::evaluateKernelParams()
  { 
    switch(mLearnType)
      {
      case L_ML:
	return negativeTotalLogLikelihood();
      case L_FIXED:
	return negativeLogLikelihood();
      case L_MAP:
	// It is a minus because the prior is the positive and we want
	// the negative.
	return negativeLogLikelihood()-mKernel.kernelLogPrior();
      case L_LOO:
	return negativeCrossValidation(); 
      default:
	FILE_LOG(logERROR) << "Learning type not supported";
	throw std::invalid_argument("Learning type not supported");
      }	  
  }


  double EmpiricalBayesProcess::negativeCrossValidation()
  {
    // This is highly ineffient implementation for comparison purposes.
    Dataset data(mData);

    size_t n = data.getNSamples();
    size_t last = n-1;
    int error = 0;
    double sum = 0.0;

    matrixd tempF(mMean.mFeatM);


    // We take the first element, use it for validation and then paste
    // it at the end. Thus, after every iteration, the first element
    // is different and, at the end, all the elements should have
    // rotated.
    for(size_t i = 0; i<n; ++i)
      {
	// Take the first element
	const double y = data.getSampleY(0);
	const vectord x = data.getSampleX(0);

	// Remove it for cross validation
	data.mX.erase(data.mX.begin()); 
	utils::erase(data.mY,data.mY.begin());
	utils::erase_column(mMean.mFeatM,0);

	// Compute the cross validation
	computeCholeskyCorrelation();
	precomputePrediction(); 
	ProbabilityDistribution* pd = prediction(x);
	sum += log(pd->pdf(y));

	//Paste it back at the end
	data.addSample(x,y);
	mMean.mFeatM.resize(mMean.mFeatM.size1(),mMean.mFeatM.size2()+1);  
	mMean.mFeatM = tempF;
      }
    std::cout << "End" << data.getNSamples();
    return -sum;   //Because we are minimizing.
  }

} // namespace bayesopt
