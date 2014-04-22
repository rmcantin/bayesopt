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

#include "bayesoptbase.hpp"

#include <ctime>
#include <cstdlib>
#include "log.hpp"
#include "posteriormodel.hpp"

namespace bayesopt
{


  BayesOptBase::BayesOptBase(size_t dim, bopt_params parameters):
    mParameters(parameters), mDims(dim)
  {
    if (mParameters.use_random_seed) mEngine.seed(std::time(0));

    mModel.reset(PosteriorModel::create(dim,parameters,mEngine));
    
    size_t verbose = mParameters.verbose_level;
    if (verbose>=3)
      {
	FILE* log_fd = fopen( mParameters.log_filename , "w" );
	Output2FILE::Stream() = log_fd; 
	verbose -= 3;
      }

    switch(verbose)
      {
      case 0: FILELog::ReportingLevel() = logWARNING; break;
      case 1: FILELog::ReportingLevel() = logINFO; break;
      case 2: FILELog::ReportingLevel() = logDEBUG4; break;
      default:
	FILELog::ReportingLevel() = logERROR; break;
      }

    // Configure iteration parameters
    if (mParameters.n_init_samples <= 0)
      mParameters.n_init_samples = 
	static_cast<size_t>(ceil(0.1*mParameters.n_iterations));

  }

  BayesOptBase::~BayesOptBase()
  { } // Default destructor


  void BayesOptBase::stepOptimization(size_t ii)
  {
    // Find what is the next point.
    const vectord xNext = nextPoint(); 
    const double yNext = evaluateSampleInternal(xNext);

    if (yNext == HUGE_VAL)
      {
	throw std::runtime_error("Function evaluation out of range");
      }

    mModel->addSample(xNext,yNext);

    // Update surrogate model
    bool retrain = ((mParameters.n_iter_relearn > 0) && 
		    ((ii + 1) % mParameters.n_iter_relearn == 0));
    
    if (retrain)  // Full update
      {
	mModel->updateHyperParameters();
	mModel->fitSurrogateModel();
      }
    else          // Incremental update
      {
	mModel->updateSurrogateModel();
      } 
    plotStepData(ii,xNext,yNext);
  }

  void BayesOptBase::initializeOptimization()
  {
    size_t nSamples = mParameters.n_init_samples;

    matrixd xPoints(nSamples,mDims);
    vectord yPoints(nSamples);

    sampleInitialPoints(xPoints,yPoints);
    mModel->setSamples(xPoints,yPoints);
 
    if(mParameters.verbose_level > 0)
      {
	mModel->plotDataset(logDEBUG);
      }
    
    mModel->updateHyperParameters();
    mModel->fitSurrogateModel();
  }


  int BayesOptBase::optimize(vectord &bestPoint)
  {
    initializeOptimization();
    assert(mDims == bestPoint.size());
    
    for (size_t ii = 0; ii < mParameters.n_iterations; ++ii)
      {      
	stepOptimization(ii);
      }
   
    bestPoint = getFinalResult();

    return 0;
  } // optimize
  

  vectord BayesOptBase::nextPoint()
  {

    //Epsilon-Greedy exploration (see Bull 2011)
    if ((mParameters.epsilon > 0.0) && (mParameters.epsilon < 1.0))
      {
	randFloat drawSample(mEngine,realUniformDist(0,1));
	double result = drawSample();
	FILE_LOG(logINFO) << "Trying random jump with prob:" << result;
	if (mParameters.epsilon > result)
	  {
	    FILE_LOG(logINFO) << "Epsilon-greedy random query!";
	    return samplePoint();
	  }
      }

    //TODO: Try to solve this without bringing the pointer to the
    //criteria. Right now it does not work with MCMC.
    vectord Xnext(mDims);    

    // GP-Hedge and related algorithms
    if (mModel->criteriaRequiresComparison())
      {
	bool changed = true;

	mModel->setFirstCriterium();
	while (changed)
	  {
	    findOptimal(Xnext);
	    changed = mModel->setNextCriterium(Xnext);
	  }
	std::string name = mModel->getBestCriteria(Xnext);
	FILE_LOG(logINFO) << name << " was selected.";
      }
    else  // Standard "Bayesian optimization"
      {
	FILE_LOG(logDEBUG) << "------ Optimizing criteria ------";
	findOptimal(Xnext);
      }
    return Xnext;
  }


  // Potential inline functions. Moved here to simplify API and header
  // structure.
  double BayesOptBase::evaluateCriteria(const vectord& query)
  {
    if (checkReachability(query)) return mModel->evaluateCriteria(query);
    else return 0.0;
  }

  vectord BayesOptBase::getPointAtMinimum() 
  { return mModel->getPointAtMinimum(); };
  
  double BayesOptBase::getValueAtMinimum()
  { return mModel->getValueAtMinimum(); };

  ProbabilityDistribution* BayesOptBase::getPrediction(const vectord& query)
  { return mModel->getPrediction(query); };

   const Dataset* BayesOptBase::getData()
  { return mModel->getData(); };

  bopt_params* BayesOptBase::getParameters() 
  {return &mParameters;};


} //namespace bayesopt

