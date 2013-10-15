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

#include <cstdlib>
#include "log.hpp"
#include "bayesoptbase.hpp"

namespace bayesopt
{
  
  BayesOptBase::BayesOptBase():
    mGP(NULL), mCrit(NULL), mEngine()
  {
    mParameters = initialize_parameters_to_default();
    __init__();
  }


  BayesOptBase::BayesOptBase(size_t dim, bopt_params parameters):
    mGP(NULL), mCrit(NULL),  mParameters(parameters), mDims(dim)
  {
    __init__();
  }

  int BayesOptBase::__init__()
  { 
    // Configure logging
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
    if ((mParameters.n_iterations <= 0) || 
	(mParameters.n_iterations > MAX_ITERATIONS))
      mParameters.n_iterations = MAX_ITERATIONS;

    if (mParameters.n_init_samples <= 0)
      mParameters.n_init_samples = 
	static_cast<size_t>(ceil(0.1*mParameters.n_iterations));

    // Configure Surrogate and Criteria Functions
    setSurrogateModel();
    setCriteria();

    return 0;
  } // __init__

  BayesOptBase::~BayesOptBase()
  {} // Default destructor

  int BayesOptBase::setSurrogateModel()
  {
    // Configure Surrogate and Criteria Functions
    mGP.reset(NonParametricProcess::create(mDims,mParameters));
    if (mGP == NULL) 
      {
	FILE_LOG(logERROR) << "Error setting the surrogate function"; 
	exit(EXIT_FAILURE);
      } 
    return 0;
  } // setSurrogateModel

  int BayesOptBase::setCriteria()
  {
    mCrit.reset(mCFactory.create(mParameters.crit_name,mGP.get()));
    if (mCrit == NULL)
      {
	FILE_LOG(logERROR) << "Error in criterium"; 
	exit(EXIT_FAILURE);
      }
    
    if (mCrit->nParameters() != mParameters.n_crit_params)
      {
	if (mParameters.n_crit_params != 0)
	  {
	    FILE_LOG(logERROR) << "Expected " << mCrit->nParameters() 
			       << " parameters. Got " 
			       << mParameters.n_crit_params << " instead.";
	  }
	FILE_LOG(logINFO) << "Usign default parameters for criteria.";
	return 0;
      }
      
    // If we have the correct number of parameters.
    vectord critParams = utils::array2vector(mParameters.crit_params,
					       mParameters.n_crit_params);
    mCrit->setParameters(critParams);
    return 0;
  } // setCriteria

  int BayesOptBase::stepOptimization(size_t ii)
  {
    vectord xNext(mDims);
    nextPoint(xNext); // Find what is the next point.
    
    double yNext = evaluateSampleInternal(xNext);

    // Update surrogate model
    if ((mParameters.n_iter_relearn > 0) && 
	((ii + 1) % mParameters.n_iter_relearn == 0))
      mGP->fullUpdateSurrogateModel(xNext,yNext); 
    else
      mGP->updateSurrogateModel(xNext,yNext); 
    
    plotStepData(ii,xNext,yNext);
    return 0;
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
  

  int BayesOptBase::nextPoint(vectord &Xnext)
  {
    int error = 0;
    
    //Epsilon-Greedy exploration (see Bull 2011)
    if ((mParameters.epsilon > 0.0) && (mParameters.epsilon < 1.0))
      {
	randFloat drawSample(mEngine,realUniformDist(0,1));
	double result = drawSample();
	FILE_LOG(logINFO) << "Trying random jump with prob:" << result;
	if (mParameters.epsilon > result)
	  {
	    for (size_t i = 0; i <Xnext.size(); ++i)
	      {
		 Xnext(i) = drawSample();
	      } 
	    FILE_LOG(logINFO) << "Epsilon-greedy random query!";
	    return 0;
	  }
      }

    if (mCrit->requireComparison())
      {
	bool check = false;
	std::string name;

	mCrit->reset();
	while (!check)
	  {
	    findOptimal(Xnext);
	    check = mCrit->checkIfBest(Xnext,name,error);
	  }
	FILE_LOG(logINFO) << name << " was selected.";
      }
    else
      {
	findOptimal(Xnext);
      }
    return error;
  }


} //namespace bayesopt

