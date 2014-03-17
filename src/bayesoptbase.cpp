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

  void BayesOptBase::__init__()
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

    // mData.reset(new Dataset());
    // mGP->setData(mData.get());
  } // __init__

  BayesOptBase::~BayesOptBase()
  {} // Default destructor

  void BayesOptBase::setSurrogateModel()
  {
    mGP.reset(NonParametricProcess::create(mDims,mParameters,mData));
  } // setSurrogateModel

  void BayesOptBase::setCriteria()
  {
    mCrit.reset(mCFactory.create(mParameters.crit_name,mGP.get()));
    
    if (mCrit->nParameters() == mParameters.n_crit_params)
      {
	mCrit->setParameters(utils::array2vector(mParameters.crit_params,
					       mParameters.n_crit_params));
      }
    else // If the number of paramerters is different, use default.
      {
	if (mParameters.n_crit_params != 0)
	  {
	    FILE_LOG(logERROR) << "Expected " << mCrit->nParameters() 
			       << " parameters. Got " 
			       << mParameters.n_crit_params << " instead.";
	  }
	FILE_LOG(logINFO) << "Usign default parameters for criteria.";
      }
  } // setCriteria

  void BayesOptBase::stepOptimization(size_t ii)
  {
    vectord xNext(mDims);
    nextPoint(xNext); // Find what is the next point.
    
    double yNext = evaluateSampleInternal(xNext);

    // Update surrogate model
    bool retrain = ((mParameters.n_iter_relearn > 0) && 
		    ((ii + 1) % mParameters.n_iter_relearn == 0));
    
    addSample(xNext,yNext);
    mGP->updateSurrogateModel(xNext,yNext,retrain); 
    plotStepData(ii,xNext,yNext);
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
  

  void BayesOptBase::nextPoint(vectord &Xnext)
  {
    //Epsilon-Greedy exploration (see Bull 2011)
    if ((mParameters.epsilon > 0.0) && (mParameters.epsilon < 1.0))
      {
	randFloat drawSample(mEngine,realUniformDist(0,1));
	double result = drawSample();
	FILE_LOG(logINFO) << "Trying random jump with prob:" << result;
	if (mParameters.epsilon > result)
	  {
	    for (size_t i = 0; i<Xnext.size(); ++i)
	      {
		 Xnext(i) = drawSample();
	      } 
	    FILE_LOG(logINFO) << "Epsilon-greedy random query!";
	  }
      }
    else
      {
	// GP-Hedge and related algorithms
	if (mCrit->requireComparison())
	  {
	    bool check = false;
	    std::string name;
	    
	    mCrit->reset();
	    while (!check)
	      {
		findOptimal(Xnext);
		check = mCrit->checkIfBest(Xnext,name);
	      }
	    FILE_LOG(logINFO) << name << " was selected.";
	  }
	else  // Standard "Bayesian optimization"
	  {
	    findOptimal(Xnext);
	  }
      }
  }


} //namespace bayesopt

