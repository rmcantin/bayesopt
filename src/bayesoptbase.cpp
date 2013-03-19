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

#include "log.hpp"
#include "bayesoptbase.hpp"


BayesOptBase::BayesOptBase():
  mGP(NULL), mCrit(NULL)
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
  mGP.reset(NonParametricProcess::create(mDims,mParameters));
  if (mGP == NULL) 
    {
      FILE_LOG(logERROR) << "Error setting the surrogate function"; 
      return -1;
    } 

  mCrit.reset(MetaCriteria::create(mParameters.c_name,mGP.get()));
  if (mCrit == NULL)       
    {
      FILE_LOG(logERROR) << "Error in criterium"; 
      return -1;
    }

  return 0;
} // __init__

BayesOptBase::~BayesOptBase()
{} // Default destructor

int BayesOptBase::nextPoint(vectord &Xnext)
{
  bool check = false;
  criterium_name name;
  int error = 0;
  mCrit->initializeSearch();
  while (!check)
    {
      findOptimal(Xnext);
      check = mCrit->checkIfBest(Xnext,name,error);
    }

  if ((mParameters.c_name == C_GP_HEDGE) || 
      (mParameters.c_name == C_GP_HEDGE_RANDOM))
    {
      FILE_LOG(logINFO) << crit2str(name) << " was selected.";
    }

  return error;
}



