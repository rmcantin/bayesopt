/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2012 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
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
#include "gaussian_process_ign.hpp"
#include "gaussian_process.hpp"
#include "student_t_process.hpp"


BayesOptBase::BayesOptBase( bopt_params parameters,
			    bool uselogfile,
			    const char* logfilename):

  mGP(NULL), mCrit(NULL),  mParameters(parameters)
{ 
  switch(mParameters.verbose_level)
    {
    case 0: FILELog::ReportingLevel() = logWARNING; break;
    case 1: FILELog::ReportingLevel() = logINFO; break;
    case 2: FILELog::ReportingLevel() = logDEBUG4; break;
    default:
      FILELog::ReportingLevel() = logERROR; break;
    }

  
  if (uselogfile)
    {
      FILE* log_fd = fopen( logfilename , "w" );
      Output2FILE::Stream() = log_fd; 
    }

  setInitSet();
  setNumberIterations();
  setSurrogateFunction();
  setCriteriumFunction();
} // Constructor

BayesOptBase::~BayesOptBase()
{} // Default destructor

int BayesOptBase::setCriteriumFunction()
{
  // We need a valid surrogate function pointer for the criteria
  if(mGP == NULL)
    setSurrogateFunction();

  mCrit.reset(MetaCriteria::create(mParameters.c_name,mGP.get()));
  if (mCrit == NULL)       
    {
      FILE_LOG(logERROR) << "Error in criterium"; 
      return -1;
    }

  return 1;
}

int BayesOptBase::setSurrogateFunction()
{
  mGP.reset(NonParametricProcess::create(mParameters));
  if (mGP == NULL)  return -1;
  else              return 0;
}

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



