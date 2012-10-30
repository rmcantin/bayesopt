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
  FILELog::ReportingLevel() = logINFO; //mParameters.verbose_level;
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
{
  if (mCrit != NULL)
    delete mCrit;

} // Default destructor

int BayesOptBase::setCriteriumFunction()
{
  if (mCrit != NULL)
    delete mCrit;

  // We need a valid surrogate function pointer for the criteria
  if(mGP == NULL)
    setSurrogateFunction();

  switch(mParameters.c_name)
    {
    case C_GP_HEDGE:        mCrit = new GP_Hedge(mGP.get()); break;
    case C_GP_HEDGE_RANDOM: mCrit = new GP_Hedge_Random(mGP.get()); break;
    case C_EI:
    case C_EI_A:
    case C_LCB: 
    case C_LCB_A:
    case C_POI:
    case C_GREEDY_A_OPTIMALITY: 
    case C_EXPECTED_RETURN:   
    case C_OPTIMISTIC_SAMPLING:
      mCrit = new SingleCriteria(mGP.get());
      mCrit->setCriterium(mParameters.c_name);
      break;
    case C_ERROR:
    default:
      FILE_LOG(logERROR) << "Error in criterium"; return -1;
    }

  if (mCrit == NULL)       
    {
      FILE_LOG(logERROR) << "Error in criterium"; return -1;
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



