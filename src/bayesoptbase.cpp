/*
-----------------------------------------------------------------------------
   This file is part of BayesOptimization, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOptimization is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOptimization is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with BayesOptimization.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/
#include "bayesoptbase.hpp"
#include "gaussprocess.hpp"
#include "basicgaussprocess.hpp"
#include "studenttprocess.hpp"


SKO_BASE::SKO_BASE( sko_params parameters,
       bool uselogfile,
       const char* logfilename):
  Logger(uselogfile,logfilename),
  mGP(NULL)
{ 
  mParameters = parameters;
  setNumberIterations();
  setSurrogateFunction();
} // Constructor

SKO_BASE::~SKO_BASE()
{
  if (mGP != NULL)
    delete mGP;
} // Default destructor

int SKO_BASE::setSurrogateFunction()
{
  if (mGP != NULL)
    delete mGP;
 
  switch(mParameters.s_name)
    {
    case S_GAUSSIAN_PROCESS: 
      mGP = new BasicGaussianProcess(mParameters.noise); break;

    case S_GAUSSIAN_PROCESS_INV_GAMMA_NORMAL:
      mGP = new GaussianProcess(mParameters.noise, mParameters.alpha,
				mParameters.beta,mParameters.delta);  break;

    case S_STUDENT_T_PROCESS_JEFFREYS:
      mGP = new StudentTProcess(mParameters.noise); break;

    default:
      std::cout << "Error: surrogate function not supported." << std::endl;
      return -1;
    }
  
  mGP->setKernel(mParameters.theta,mParameters.k_name);
  return 0;
}



int SKO_BASE::nextPoint(vectord &Xnext)
{
  crit.resetAnnealValues();
  if (mParameters.c_name == C_GP_HEDGE)
    {
      vectord best_ei(Xnext);
      vectord best_lcb(Xnext);
      vectord best_poi(Xnext);
      double l_ei,l_lcb,l_poi,foo;

      crit.setCriterium(C_EI);
      findOptimal(best_ei);
      mGP->prediction(best_ei,l_ei,foo);
      
      crit.setCriterium(C_LCB);
      findOptimal(best_lcb);
      mGP->prediction(best_lcb,l_lcb,foo);

      crit.setCriterium(C_POI);
      findOptimal(best_poi);
      mGP->prediction(best_poi,l_poi,foo);

      // Since we want to find the minimum, the predicted value is loss value, not a
      // reward value.
      criterium_name better = crit.update_hedge(l_ei,l_lcb,l_poi);
      switch(better)
	{
	case C_EI: 
	  Xnext = best_ei;
	  if (mParameters.verbose_level > 0) mOutput << "EI used." << std::endl;
	  break;
	case C_LCB: 
	  Xnext = best_lcb; 
	  if (mParameters.verbose_level > 0) mOutput << "LCB used." << std::endl;
	  break;
	case C_POI: 
	  Xnext = best_poi; 
	  if (mParameters.verbose_level > 0) mOutput << "POI used." << std::endl;
	  break;
	default: return -1;
	}
      return 1;
    }
  else
    {
      crit.setCriterium(mParameters.c_name);
      return findOptimal(Xnext);
    }
}



