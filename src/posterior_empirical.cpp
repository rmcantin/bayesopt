/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2014 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
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
#include "inneroptimization.hpp"
#include "posterior_empirical.hpp"

namespace bayesopt
{


  EmpiricalBayes::EmpiricalBayes(size_t dim, bopt_params parameters, 
				 randEngine& eng):
    PosteriorModel(dim,parameters,eng)
  {
    // Configure Surrogate and Criteria Functions
    setSurrogateModel(eng);
    setCriteria(eng);

    // Seting kernel optimization
    size_t nhp = mGP->nHyperParameters();
    kOptimizer.reset(new NLOPT_Optimization(mGP.get(),nhp));

    //TODO: Generalize
    if (mParameters.sc_type == SC_ML)
      {
	// local search to avoid underfitting
	kOptimizer->setAlgorithm(BOBYQA); 
	kOptimizer->setMaxEvals(20*nhp);
      }
    else
      {
	kOptimizer->setAlgorithm(COMBINED);
	kOptimizer->setMaxEvals(50*nhp);
      }
    //Limits in log space
    kOptimizer->setLimits(svectord(nhp,-6.0),svectord(nhp,1.0));
  } 

  EmpiricalBayes::~EmpiricalBayes()
  { } // Default destructor


  void EmpiricalBayes::updateHyperParameters()
  {
    FILE_LOG(logDEBUG) << "------ Optimizing hyperparameters ------";
    vectord optimalTheta = mGP->getHyperParameters();

    FILE_LOG(logDEBUG) << "Initial hyper parameters: " << optimalTheta;
    kOptimizer->run(optimalTheta);
    mGP->setHyperParameters(optimalTheta);
    FILE_LOG(logDEBUG) << "Final hyper parameters: " << optimalTheta;	
  };


  void EmpiricalBayes::setSurrogateModel(randEngine& eng)
  {
    mGP.reset(NonParametricProcess::create(mDims,mParameters,
					   mData,mMean,eng));
  } // setSurrogateModel

  void EmpiricalBayes::setCriteria(randEngine& eng)
  {
    CriteriaFactory mCFactory;

    mCrit.reset(mCFactory.create(mParameters.crit_name,mGP.get()));
    mCrit->setRandomEngine(eng);

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


} //namespace bayesopt

